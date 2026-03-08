#!/usr/bin/env python3
"""general-purpose ASR + convergence pipeline for any Pfam family.

downloads pfam seed, builds diverse dataset, runs MAFFT + IQ-TREE with
ancestral reconstruction, computes per-residue NA contacts from PDB
structures, maps contacts to alignment columns, and generates convergence
summary showing ancestral states at specificity-determining positions.

usage:
  python asr_run_family.py PF02272 5O58 7BJQ
  python asr_run_family.py PF00575 7DID 9HVQ --dna-na-type DNA  # for split structures
"""

import argparse
import os
import sys
import re
import gzip
import json
import time
import subprocess
import tempfile
import requests
import numpy as np
from pathlib import Path
from collections import defaultdict
from io import StringIO

# ── constants ──
SCRIPT_DIR = Path(__file__).parent
PROJECT_DIR = SCRIPT_DIR.parent
CONDA_ENV = "phylo_asr"
AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")
AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "MSE": "M",
}
AA_RES = set(AA3TO1.keys())
NA_RES = {"A", "C", "G", "U", "DA", "DC", "DG", "DT", "DU", "I",
           "ADE", "CYT", "GUA", "URA", "THY", "PSU", "5MU", "H2U",
           "OMC", "OMG", "YYG", "1MA", "5MC", "7MG", "2MG", "M2G"}
RNA_RES = {"A", "C", "G", "U", "I", "ADE", "CYT", "GUA", "URA",
            "PSU", "5MU", "H2U", "OMC", "OMG", "YYG", "1MA", "5MC", "7MG", "2MG", "M2G"}
DNA_RES = {"DA", "DC", "DG", "DT", "DU", "THY"}


def log(msg):
    print(f"[asr_pipeline] {msg}", flush=True)


# ═══════════════════════════════════════════════════════════════════
# STEP 1: download and parse Pfam seed
# ═══════════════════════════════════════════════════════════════════

def download_pfam_seed(pfam_id, outdir):
    """download pfam seed alignment in stockholm format."""
    seed_path = outdir / f"{pfam_id}_seed.sto"
    if seed_path.exists() and seed_path.stat().st_size > 100:
        log(f"seed already exists: {seed_path}")
        return seed_path

    # try InterPro API (returns gzipped Stockholm)
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}/?annotation=alignment:seed"
    log(f"downloading seed from {url}")
    try:
        r = requests.get(url, timeout=120)
        if r.status_code == 200 and len(r.content) > 50:
            # always try gzip decompress first (InterPro returns gzipped)
            try:
                import gzip as gz
                text = gz.decompress(r.content).decode("utf-8", errors="replace")
            except Exception:
                text = r.text
            if "STOCKHOLM" in text or "#=GF" in text:
                seed_path.write_text(text)
                log(f"saved seed: {seed_path} ({len(text)} bytes)")
                return seed_path
            else:
                log(f"  response doesn't look like Stockholm format")
    except Exception as e:
        log(f"  failed: {e}")

    raise RuntimeError(f"could not download seed for {pfam_id}")


def parse_stockholm(sto_path):
    """parse stockholm alignment → list of (seq_id, accession, start, end, sequence).
    extracts domain boundaries from sequence headers like 'Q9YAX2_AERPE/7-71'."""
    sequences = []
    seq_data = {}  # id → sequence fragments
    ac_map = {}    # id → accession (from #=GS lines)

    with open(sto_path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#=GF") or line.startswith("#=GC") or line == "//":
                continue
            if line.startswith("#=GS"):
                # #=GS Q9YAX2_AERPE/7-71 AC Q9YAX2
                parts = line.split()
                if len(parts) >= 4 and parts[2] == "AC":
                    ac_map[parts[1]] = parts[3]
                continue
            if line.startswith("#"):
                continue

            # sequence line: "seq_id  sequence"
            parts = line.split()
            if len(parts) >= 2:
                sid = parts[0]
                seq_frag = parts[1]
                if sid not in seq_data:
                    seq_data[sid] = []
                seq_data[sid].append(seq_frag)

    for sid, frags in seq_data.items():
        full_seq = "".join(frags)
        # remove gap characters (., -, lowercase for inserts)
        clean_seq = re.sub(r'[.\-]', '', full_seq)
        # also remove lowercase (insert states in Stockholm)
        clean_seq_upper = re.sub(r'[a-z]', '', clean_seq)

        # parse boundaries from header: "ENTRY_SPECIES/start-end"
        acc = ac_map.get(sid, "")
        start, end = 0, 0
        m = re.search(r'/(\d+)-(\d+)$', sid)
        if m:
            start, end = int(m.group(1)), int(m.group(2))

        if len(clean_seq_upper) > 20:  # skip very short fragments
            sequences.append({
                "sid": sid,
                "accession": acc if acc else sid.split("/")[0].split("_")[0],
                "start": start,
                "end": end,
                "sequence": clean_seq_upper,
            })

    return sequences


# ═══════════════════════════════════════════════════════════════════
# STEP 2: curate dataset
# ═══════════════════════════════════════════════════════════════════

def curate_dataset(seed_seqs, pfam_id, outdir, cdhit_thresh=0.85, max_seqs=200):
    """select diverse sequences from seed. runs CD-HIT to reduce redundancy.
    if after CD-HIT the count exceeds max_seqs, runs a second round at lower threshold."""
    if len(seed_seqs) < 5:
        log(f"WARNING: only {len(seed_seqs)} seed sequences, using all")
        final_seqs = seed_seqs
    else:
        # write all seed sequences to FASTA
        raw_fasta = outdir / f"{pfam_id}_seed_raw.fasta"
        with open(raw_fasta, "w") as f:
            for s in seed_seqs:
                f.write(f">{s['sid']}|{s['accession']}|{s['start']}-{s['end']}\n{s['sequence']}\n")

        # run CD-HIT
        cdhit_out = outdir / f"{pfam_id}_cdhit.fasta"
        cmd = f"conda run -n {CONDA_ENV} cd-hit -i {raw_fasta} -o {cdhit_out} -c {cdhit_thresh} -n 5 -T 4 -M 2000"
        log(f"running CD-HIT at {cdhit_thresh}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            log(f"CD-HIT failed: {result.stderr[:200]}")
            # fallback: use all sequences
            final_seqs = seed_seqs
        else:
            # parse CD-HIT output — collect full headers of representative sequences
            cdhit_headers = set()
            with open(cdhit_out) as f:
                for line in f:
                    if line.startswith(">"):
                        cdhit_headers.add(line[1:].strip().split()[0])

            log(f"  CD-HIT kept {len(cdhit_headers)} representatives")

            # match back to original seed sequences by header
            final_seqs = []
            for s in seed_seqs:
                header = f"{s['sid']}|{s['accession']}|{s['start']}-{s['end']}"
                if header in cdhit_headers:
                    final_seqs.append(s)

            if len(final_seqs) < 5:
                # try partial matching (header might be truncated by CD-HIT)
                for s in seed_seqs:
                    header = f"{s['sid']}|{s['accession']}|{s['start']}-{s['end']}"
                    for ch in cdhit_headers:
                        if ch.startswith(header[:30]) or header.startswith(ch[:30]):
                            final_seqs.append(s)
                            break
                # deduplicate
                seen = set()
                deduped = []
                for s in final_seqs:
                    key = s['accession']
                    if key not in seen:
                        seen.add(key)
                        deduped.append(s)
                final_seqs = deduped

            if len(final_seqs) < 5:
                log(f"CD-HIT matching recovered only {len(final_seqs)} — using all seed instead")
                final_seqs = seed_seqs

    # if still too many, run CD-HIT at progressively lower thresholds
    if len(final_seqs) > max_seqs:
        for lower_thresh in [0.70, 0.60, 0.50, 0.40]:
            log(f"  {len(final_seqs)} > {max_seqs} — re-running CD-HIT at {lower_thresh}")
            tmp_fasta = outdir / f"{pfam_id}_cdhit_input.fasta"
            tmp_out = outdir / f"{pfam_id}_cdhit2.fasta"
            with open(tmp_fasta, "w") as f:
                for s in final_seqs:
                    f.write(f">{s['sid']}|{s['accession']}|{s['start']}-{s['end']}\n{s['sequence']}\n")
            # choose word size: n=5 for >0.7, n=4 for >0.6, n=3 for >0.5, n=2 for >0.4
            nword = 5 if lower_thresh > 0.7 else (4 if lower_thresh > 0.6 else (3 if lower_thresh > 0.5 else 2))
            cmd = f"conda run -n {CONDA_ENV} cd-hit -i {tmp_fasta} -o {tmp_out} -c {lower_thresh} -n {nword} -T 4 -M 2000"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                new_headers = set()
                with open(tmp_out) as f:
                    for line in f:
                        if line.startswith(">"):
                            new_headers.add(line[1:].strip().split()[0])
                new_seqs = [s for s in final_seqs
                            if f"{s['sid']}|{s['accession']}|{s['start']}-{s['end']}" in new_headers]
                log(f"  CD-HIT@{lower_thresh} reduced to {len(new_seqs)} sequences")
                final_seqs = new_seqs
                if len(final_seqs) <= max_seqs:
                    break

    # write final dataset
    final_fasta = outdir / f"{pfam_id}_sequences.fasta"
    with open(final_fasta, "w") as f:
        for s in final_seqs:
            header = f"{s['sid']}|{s['accession']}"
            f.write(f">{header}\n{s['sequence']}\n")

    log(f"curated dataset: {len(final_seqs)} sequences")
    return final_fasta, final_seqs


# ═══════════════════════════════════════════════════════════════════
# STEP 3: MAFFT alignment
# ═══════════════════════════════════════════════════════════════════

def run_mafft(fasta_path, outdir, pfam_id):
    """align sequences with MAFFT L-INS-i."""
    aln_path = outdir / f"{pfam_id}_aligned.fasta"
    if aln_path.exists() and aln_path.stat().st_size > 100:
        log(f"alignment exists: {aln_path}")
        return aln_path

    cmd = f"conda run -n {CONDA_ENV} mafft --maxiterate 1000 --localpair --thread 4 {fasta_path}"
    log("running MAFFT L-INS-i")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        # fallback to default MAFFT
        log(f"L-INS-i failed, trying default MAFFT: {result.stderr[:200]}")
        cmd = f"conda run -n {CONDA_ENV} mafft --thread 4 {fasta_path}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=300)
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr[:300]}")

    aln_path.write_text(result.stdout)
    # count sequences and true alignment length (not just first line width)
    n_seqs = result.stdout.count(">")
    # read first full sequence to get actual alignment length
    seq_chars = []
    in_seq = False
    for line in result.stdout.split("\n"):
        if line.startswith(">"):
            if in_seq:
                break  # finished first sequence
            in_seq = True
            continue
        if in_seq and line.strip():
            seq_chars.append(line.strip())
    aln_len = len("".join(seq_chars))
    log(f"alignment: {n_seqs} sequences, {aln_len} columns")
    return aln_path


# ═══════════════════════════════════════════════════════════════════
# STEP 4: IQ-TREE with ancestral reconstruction
# ═══════════════════════════════════════════════════════════════════

def run_iqtree(aln_path, outdir, pfam_id, threads=8):
    """run IQ-TREE with model selection, bootstrapping, and ASR."""
    prefix = outdir / f"{pfam_id}_asr"
    state_path = Path(f"{prefix}.state")
    tree_path = Path(f"{prefix}.treefile")

    if state_path.exists() and tree_path.exists():
        log(f"IQ-TREE results exist: {prefix}")
        return tree_path, state_path

    cmd = (f"conda run -n {CONDA_ENV} iqtree -s {aln_path} -m TEST "
           f"-bb 1000 --ancestral -nt {threads} -pre {prefix}")
    log(f"running IQ-TREE (threads={threads})")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=7200)
    if result.returncode != 0:
        # check if it's just a warning
        if state_path.exists():
            log("IQ-TREE completed with warnings")
        else:
            raise RuntimeError(f"IQ-TREE failed: {result.stderr[:500]}")

    # parse model from log
    log_path = Path(f"{prefix}.log")
    if log_path.exists():
        for line in open(log_path):
            if "Best-fit model:" in line:
                log(f"  {line.strip()}")
                break

    return tree_path, state_path


# ═══════════════════════════════════════════════════════════════════
# STEP 5: extract ancestral sequences
# ═══════════════════════════════════════════════════════════════════

def extract_ancestors(state_path, aln_path, outdir, pfam_id):
    """parse IQ-TREE .state file → ML ancestral sequences with PP."""
    anc_fasta = outdir / f"{pfam_id}_ancestors.fasta"
    anc_pp = outdir / f"{pfam_id}_ancestors_pp.tsv"

    nodes = defaultdict(dict)
    with open(state_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Node\t"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 23:
                continue
            node, site, state = parts[0], int(parts[1]), parts[2]
            probs = {}
            for i, aa in enumerate(AA_ORDER):
                probs[aa] = float(parts[3 + i])
            if state == "-":
                pp = 1.0
            else:
                pp = probs.get(state, 0.0)
            nodes[node][site] = (state, pp, probs)

    # build sequences and stats
    with open(anc_fasta, "w") as ff, open(anc_pp, "w") as fp:
        fp.write("node\tsite\tstate\tpp\n")
        for node in sorted(nodes.keys(), key=lambda x: int(x.replace("Node", ""))):
            sites = nodes[node]
            seq = ""
            pps = []
            for site_num in sorted(sites.keys()):
                state, pp, probs = sites[site_num]
                seq += state if state != "-" else "-"
                pps.append(pp)
                fp.write(f"{node}\t{site_num}\t{state}\t{pp:.4f}\n")

            # remove trailing gaps
            seq = seq.rstrip("-")
            mean_pp = np.mean(pps) if pps else 0
            n_high = sum(1 for p in pps if p >= 0.7)

            ff.write(f">{node}\n{seq}\n")

    n_nodes = len(nodes)
    if n_nodes > 0:
        root = f"Node1"
        if root in nodes:
            root_pps = [nodes[root][s][1] for s in sorted(nodes[root].keys())]
            root_mean = np.mean(root_pps)
            n_high = sum(1 for p in root_pps if p >= 0.7)
            log(f"extracted {n_nodes} ancestors. root mean PP={root_mean:.3f}, {n_high}/{len(root_pps)} sites PP>=0.7")

    return nodes


# ═══════════════════════════════════════════════════════════════════
# STEP 6: per-residue contacts
# ═══════════════════════════════════════════════════════════════════

def classify_na_atom(atom_name, resname):
    """classify nucleic acid atom into backbone/base/2oh/sugar categories."""
    a = atom_name.strip().upper()
    rn = resname.strip().upper()

    bb_atoms = {"P", "OP1", "OP2", "O5'", "C5'", "O3'", "C3'"}
    oh2_atoms = {"O2'"}
    sr_atoms = {"C1'", "C2'", "C4'", "O4'"}

    if a in bb_atoms:
        return "bb"
    elif a in oh2_atoms:
        return "2oh"
    elif a in sr_atoms:
        return "sr"
    else:
        return "base"


def download_pdb(pdb_id, outdir):
    """download PDB file."""
    pdb_path = outdir / f"{pdb_id.upper()}.pdb"
    if pdb_path.exists() and pdb_path.stat().st_size > 1000:
        return pdb_path

    pdb_id_lower = pdb_id.lower()
    urls = [
        f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
        f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb.gz",
    ]
    for url in urls:
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                if url.endswith(".gz"):
                    import gzip as gz
                    content = gz.decompress(r.content).decode()
                else:
                    content = r.text
                pdb_path.write_text(content)
                log(f"downloaded PDB: {pdb_id}")
                return pdb_path
        except Exception as e:
            continue

    # try mmCIF
    cif_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    try:
        r = requests.get(cif_url, timeout=60)
        if r.status_code == 200:
            cif_path = outdir / f"{pdb_id.upper()}.cif"
            cif_path.write_text(r.text)
            log(f"downloaded mmCIF: {pdb_id} (PDB format not available)")
            return cif_path
    except:
        pass

    raise RuntimeError(f"could not download structure for {pdb_id}")


def compute_contacts(pdb_id, na_type_filter, outdir, cutoff=4.0):
    """compute per-residue protein-NA contacts from a PDB structure.
    na_type_filter: 'RNA' or 'DNA' — only count contacts to this NA type."""
    from Bio.PDB import PDBParser, MMCIFParser
    from scipy.spatial import cKDTree

    pdb_path = download_pdb(pdb_id, outdir)

    if str(pdb_path).endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    s = parser.get_structure(pdb_id, str(pdb_path))
    m = s[0]

    # detect chains
    pc, nc = set(), set()
    for ch in m:
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn in AA_RES:
                pc.add(ch.id)
                break
            elif rn in NA_RES:
                nc.add(ch.id)
                break

    if not pc or not nc:
        log(f"WARNING: {pdb_id} has no protein({len(pc)}) or NA({len(nc)}) chains")
        return []

    # collect protein atoms
    prot_coords, prot_info = [], []
    for ch in m:
        if ch.id not in pc:
            continue
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in AA_RES:
                continue
            if r.get_id()[0] not in (' ', 'A'):
                continue
            for a in r:
                if a.element == 'H':
                    continue
                prot_coords.append(a.get_vector().get_array())
                prot_info.append((ch.id, r.get_id()[1], rn, a.get_name()))

    # collect NA atoms — filter by na_type_filter
    na_coords, na_info = [], []
    for ch in m:
        if ch.id not in nc:
            continue
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in NA_RES:
                continue
            is_rna = rn in RNA_RES
            is_dna = rn in DNA_RES
            na_tag = "RNA" if is_rna else "DNA"
            # filter: only include atoms of the requested NA type
            if na_type_filter == "RNA" and not is_rna:
                continue
            if na_type_filter == "DNA" and not is_dna:
                continue
            for a in r:
                if a.element == 'H':
                    continue
                cls = classify_na_atom(a.get_name(), rn)
                if cls is None:
                    continue
                na_coords.append(a.get_vector().get_array())
                na_info.append((ch.id, r.get_id()[1], rn, a.get_name(), cls, na_tag))

    if not prot_coords or not na_coords:
        log(f"WARNING: {pdb_id} no protein or filtered NA atoms")
        return []

    prot_arr = np.array(prot_coords, dtype=np.float32)
    na_arr = np.array(na_coords, dtype=np.float32)

    tree = cKDTree(na_arr)

    res_counts = defaultdict(lambda: {"bb": 0, "base": 0, "2oh": 0, "sr": 0})
    res_identity = {}

    for i, coord in enumerate(prot_arr):
        neighbors = tree.query_ball_point(coord, cutoff)
        if not neighbors:
            continue
        p_chain, p_resnum, p_resname, p_aname = prot_info[i]
        res_key = (p_chain, p_resnum)
        res_identity[res_key] = p_resname

        for j in neighbors:
            _, _, _, _, na_cls, _ = na_info[j]
            res_counts[res_key][na_cls] += 1

    # build contact list
    contacts = []
    for (chain, resnum), counts in sorted(res_counts.items()):
        resname = res_identity[(chain, resnum)]
        aa1 = AA3TO1.get(resname, "X")
        n_total = sum(counts.values())
        dominant = max(counts, key=counts.get)
        # classify
        has_2oh = counts["2oh"] > 0
        has_base = counts["base"] > 0
        has_bb = counts["bb"] > 0
        if has_2oh:
            classification = "2OH-contacting"
        elif has_base and has_bb:
            classification = "base+backbone"
        elif has_base:
            classification = "base-only"
        elif has_bb:
            classification = "backbone-only"
        else:
            classification = "sugar-ring-only"

        contacts.append({
            "chain": chain,
            "resnum": resnum,
            "resname": resname,
            "aa1": aa1,
            "n_bb": counts["bb"],
            "n_base": counts["base"],
            "n_2oh": counts["2oh"],
            "n_sr": counts["sr"],
            "n_total": n_total,
            "dominant": dominant,
            "classification": classification,
        })

    # save to TSV
    tag = na_type_filter.lower()
    tsv_path = outdir / f"{pdb_id}_{tag}_contacts.tsv"
    with open(tsv_path, "w") as f:
        if contacts:
            f.write("\t".join(contacts[0].keys()) + "\n")
            for c in contacts:
                f.write("\t".join(str(v) for v in c.values()) + "\n")
    log(f"  {pdb_id} {na_type_filter}: {len(contacts)} contacting residues")
    return contacts


# ═══════════════════════════════════════════════════════════════════
# STEP 7: domain mapping (PDB residue → alignment column)
# ═══════════════════════════════════════════════════════════════════

def get_sifts_mapping(pdb_id):
    """get PDB → UniProt mapping from SIFTS API.
    returns: {chain_id: (uniprot_id, [(pdb_start, pdb_end, unp_start, unp_end), ...])}"""
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            return {}
        data = r.json()
        result = {}
        pdb_data = data.get(pdb_id.lower(), {}).get("UniProt", {})
        for uniprot_id, info in pdb_data.items():
            for mapping in info.get("mappings", []):
                chain = mapping.get("chain_id", "A")
                seg = (mapping["start"]["residue_number"],
                       mapping["end"]["residue_number"],
                       mapping["unp_start"],
                       mapping["unp_end"])
                if chain not in result:
                    result[chain] = (uniprot_id, [])
                result[chain][1].append(seg)
        return result
    except Exception as e:
        log(f"  SIFTS failed for {pdb_id}: {e}")
        return {}


def pdb_resnum_to_uniprot(resnum, sifts_segments):
    """convert PDB residue number to UniProt position using SIFTS segments."""
    for pdb_start, pdb_end, unp_start, unp_end in sifts_segments:
        if pdb_start <= resnum <= pdb_end:
            offset = resnum - pdb_start
            return unp_start + offset
    return None


def find_domain_in_seed(uniprot_id, seed_seqs, pfam_id):
    """find domain boundaries for a UniProt ID in the Pfam seed.
    returns (domain_start, domain_end) in UniProt coordinates, or None."""
    # search by accession
    for s in seed_seqs:
        acc = s["accession"].split(".")[0]  # remove version
        if acc == uniprot_id and s["start"] > 0:
            return s["start"], s["end"]
    # search by ID prefix (e.g., Q9YAX2 in Q9YAX2_AERPE)
    for s in seed_seqs:
        if uniprot_id in s["sid"] and s["start"] > 0:
            return s["start"], s["end"]
    return None


def get_interpro_domain(pfam_id, uniprot_id):
    """get domain boundaries from InterPro API."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}/protein/uniprot/{uniprot_id}/"
    try:
        r = requests.get(url, timeout=30, headers={"Accept": "application/json"})
        if r.status_code != 200:
            return None
        data = r.json()
        # extract domain locations
        for result in data.get("results", []):
            for entry in result.get("entries", []):
                for loc in entry.get("entry_protein_locations", []):
                    for frag in loc.get("fragments", []):
                        return frag["start"], frag["end"]
    except Exception as e:
        log(f"  InterPro API failed: {e}")
    return None


def build_dompos_to_alncol(aligned_seq):
    """map ungapped domain position (1-based) → alignment column (1-based)."""
    mapping = {}
    dom_pos = 0
    for col_idx, char in enumerate(aligned_seq):
        if char != "-":
            dom_pos += 1
            mapping[dom_pos] = col_idx + 1
    return mapping


def extract_pdb_chain_sequence(pdb_id, outdir):
    """extract protein chain sequence from PDB with residue numbering.
    returns: list of (chain_id, resnum, aa1) for each protein residue."""
    from Bio.PDB import PDBParser, MMCIFParser

    pdb_path = outdir / f"{pdb_id.upper()}.pdb"
    cif_path = outdir / f"{pdb_id.upper()}.cif"

    if cif_path.exists():
        parser = MMCIFParser(QUIET=True)
        s = parser.get_structure(pdb_id, str(cif_path))
    elif pdb_path.exists():
        parser = PDBParser(QUIET=True)
        s = parser.get_structure(pdb_id, str(pdb_path))
    else:
        return []

    residues = []
    for ch in s[0]:
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn in AA_RES and r.get_id()[0] in (' ', 'A'):
                aa1 = AA3TO1.get(rn, "X")
                residues.append((ch.id, r.get_id()[1], aa1))
    return residues


def build_hmm_and_search(aln_path, pdb_chain_seq, pdb_id, outdir):
    """use HMMER to map PDB chain sequence to alignment columns.

    approach:
    1. hmmbuild: alignment → HMM profile
    2. hmmsearch: PDB chain → HMM hits
    3. parse: PDB_resnum → HMM_position → alignment_column

    returns: dict mapping PDB_resnum → alignment_column (1-based).
    """
    hmm_path = outdir / "model.hmm"
    chain_fasta = outdir / f"{pdb_id}_chain.fasta"

    # step 1: build HMM (reuse if exists)
    if not hmm_path.exists():
        cmd = f"conda run -n {CONDA_ENV} hmmbuild {hmm_path} {aln_path}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            log(f"  hmmbuild failed: {result.stderr[:200]}")
            return {}

    # parse HMM to get match_state → alignment_column mapping
    hmm_to_alncol = parse_hmm_map(hmm_path)

    # step 2: write ALL protein chains to FASTA and let hmmsearch find the right one
    chains = defaultdict(list)
    for chain_id, resnum, aa1 in pdb_chain_seq:
        chains[chain_id].append((resnum, aa1))

    if not chains:
        return {}

    # write all chains as separate FASTA entries
    with open(chain_fasta, "w") as f:
        for cid in sorted(chains):
            seq = "".join(aa for _, aa in chains[cid])
            f.write(f">{pdb_id}_chain{cid}\n{seq}\n")

    # step 3: run hmmsearch against all chains
    cmd = f"conda run -n {CONDA_ENV} hmmsearch {hmm_path} {chain_fasta}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        log(f"  hmmsearch failed: {result.stderr[:200]}")
        return {}

    # find best-scoring chain from hmmsearch output
    best_chain = None
    best_score = -float("inf")
    for line in result.stdout.split("\n"):
        parts = line.split()
        # score table line: "E-value score bias E-value score bias exp N Sequence"
        if len(parts) >= 9 and parts[-1].startswith(f"{pdb_id}_chain"):
            try:
                score = float(parts[1])
                chain_id = parts[-1].replace(f"{pdb_id}_chain", "")
                if score > best_score:
                    best_score = score
                    best_chain = chain_id
            except ValueError:
                pass

    if best_chain is None:
        log(f"  hmmsearch found no matching chain in {pdb_id}")
        return {}

    log(f"  hmmsearch best chain: {best_chain} (score={best_score:.1f})")
    chain_residues = chains[best_chain]

    # step 4: parse hmmsearch alignment to get PDB_pos → HMM_pos mapping
    pdb_to_hmm = parse_hmmsearch_alignment(result.stdout, pdb_id, best_chain)

    # step 5: chain: PDB_resnum → HMM_pos → alignment_column
    pdb_to_aln = {}
    for pdb_pos, hmm_pos in pdb_to_hmm.items():
        aln_col = hmm_to_alncol.get(hmm_pos)
        if aln_col:
            pdb_to_aln[pdb_pos] = aln_col

    # map 0-based sequence positions back to PDB residue numbers
    resnum_to_aln = {}
    for pdb_pos, aln_col in pdb_to_aln.items():
        if 0 <= pdb_pos < len(chain_residues):
            resnum = chain_residues[pdb_pos][0]
            resnum_to_aln[resnum] = aln_col

    return resnum_to_aln, best_chain


def parse_hmm_map(hmm_path):
    """parse HMM file to get match_state_number → alignment_column mapping.
    hmmbuild stores MAP annotations: each match state line has the alignment column."""
    mapping = {}
    in_model = False
    match_num = 0
    with open(hmm_path) as f:
        for line in f:
            if line.startswith("HMM "):
                in_model = True
                continue
            if not in_model:
                continue
            if line.startswith("//"):
                break
            parts = line.split()
            if not parts:
                continue
            # match state lines start with a number (the match state index)
            try:
                state_idx = int(parts[0])
                # the last field on match state lines is the MAP value (alignment column)
                # format: idx aa1 aa2 ... aa20 MAP RF CS
                # MAP is at position -3 (third from end) or -1 (last)
                # in HMMER3 format, MAP is the field after the 20 emission probs
                # typical: "1  2.34 ... 0.12  1  -  -"  where 1 is the MAP column
                map_val = int(parts[21])  # 20 emission log-probs + 1 for index, MAP at 22nd
                mapping[state_idx] = map_val
            except (ValueError, IndexError):
                pass
    return mapping


def parse_hmmsearch_alignment(output_text, pdb_id, chain_id):
    """parse hmmsearch text output to extract PDB_sequence_position → HMM_position mapping.

    hmmsearch alignment blocks have 4 lines each:
      1. model line:  "  ModelName  1 ngkvlvakle......seq 87"
      2. match line:  "              ngk+++++++      ..."
      3. seq line:    "  SeqName  215 NGKIAYSYIDydtylr... 302"
      4. PP line:     "              6799999999... PP"
    blocks are separated by blank lines within a domain section.

    only parses domains for the specified chain_id (matching ">>" header).
    returns: dict {0-based_seq_position: hmm_position}"""
    mapping = {}

    # target sequence header in hmmsearch output: ">> PDB_chainX"
    target_header = f"{pdb_id}_chain{chain_id}"

    lines = output_text.split("\n")
    i = 0
    in_target_chain = False
    while i < len(lines):
        stripped = lines[i].strip()

        # track which chain's domain section we're in
        if stripped.startswith(">>"):
            # extract sequence name: ">> 8CSQ_chain5" → "8CSQ_chain5"
            seq_name = stripped.replace(">>", "").strip().split()[0] if stripped.replace(">>", "").strip() else ""
            in_target_chain = (seq_name == target_header)
            i += 1
            continue

        # only parse domains for our target chain
        if "== domain" not in lines[i] or not in_target_chain:
            i += 1
            continue
        i += 1  # move past "== domain" line

        # process alignment blocks for this domain
        while i < len(lines):
            # skip blank lines between blocks
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i >= len(lines):
                break

            # check for end of domain alignment
            stripped = lines[i].strip()
            if stripped.startswith(">>") or "== domain" in stripped or \
               stripped.startswith("Internal") or stripped.startswith("//"):
                break

            # expect 4 lines per block: model, match, seq, PP
            if i + 3 >= len(lines):
                break

            model_parts = lines[i].split()
            seq_parts = lines[i + 2].split()

            if len(model_parts) >= 4 and len(seq_parts) >= 4:
                try:
                    hmm_start = int(model_parts[1])
                    hmm_seq = model_parts[2]
                    seq_start = int(seq_parts[1])
                    seq_seq = seq_parts[2]

                    # align positions character by character
                    # model: letter = match state (advance HMM), '.' = insert (don't advance)
                    # seq: letter = residue (advance seq), '-' = gap (don't advance)
                    h_pos = hmm_start
                    s_pos = seq_start - 1  # 0-based
                    for h_char, s_char in zip(hmm_seq, seq_seq):
                        if h_char != '.' and s_char != '-':
                            # both present: map seq position to HMM match state
                            mapping[s_pos] = h_pos
                            h_pos += 1
                            s_pos += 1
                        elif h_char == '.' and s_char != '-':
                            # insert in sequence (no HMM match state)
                            s_pos += 1
                        elif s_char == '-' and h_char != '.':
                            # deletion in sequence (HMM has match state, seq has gap)
                            h_pos += 1
                except ValueError:
                    pass

            i += 4  # advance past all 4 lines of this block

        # continue outer loop — i is already positioned at the end marker
        continue

    return mapping


def map_contacts_to_alignment(contacts, aln_seqs, seed_seqs, pfam_id, pdb_id, outdir):
    """map PDB contact residue numbers to alignment columns using HMMER.

    approach:
    1. extract PDB chain sequence
    2. hmmbuild from alignment → HMM profile
    3. hmmsearch PDB chain against HMM → PDB_resnum → HMM_pos → alignment_col
    4. for each contact, look up alignment column

    returns: list of (aln_col, contact_dict) tuples.
    """
    if not contacts:
        return []

    # extract PDB chain sequence
    pdb_residues = extract_pdb_chain_sequence(pdb_id, outdir)
    if not pdb_residues:
        log(f"  WARNING: could not extract chain sequence from {pdb_id}")
        return []

    # get alignment path
    aln_path = None
    for fname in outdir.iterdir():
        if fname.name.endswith("_aligned.fasta"):
            aln_path = fname
            break

    if aln_path is None:
        log(f"  WARNING: no alignment file found in {outdir}")
        return []

    # run HMMER mapping
    hmmer_result = build_hmm_and_search(aln_path, pdb_residues, pdb_id, outdir)

    if not hmmer_result or not hmmer_result[0]:
        log(f"  WARNING: HMMER mapping failed for {pdb_id}, trying SIFTS fallback")
        return _sifts_fallback(contacts, aln_seqs, seed_seqs, pfam_id, pdb_id, outdir)

    resnum_to_aln, matched_chain = hmmer_result

    # filter contacts to only those on the domain-containing chain
    # (in multi-chain complexes, contacts from other proteins shouldn't be mapped)
    chain_contacts = [c for c in contacts if c["chain"] == matched_chain]
    other_chain = len(contacts) - len(chain_contacts)
    if other_chain > 0:
        log(f"  filtered to {len(chain_contacts)} contacts on chain {matched_chain} "
            f"({other_chain} on other chains)")

    # map each contact
    mapped = []
    unmapped = 0
    mismatches = 0
    for c in chain_contacts:
        pdb_resnum = c["resnum"]
        aln_col = resnum_to_aln.get(pdb_resnum)
        if aln_col is None:
            unmapped += 1
            continue

        # verify amino acid match (find the alignment sequence and check)
        # just check first non-gap character at this column
        verified = False
        for sid, seq in aln_seqs.items():
            if aln_col - 1 < len(seq) and seq[aln_col - 1] == c["aa1"]:
                verified = True
                break

        mapped.append((aln_col, c))

    if unmapped > 0:
        log(f"  {unmapped} contacts could not be mapped (outside domain or gaps)")
    log(f"  {len(mapped)} contacts mapped to alignment columns")
    return mapped


def _sifts_fallback(contacts, aln_seqs, seed_seqs, pfam_id, pdb_id, outdir):
    """fallback: use SIFTS + seed boundaries for mapping."""
    sifts = get_sifts_mapping(pdb_id)
    if not sifts:
        return []

    contact_chains = set(c["chain"] for c in contacts)
    protein_chain = sorted(contact_chains)[0]

    uniprot_id = None
    sifts_segments = []
    if protein_chain in sifts:
        uniprot_id, sifts_segments = sifts[protein_chain]
    elif sifts:
        for ch, (uid, segs) in sifts.items():
            uniprot_id, sifts_segments = uid, segs
            break

    if not uniprot_id:
        return []

    bounds = find_domain_in_seed(uniprot_id, seed_seqs, pfam_id)
    if not bounds:
        bounds = get_interpro_domain(pfam_id, uniprot_id)
    if not bounds:
        return []

    domain_start, domain_end = bounds
    aln_seq = None
    for sid, seq in aln_seqs.items():
        if uniprot_id in sid:
            aln_seq = seq
            break
    if not aln_seq:
        return []

    dompos_to_col = build_dompos_to_alncol(aln_seq)

    mapped = []
    for c in contacts:
        uniprot_pos = pdb_resnum_to_uniprot(c["resnum"], sifts_segments) if sifts_segments else c["resnum"]
        if uniprot_pos is None:
            continue
        domain_pos = uniprot_pos - domain_start + 1
        if domain_pos < 1 or domain_pos > (domain_end - domain_start + 1):
            continue
        aln_col = dompos_to_col.get(domain_pos)
        if aln_col:
            mapped.append((aln_col, c))

    log(f"  SIFTS fallback: {len(mapped)} contacts mapped")
    return mapped


# ═══════════════════════════════════════════════════════════════════
# STEP 8: convergence analysis
# ═══════════════════════════════════════════════════════════════════

def parse_alignment(fasta_path):
    """parse FASTA alignment → {seq_id: aligned_sequence}."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(current_seq)
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_id:
        seqs[current_id] = "".join(current_seq)
    return seqs


def run_convergence(rna_mapped, dna_mapped, ancestor_nodes, aln_path, outdir, pfam_id):
    """cross-reference ancestral states with specificity-determining positions."""
    # categorize alignment columns
    rna_cols = {}  # col → contact_dict
    dna_cols = {}
    for col, c in rna_mapped:
        rna_cols[col] = c
    for col, c in dna_mapped:
        dna_cols[col] = c

    rna_only = set(rna_cols.keys()) - set(dna_cols.keys())
    dna_only = set(dna_cols.keys()) - set(rna_cols.keys())
    shared = set(rna_cols.keys()) & set(dna_cols.keys())
    all_cols = sorted(rna_only | dna_only | shared)

    if not all_cols:
        log("WARNING: no contacts mapped to alignment — skipping convergence")
        # write empty summary
        summary_path = outdir / f"{pfam_id}_convergence_summary.tsv"
        summary_path.write_text("# no contacts could be mapped to alignment columns\n")
        return

    log(f"\nCONVERGENCE: {len(rna_only)} RNA-only, {len(dna_only)} DNA-only, {len(shared)} shared columns")

    # find key ancestor nodes
    key_nodes = []
    for node_name in sorted(ancestor_nodes.keys(), key=lambda x: int(x.replace("Node", ""))):
        node_num = int(node_name.replace("Node", ""))
        if node_num <= 5:  # include root + early nodes
            key_nodes.append(node_name)

    # also add root explicitly
    if "Node1" not in [n for n in key_nodes]:
        if "Node1" in ancestor_nodes:
            key_nodes.insert(0, "Node1")

    # limit to first 5 key nodes
    key_nodes = key_nodes[:5]

    # header
    print(f"\n{'='*80}")
    print(f"CONVERGENCE ANALYSIS FOR {pfam_id}")
    print(f"{'='*80}")
    print(f"RNA-only columns: {sorted(rna_only)}")
    print(f"DNA-only columns: {sorted(dna_only)}")
    print(f"Shared columns:   {sorted(shared)}")

    # print detailed table
    print(f"\n{'col':>4} {'rna_aa':>6} {'dna_aa':>6} {'type':>10}", end="")
    for node in key_nodes:
        print(f" | {node:>10}", end="")
    print()
    print("-" * (30 + 14 * len(key_nodes)))

    for col in all_cols:
        rna_c = rna_cols.get(col)
        dna_c = dna_cols.get(col)
        rna_aa = rna_c["aa1"] if rna_c else "-"
        dna_aa = dna_c["aa1"] if dna_c else "-"

        if col in shared:
            col_type = "shared"
        elif col in rna_only:
            col_type = "RNA"
            if rna_c and "2OH" in rna_c.get("classification", ""):
                col_type += "*"
        else:
            col_type = "DNA"

        print(f"{col:>4} {rna_aa:>6} {dna_aa:>6} {col_type:>10}", end="")

        for node in key_nodes:
            if node in ancestor_nodes and col in ancestor_nodes[node]:
                state, pp, probs = ancestor_nodes[node][col]
                match = ""
                if state == rna_aa and rna_c:
                    match = "=R"
                elif state == dna_aa and dna_c:
                    match = "=D"
                print(f" | {state}({pp:.2f}){match:>3}", end="")
            else:
                print(f" | {'?':>10}", end="")
        print()

    # summary per node
    print(f"\n{'='*80}")
    print(f"SUMMARY: ANCESTRAL RESIDUE MATCHES")
    print(f"{'='*80}")

    summary_rows = []
    for node in key_nodes:
        if node not in ancestor_nodes:
            continue
        rna_match, dna_match, neither, total, high_pp = 0, 0, 0, 0, 0

        for col in rna_only:
            if col in ancestor_nodes[node]:
                state, pp, _ = ancestor_nodes[node][col]
                total += 1
                if pp >= 0.7: high_pp += 1
                if state == rna_cols[col]["aa1"]:
                    rna_match += 1
                else:
                    neither += 1

        for col in dna_only:
            if col in ancestor_nodes[node]:
                state, pp, _ = ancestor_nodes[node][col]
                total += 1
                if pp >= 0.7: high_pp += 1
                if state == dna_cols[col]["aa1"]:
                    dna_match += 1
                else:
                    neither += 1

        for col in shared:
            if col in ancestor_nodes[node]:
                state, pp, _ = ancestor_nodes[node][col]
                total += 1
                if pp >= 0.7: high_pp += 1
                if state == rna_cols[col]["aa1"]:
                    rna_match += 1
                elif state == dna_cols[col]["aa1"]:
                    dna_match += 1
                else:
                    neither += 1

        if total > 0:
            pct_rna = 100 * rna_match / total
            pct_dna = 100 * dna_match / total
            print(f"\n{node}:")
            print(f"  RNA matches: {rna_match}/{total} ({pct_rna:.0f}%)")
            print(f"  DNA matches: {dna_match}/{total} ({pct_dna:.0f}%)")
            print(f"  neither:     {neither}/{total}")
            print(f"  PP>=0.7:     {high_pp}/{total}")

            summary_rows.append({
                "node": node, "rna_match": rna_match, "dna_match": dna_match,
                "neither": neither, "total": total, "high_pp": high_pp,
                "pct_rna": pct_rna, "pct_dna": pct_dna,
            })

    # 2'OH residue analysis
    print(f"\n{'='*80}")
    print(f"2'OH RESIDUE ANALYSIS")
    print(f"{'='*80}")
    oh2_found = False
    for col, c in rna_mapped:
        if "2OH" in c.get("classification", ""):
            oh2_found = True
            print(f"\n2'OH contact: {c['aa1']} at PDB resnum {c['resnum']}, aln col {col}")
            for node in key_nodes:
                if node in ancestor_nodes and col in ancestor_nodes[node]:
                    state, pp, probs = ancestor_nodes[node][col]
                    p_modern = probs.get(c["aa1"], 0)
                    print(f"  {node}: {state} (PP={pp:.3f}), P({c['aa1']})={p_modern:.3f}")
    if not oh2_found:
        print("  no 2'OH-contacting residues found")

    # write convergence table
    summary_path = outdir / f"{pfam_id}_convergence_summary.tsv"
    with open(summary_path, "w") as f:
        header = ["aln_col", "rna_aa", "dna_aa", "col_type", "rna_classification", "dna_classification"]
        for node in key_nodes:
            header.extend([f"{node}_state", f"{node}_pp"])
        f.write("\t".join(header) + "\n")

        for col in all_cols:
            rna_c = rna_cols.get(col)
            dna_c = dna_cols.get(col)
            rna_aa = rna_c["aa1"] if rna_c else ""
            dna_aa = dna_c["aa1"] if dna_c else ""
            rna_cls = rna_c["classification"] if rna_c else ""
            dna_cls = dna_c["classification"] if dna_c else ""
            if col in shared: ct = "shared"
            elif col in rna_only: ct = "RNA-only"
            else: ct = "DNA-only"

            row = [str(col), rna_aa, dna_aa, ct, rna_cls, dna_cls]
            for node in key_nodes:
                if node in ancestor_nodes and col in ancestor_nodes[node]:
                    state, pp, _ = ancestor_nodes[node][col]
                    row.extend([state, f"{pp:.4f}"])
                else:
                    row.extend(["", ""])
            f.write("\t".join(row) + "\n")

    # write node summary
    node_summary_path = outdir / f"{pfam_id}_node_summary.tsv"
    with open(node_summary_path, "w") as f:
        f.write("node\trna_match\tdna_match\tneither\ttotal\thigh_pp\tpct_rna\tpct_dna\n")
        for row in summary_rows:
            f.write("\t".join(str(row[k]) for k in
                ["node", "rna_match", "dna_match", "neither", "total", "high_pp", "pct_rna", "pct_dna"]) + "\n")

    log(f"convergence summary saved to {summary_path}")
    return summary_rows


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="ASR + convergence pipeline for a Pfam family")
    parser.add_argument("pfam_id", help="Pfam ID (e.g., PF02272)")
    parser.add_argument("rna_pdb", help="PDB ID for RNA-bound structure")
    parser.add_argument("dna_pdb", help="PDB ID for DNA-bound structure")
    parser.add_argument("--cdhit", type=float, default=0.85, help="CD-HIT identity threshold")
    parser.add_argument("--max-seqs", type=int, default=200, help="max sequences for IQ-TREE (caps via iterative CD-HIT)")
    parser.add_argument("--threads", type=int, default=8, help="IQ-TREE threads")
    args = parser.parse_args()

    pfam_id = args.pfam_id.upper()
    outdir = PROJECT_DIR / "results" / "asr" / pfam_id.lower()
    outdir.mkdir(parents=True, exist_ok=True)

    log(f"starting pipeline for {pfam_id}")
    log(f"  RNA PDB: {args.rna_pdb}, DNA PDB: {args.dna_pdb}")
    log(f"  output: {outdir}")

    # step 1: download pfam seed
    seed_path = download_pfam_seed(pfam_id, outdir)

    # step 2: parse and curate
    seed_seqs = parse_stockholm(seed_path)
    log(f"parsed {len(seed_seqs)} seed sequences")
    fasta_path, final_seqs = curate_dataset(seed_seqs, pfam_id, outdir, args.cdhit, args.max_seqs)

    # step 3: MAFFT alignment
    aln_path = run_mafft(fasta_path, outdir, pfam_id)

    # step 4: IQ-TREE ASR
    tree_path, state_path = run_iqtree(aln_path, outdir, pfam_id, args.threads)

    # step 5: extract ancestors
    ancestor_nodes = extract_ancestors(state_path, aln_path, outdir, pfam_id)

    # step 6: per-residue contacts
    log(f"\ncomputing RNA contacts from {args.rna_pdb}")
    rna_contacts = compute_contacts(args.rna_pdb, "RNA", outdir)
    log(f"computing DNA contacts from {args.dna_pdb}")
    dna_contacts = compute_contacts(args.dna_pdb, "DNA", outdir)

    # step 7: map contacts to alignment
    aln_seqs = parse_alignment(aln_path)

    log(f"\nmapping RNA contacts to alignment")
    rna_mapped = map_contacts_to_alignment(rna_contacts, aln_seqs, seed_seqs,
                                            pfam_id, args.rna_pdb, outdir)
    log(f"mapping DNA contacts to alignment")
    dna_mapped = map_contacts_to_alignment(dna_contacts, aln_seqs, seed_seqs,
                                            pfam_id, args.dna_pdb, outdir)

    # step 8: convergence
    run_convergence(rna_mapped, dna_mapped, ancestor_nodes, aln_path, outdir, pfam_id)

    log(f"\npipeline complete for {pfam_id}")


if __name__ == "__main__":
    main()
