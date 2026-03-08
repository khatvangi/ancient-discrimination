#!/usr/bin/env python3
"""
NAS-Bench v0.2 — Nucleic Acid Specificity Index
Memory-efficient version using scipy KDTree.
"""

import sys, os, gc, warnings, tempfile, urllib.request
from collections import defaultdict
import numpy as np

warnings.filterwarnings("ignore")

# === Constants ===
BACKBONE_ATOMS = {"P","OP1","OP2","O5'","O3'","O5*","O3*"}
SUGAR_RING_ATOMS = {"C1'","C2'","C3'","C4'","C5'","O4'","C1*","C2*","C3*","C4*","C5*","O4*"}
SUGAR_2OH = {"O2'","O2*"}
RNA_RES = {"A","C","G","U","RA","RC","RG","RU","ADE","CYT","GUA","URA",
           "PSU","5MU","H2U","OMC","OMG","YG","M2G","2MG","7MG","5MC","1MA","MA6"}
DNA_RES = {"DA","DC","DG","DT"}
NA_RES = RNA_RES | DNA_RES
AA_RES = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
          "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","MSE","SEC"}

# 3-letter to 1-letter amino acid mapping (includes selenomethionine)
AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U",
}

def classify(aname, rname):
    n = aname.strip()
    if n[0] in "H123" and len(n) > 1: return None
    if n in BACKBONE_ATOMS: return "bb"
    if n in SUGAR_2OH: return "2oh"
    if n in SUGAR_RING_ATOMS: return "sr"
    if rname.strip().upper() in NA_RES: return "base"
    return None

def fetch(pid):
    d = os.path.join(tempfile.gettempdir(), 'nasbench')
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, f'{pid}.pdb')
    if not os.path.exists(p):
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pid}.pdb', p)
    return p

def compute_si(pdb_path, pdb_id="X", cutoff=4.0):
    from Bio.PDB import PDBParser
    from scipy.spatial import cKDTree

    s = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    m = s[0]

    # Detect chains
    pc, nc = set(), set()
    for ch in m:
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn in AA_RES: pc.add(ch.id); break
            elif rn in NA_RES: nc.add(ch.id); break

    if not pc or not nc:
        raise ValueError(f"No prot({pc}) or NA({nc}) chains")

    # Collect atoms as lightweight lists
    prot_coords, prot_info = [], []
    for ch in m:
        if ch.id not in pc: continue
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in AA_RES: continue
            if r.get_id()[0] not in (' ','A'): continue
            for a in r:
                if a.element == 'H': continue
                prot_coords.append(a.get_vector().get_array())
                prot_info.append((ch.id, r.get_id()[1]))

    na_coords, na_info = [], []
    na_rna, na_dna = False, False
    for ch in m:
        if ch.id not in nc: continue
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in NA_RES: continue
            if rn in RNA_RES: na_rna = True
            if rn in DNA_RES: na_dna = True
            for a in r:
                if a.element == 'H': continue
                cls = classify(a.get_name(), rn)
                if cls is None: continue
                na_coords.append(a.get_vector().get_array())
                na_info.append((ch.id, r.get_id()[1], cls))

    del s, m  # free structure memory
    gc.collect()

    if not prot_coords or not na_coords:
        raise ValueError("Empty atom lists")

    prot_arr = np.array(prot_coords, dtype=np.float32)
    na_arr = np.array(na_coords, dtype=np.float32)
    del prot_coords, na_coords
    gc.collect()

    # KDTree query — memory efficient
    tree = cKDTree(na_arr)
    counts = {"bb": 0, "sr": 0, "2oh": 0, "base": 0}
    pres, nres = set(), set()

    for i, coord in enumerate(prot_arr):
        neighbors = tree.query_ball_point(coord, cutoff)
        for j in neighbors:
            cls = na_info[j][2]
            counts[cls] += 1
            pres.add(prot_info[i])
            nres.add((na_info[j][0], na_info[j][1]))

    del prot_arr, na_arr, tree
    gc.collect()

    na_type = "mixed" if (na_rna and na_dna) else ("RNA" if na_rna else ("DNA" if na_dna else "?"))
    n_tot = sum(counts.values())
    spec = counts["base"] + counts["2oh"]
    gen = counts["bb"]
    si = spec / (spec + gen) if (spec + gen) > 0 else float('nan')
    bbr = counts["base"] / counts["bb"] if counts["bb"] > 0 else float('inf')

    return {
        "pdb_id": pdb_id, "na_type": na_type, "n_total": n_tot,
        "n_bb": counts["bb"], "n_sr": counts["sr"], "n_2oh": counts["2oh"],
        "n_base": counts["base"], "n_pres": len(pres), "n_nres": len(nres),
        "SI": si, "bb_frac": counts["bb"]/n_tot if n_tot else 0,
        "base_frac": counts["base"]/n_tot if n_tot else 0,
        "oh_frac": counts["2oh"]/n_tot if n_tot else 0,
        "bbr": bbr,
    }


def compute_si_split(pdb_path, pdb_id="X", cutoff=4.0):
    """compute SI separately for RNA and DNA contacts within the same structure.

    for mixed RNA+DNA structures (e.g., RNAP elongation complexes, ribosomes
    with DNA spacers), the original compute_si combines all NA contacts together.
    this function separates contacts to RNA-type residues from contacts to
    DNA-type residues, computing independent SI for each.

    returns dict with 'rna' and 'dna' keys, each containing the SI result dict
    (or None if that NA type is absent). also includes 'combined' for reference.
    """
    from Bio.PDB import PDBParser
    from scipy.spatial import cKDTree

    s = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    m = s[0]

    # detect chain types
    pc = set()  # protein chains
    rna_chains = set()  # chains with RNA residues
    dna_chains = set()  # chains with DNA residues
    for ch in m:
        chain_rna, chain_dna, chain_prot = 0, 0, 0
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn in AA_RES: chain_prot += 1
            elif rn in RNA_RES: chain_rna += 1
            elif rn in DNA_RES: chain_dna += 1
        if chain_prot > 0:
            pc.add(ch.id)
        # classify NA chain by majority type
        if chain_rna > 0 and chain_rna >= chain_dna:
            rna_chains.add(ch.id)
        elif chain_dna > 0:
            dna_chains.add(ch.id)

    if not pc:
        raise ValueError(f"no protein chains found")
    if not rna_chains and not dna_chains:
        raise ValueError(f"no NA chains found")

    # collect protein atoms
    prot_coords, prot_info = [], []
    for ch in m:
        if ch.id not in pc: continue
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in AA_RES: continue
            if r.get_id()[0] not in (' ', 'A'): continue
            for a in r:
                if a.element == 'H': continue
                prot_coords.append(a.get_vector().get_array())
                prot_info.append((ch.id, r.get_id()[1]))

    # collect NA atoms, tagged by whether they're RNA or DNA
    na_coords, na_info = [], []
    for ch in m:
        if ch.id not in (rna_chains | dna_chains): continue
        na_type_tag = "RNA" if ch.id in rna_chains else "DNA"
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in NA_RES: continue
            for a in r:
                if a.element == 'H': continue
                cls = classify(a.get_name(), rn)
                if cls is None: continue
                na_coords.append(a.get_vector().get_array())
                # add na_type_tag to distinguish RNA vs DNA contacts
                na_info.append((ch.id, r.get_id()[1], cls, na_type_tag))

    del s, m
    gc.collect()

    if not prot_coords or not na_coords:
        raise ValueError("empty atom lists")

    prot_arr = np.array(prot_coords, dtype=np.float32)
    na_arr = np.array(na_coords, dtype=np.float32)
    del prot_coords, na_coords
    gc.collect()

    # KDTree query — count contacts separately for RNA and DNA
    tree = cKDTree(na_arr)
    rna_counts = {"bb": 0, "sr": 0, "2oh": 0, "base": 0}
    dna_counts = {"bb": 0, "sr": 0, "2oh": 0, "base": 0}
    rna_pres, dna_pres = set(), set()
    rna_nres, dna_nres = set(), set()

    for i, coord in enumerate(prot_arr):
        neighbors = tree.query_ball_point(coord, cutoff)
        for j in neighbors:
            cls = na_info[j][2]
            na_type_tag = na_info[j][3]
            if na_type_tag == "RNA":
                rna_counts[cls] += 1
                rna_pres.add(prot_info[i])
                rna_nres.add((na_info[j][0], na_info[j][1]))
            else:
                dna_counts[cls] += 1
                dna_pres.add(prot_info[i])
                dna_nres.add((na_info[j][0], na_info[j][1]))

    del prot_arr, na_arr, tree
    gc.collect()

    def _build_result(counts, pres, nres, pdb_id, na_type):
        n_tot = sum(counts.values())
        if n_tot == 0:
            return None
        spec = counts["base"] + counts["2oh"]
        gen = counts["bb"]
        si = spec / (spec + gen) if (spec + gen) > 0 else float('nan')
        bbr = counts["base"] / counts["bb"] if counts["bb"] > 0 else float('inf')
        return {
            "pdb_id": pdb_id, "na_type": na_type, "n_total": n_tot,
            "n_bb": counts["bb"], "n_sr": counts["sr"], "n_2oh": counts["2oh"],
            "n_base": counts["base"], "n_pres": len(pres), "n_nres": len(nres),
            "SI": si, "bb_frac": counts["bb"]/n_tot if n_tot else 0,
            "base_frac": counts["base"]/n_tot if n_tot else 0,
            "oh_frac": counts["2oh"]/n_tot if n_tot else 0,
            "bbr": bbr,
        }

    result = {
        "rna": _build_result(rna_counts, rna_pres, rna_nres, pdb_id, "RNA"),
        "dna": _build_result(dna_counts, dna_pres, dna_nres, pdb_id, "DNA"),
        "rna_chains": sorted(rna_chains),
        "dna_chains": sorted(dna_chains),
    }
    return result


def compute_residue_contacts(pdb_path, pdb_id="X", cutoff=4.0):
    """compute per-residue NA contact profiles.

    for each protein residue within cutoff distance of any NA atom, records:
    - residue identity (chain, resnum, resname, one-letter code)
    - contact counts per category (backbone, base, 2oh, sugar)
    - dominant contact category
    - detailed contact list with distances
    - NA type contacted (RNA, DNA, or mixed)

    returns:
        residue_data: list of dicts, one per contacting protein residue
        summary: dict with aggregate SI info (same as compute_si output)
    """
    from Bio.PDB import PDBParser
    from scipy.spatial import cKDTree

    s = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    m = s[0]

    # detect chain types
    pc, nc = set(), set()
    for ch in m:
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn in AA_RES: pc.add(ch.id); break
            elif rn in NA_RES: nc.add(ch.id); break

    if not pc or not nc:
        raise ValueError(f"no prot({pc}) or NA({nc}) chains")

    # collect protein atoms — store resname alongside chain/resnum
    prot_coords, prot_info = [], []
    for ch in m:
        if ch.id not in pc: continue
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in AA_RES: continue
            if r.get_id()[0] not in (' ', 'A'): continue
            for a in r:
                if a.element == 'H': continue
                prot_coords.append(a.get_vector().get_array())
                # store (chain, resnum, resname, atom_name) for protein
                prot_info.append((ch.id, r.get_id()[1], rn, a.get_name()))

    # collect NA atoms — store full detail for contact reporting
    na_coords, na_info = [], []
    na_rna, na_dna = False, False
    for ch in m:
        if ch.id not in nc: continue
        for r in ch:
            rn = r.get_resname().strip().upper()
            if rn not in NA_RES: continue
            if rn in RNA_RES: na_rna = True
            if rn in DNA_RES: na_dna = True
            na_type_tag = "RNA" if rn in RNA_RES else "DNA"
            for a in r:
                if a.element == 'H': continue
                cls = classify(a.get_name(), rn)
                if cls is None: continue
                na_coords.append(a.get_vector().get_array())
                # store (chain, resnum, resname, atom_name, atom_class, na_type_tag)
                na_info.append((ch.id, r.get_id()[1], rn, a.get_name(), cls, na_type_tag))

    if not prot_coords or not na_coords:
        raise ValueError("empty atom lists")

    prot_arr = np.array(prot_coords, dtype=np.float32)
    na_arr = np.array(na_coords, dtype=np.float32)

    # KDTree query — collect per-residue contacts with distances
    tree = cKDTree(na_arr)

    # per-residue accumulators keyed by (chain, resnum)
    res_counts = defaultdict(lambda: {"bb": 0, "base": 0, "2oh": 0, "sr": 0})
    res_contacts = defaultdict(list)  # detailed contact list
    res_na_types = defaultdict(set)   # which NA types this residue contacts
    res_identity = {}                 # (chain, resnum) -> resname

    # aggregate counts for summary
    agg_counts = {"bb": 0, "sr": 0, "2oh": 0, "base": 0}
    all_pres, all_nres = set(), set()

    for i, coord in enumerate(prot_arr):
        neighbors = tree.query_ball_point(coord, cutoff)
        if not neighbors:
            continue

        p_chain, p_resnum, p_resname, p_aname = prot_info[i]
        res_key = (p_chain, p_resnum)
        res_identity[res_key] = p_resname

        for j in neighbors:
            na_chain, na_resnum, na_resname, na_aname, na_cls, na_tag = na_info[j]

            # compute actual distance
            dist = float(np.linalg.norm(prot_arr[i] - na_arr[j]))

            # per-residue tallies
            res_counts[res_key][na_cls] += 1
            res_contacts[res_key].append(
                (na_chain, na_resnum, na_resname, na_aname, na_cls, round(dist, 2))
            )
            res_na_types[res_key].add(na_tag)

            # aggregate tallies
            agg_counts[na_cls] += 1
            all_pres.add(res_key)
            all_nres.add((na_chain, na_resnum))

    del prot_arr, na_arr, tree
    gc.collect()

    # build per-residue result list, sorted by chain then resnum
    residue_data = []
    for res_key in sorted(res_counts.keys()):
        p_chain, p_resnum = res_key
        p_resname = res_identity[res_key]
        counts = res_counts[res_key]
        n_total = sum(counts.values())

        # determine dominant contact category
        dominant = max(counts, key=counts.get)

        # determine NA type for this residue
        na_types_seen = res_na_types[res_key]
        if len(na_types_seen) > 1:
            na_type = "mixed"
        else:
            na_type = next(iter(na_types_seen))

        residue_data.append({
            'chain': p_chain,
            'resnum': p_resnum,
            'resname': p_resname,
            'aa1': AA_3TO1.get(p_resname, '?'),
            'n_bb': counts['bb'],
            'n_base': counts['base'],
            'n_2oh': counts['2oh'],
            'n_sr': counts['sr'],
            'n_total': n_total,
            'dominant': dominant,
            'contacts': res_contacts[res_key],
            'na_type': na_type,
        })

    # build aggregate summary (same format as compute_si output)
    na_type_overall = "mixed" if (na_rna and na_dna) else ("RNA" if na_rna else ("DNA" if na_dna else "?"))
    n_tot = sum(agg_counts.values())
    spec = agg_counts["base"] + agg_counts["2oh"]
    gen = agg_counts["bb"]
    si = spec / (spec + gen) if (spec + gen) > 0 else float('nan')
    bbr = agg_counts["base"] / agg_counts["bb"] if agg_counts["bb"] > 0 else float('inf')

    summary = {
        "pdb_id": pdb_id, "na_type": na_type_overall, "n_total": n_tot,
        "n_bb": agg_counts["bb"], "n_sr": agg_counts["sr"],
        "n_2oh": agg_counts["2oh"], "n_base": agg_counts["base"],
        "n_pres": len(all_pres), "n_nres": len(all_nres),
        "SI": si, "bb_frac": agg_counts["bb"]/n_tot if n_tot else 0,
        "base_frac": agg_counts["base"]/n_tot if n_tot else 0,
        "oh_frac": agg_counts["2oh"]/n_tot if n_tot else 0,
        "bbr": bbr,
    }

    return residue_data, summary


if __name__ == "__main__":
    structs = [
        # RNA specialists
        ("1URN", "U1A RRM + RNA", "RNA-specialist"),
        ("1AQ3", "MS2 coat + RNA", "RNA-specialist"),
        ("1ASY", "AspRS + tRNA", "RNA-specialist"),
        ("1EFW", "RNaseIII + dsRNA", "RNA-specialist"),
        # RNA ancient/ribosomal
        ("1DK1", "L11 + rRNA", "RNA-ribosomal"),
        ("1FJG", "S15 + 16S rRNA", "RNA-ribosomal"),
        # Dual - RNA mode (Alba domain = PF01918)
        # note: 1J1U is TyrRS (PF00579), NOT Alba. use 3IAB (Pop6=Alba+RNA) instead.
        ("3IAB", "Pop6/Alba + RNA", "DUAL-RNA"),
        # DNA sequence-specific
        ("1LMB", "Lambda rep + DNA", "DNA-seq-spec"),
        ("3CRO", "Cro rep + DNA", "DNA-seq-spec"),
        ("1HDD", "Homeodomain + DNA", "DNA-seq-spec"),
        # DNA structural/repair
        ("1TSR", "p53 + DNA", "DNA-structural"),
        ("1BPX", "UDG + DNA", "DNA-repair"),
        # Dual - DNA mode (Alba domain = PF01918)
        # note: 2BQ2 has no protein chains. 3U6Y requires symmetry expansion (see task 3).
        # 3U6Y asymmetric unit has DNA 4.35A from protein; symm mates bring it to 2.47A.
        # result file: /tmp/nasbench/3U6Y_complex.pdb (symmetry-expanded)
        ("3U6Y", "Alba2 + dsDNA (needs symm)", "DUAL-DNA"),
        # CSD family
        ("3PF4", "CspB + ssDNA", "CSD-DNA"),
        # additional dual binders
        ("2HAX", "CSD + ssDNA", "CSD-DNA"),
    ]

    print(f"\n{'PDB':<6} {'NA':>4} {'Category':<18} {'#C':>5} {'BB%':>6} {'Base%':>6} {'2OH%':>5} {'B:BB':>5} {'SI':>6}")
    print("-" * 82)

    results = {}
    for pdb_id, desc, cat in structs:
        try:
            path = fetch(pdb_id)
            r = compute_si(path, pdb_id)
            bbr_s = f"{r['bbr']:.1f}" if r['bbr'] != float('inf') else "inf"
            print(f"{pdb_id:<6} {r['na_type']:>4} {cat:<18} {r['n_total']:>5} {r['bb_frac']:>6.1%} {r['base_frac']:>6.1%} {r['oh_frac']:>5.1%} {bbr_s:>5} {r['SI']:>6.3f}")
            results[pdb_id] = (cat, r)
            gc.collect()
        except Exception as e:
            print(f"{pdb_id:<6} {'?':>4} {cat:<18} ERR: {str(e)[:45]}")

    # Summaries
    rna = [(p,c,r) for p,(c,r) in results.items() if r['na_type']=='RNA']
    dna = [(p,c,r) for p,(c,r) in results.items() if r['na_type']=='DNA']

    print(f"\n--- RNA complexes ranked by SI ---")
    for p,c,r in sorted(rna, key=lambda x: x[2]['SI'], reverse=True):
        print(f"  {p} {c:<18} SI={r['SI']:.3f}")

    print(f"\n--- DNA complexes ranked by SI ---")
    for p,c,r in sorted(dna, key=lambda x: x[2]['SI'], reverse=True):
        print(f"  {p} {c:<18} SI={r['SI']:.3f}")

    # Paired Alba (corrected: 3IAB for RNA, 3U6Y needs symmetry expansion for DNA)
    if "3IAB" in results:
        a_rna = results["3IAB"][1]
        print(f"\n--- PAIRED: Alba (PF01918) ---")
        print(f"  + RNA (3IAB): SI={a_rna['SI']:.3f} BB={a_rna['bb_frac']:.1%} Base={a_rna['base_frac']:.1%}")
        print(f"  note: 3U6Y DNA requires symmetry expansion → SI=0.370 (see results/nasbench_paired_si_v2.tsv)")

    # Paired CSD if available
    if "2HAX" in results:
        print(f"\n--- CSD + ssDNA ---")
        r = results["2HAX"][1]
        print(f"  SI={r['SI']:.3f} BB={r['bb_frac']:.1%} Base={r['base_frac']:.1%} 2OH={r['oh_frac']:.1%}")
