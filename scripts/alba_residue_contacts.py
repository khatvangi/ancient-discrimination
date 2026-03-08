#!/usr/bin/env python3
"""
alba per-residue contact analysis.

compares RNA-contacting residues (3IAB: Pop6/Alba + RNA)
vs DNA-contacting residues (3U6Y: Alba2 + dsDNA, symmetry-expanded).

uses the existing compute_residue_contacts() from nasbench.py.

output files (all in results/asr/alba/):
  - alba_rna_residue_contacts.tsv
  - alba_dna_residue_contacts.tsv
  - alba_specificity_map.tsv
  - alba_critical_columns.tsv
"""

import sys, os, csv

# add scripts directory to path for nasbench import
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from nasbench import fetch, compute_residue_contacts

# project root and output dir
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
OUT_DIR = os.path.join(PROJECT_ROOT, "results", "asr", "alba")
os.makedirs(OUT_DIR, exist_ok=True)


# ── step 1: compute residue contacts ──────────────────────────────

print("=" * 60)
print("STEP 1: compute per-residue NA contacts")
print("=" * 60)

# --- RNA structure: 3IAB (Pop6/Alba chain A + RNA chain R) ---
print("\n[3IAB] fetching Pop6/Alba + RNA structure...")
path_rna = fetch("3IAB")
rna_residues, rna_summary = compute_residue_contacts(path_rna, "3IAB")

# filter to chain A only (Pop6 = the Alba domain protein)
# chain B is Pop7 (not Alba)
rna_residues = [r for r in rna_residues if r["chain"] == "A"]
print(f"  3IAB chain A (Pop6/Alba): {len(rna_residues)} contacting residues")
print(f"  aggregate SI = {rna_summary['SI']:.3f} (both chains)")
# recompute summary for chain A only
rna_chain_a_total = sum(r["n_total"] for r in rna_residues)
rna_chain_a_bb = sum(r["n_bb"] for r in rna_residues)
rna_chain_a_base = sum(r["n_base"] for r in rna_residues)
rna_chain_a_2oh = sum(r["n_2oh"] for r in rna_residues)
rna_chain_a_sr = sum(r["n_sr"] for r in rna_residues)
spec_a = rna_chain_a_base + rna_chain_a_2oh
gen_a = rna_chain_a_bb
si_chain_a = spec_a / (spec_a + gen_a) if (spec_a + gen_a) > 0 else float("nan")
print(f"  chain A only: {rna_chain_a_total} contacts, SI = {si_chain_a:.3f}")

# --- DNA structure: 3U6Y (Alba2 + dsDNA, symmetry-expanded) ---
symm_path = "/tmp/nasbench/3U6Y_complex.pdb"
if os.path.exists(symm_path):
    print(f"\n[3U6Y] using symmetry-expanded file: {symm_path}")
    path_dna = symm_path
else:
    print("\n[3U6Y] WARNING: symmetry-expanded file not found, using raw PDB")
    print("  DNA contacts may be incomplete (asymmetric unit only)")
    path_dna = fetch("3U6Y")

dna_residues, dna_summary = compute_residue_contacts(path_dna, "3U6Y")
print(f"  3U6Y: {len(dna_residues)} contacting residues across all protein chains")
print(f"  aggregate SI = {dna_summary['SI']:.3f}")

# for 3U6Y, both chains A and C are Alba2 copies in the dimer.
# keep both — they form a functional unit around the DNA.
# but note which chain each residue belongs to.
dna_chains_seen = set(r["chain"] for r in dna_residues)
for ch in sorted(dna_chains_seen):
    n = sum(1 for r in dna_residues if r["chain"] == ch)
    print(f"  chain {ch}: {n} contacting residues")


# ── save residue contact tables ───────────────────────────────────

def save_residue_table(residues, filepath, pdb_id):
    """write per-residue contact table as TSV."""
    fields = [
        "pdb_id", "chain", "resnum", "resname", "aa1",
        "n_bb", "n_base", "n_2oh", "n_sr", "n_total",
        "dominant", "na_type",
    ]
    with open(filepath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in residues:
            row = {k: r.get(k, "") for k in fields}
            row["pdb_id"] = pdb_id
            w.writerow(row)
    print(f"  saved: {filepath} ({len(residues)} rows)")

rna_out = os.path.join(OUT_DIR, "alba_rna_residue_contacts.tsv")
dna_out = os.path.join(OUT_DIR, "alba_dna_residue_contacts.tsv")

save_residue_table(rna_residues, rna_out, "3IAB")
save_residue_table(dna_residues, dna_out, "3U6Y")


# ── step 2: identify specificity residues ─────────────────────────

print("\n" + "=" * 60)
print("STEP 2: identify specificity residues")
print("=" * 60)

# build lookup: (aa1, resnum) -> contact profile for each structure
# NOTE: these are DIFFERENT proteins (Pop6 vs Ape10b2), so direct residue
# number comparison is meaningless. we need sequence alignment (step 3).
# for now, classify residues within each structure independently.

def classify_residue(r):
    """classify a residue's contact type."""
    if r["n_2oh"] > 0:
        return "2OH-contacting"  # RNA-specific (contacts ribose 2'OH)
    elif r["n_base"] > 0 and r["n_bb"] == 0:
        return "base-only"  # sequence-reading
    elif r["n_bb"] > 0 and r["n_base"] == 0:
        return "backbone-only"  # non-specific
    elif r["n_base"] > 0 and r["n_bb"] > 0:
        return "base+backbone"  # mixed
    elif r["n_sr"] > 0:
        return "sugar-ring-only"
    else:
        return "other"

print("\n--- 3IAB (RNA) residue classifications ---")
rna_classes = {}
for r in rna_residues:
    cls = classify_residue(r)
    rna_classes[cls] = rna_classes.get(cls, 0) + 1
    r["classification"] = cls
for cls, n in sorted(rna_classes.items(), key=lambda x: -x[1]):
    print(f"  {cls}: {n}")

print("\n--- 3U6Y (DNA) residue classifications ---")
dna_classes = {}
for r in dna_residues:
    cls = classify_residue(r)
    dna_classes[cls] = dna_classes.get(cls, 0) + 1
    r["classification"] = cls
for cls, n in sorted(dna_classes.items(), key=lambda x: -x[1]):
    print(f"  {cls}: {n}")

# show 2'OH-contacting residues from 3IAB (these are RNA-specific contacts)
print("\n--- RNA-specific residues (2'OH contacts in 3IAB) ---")
oh_residues = [r for r in rna_residues if r["n_2oh"] > 0]
if oh_residues:
    for r in oh_residues:
        print(f"  {r['chain']}:{r['aa1']}{r['resnum']} — "
              f"2OH={r['n_2oh']} base={r['n_base']} bb={r['n_bb']} sr={r['n_sr']}")
else:
    print("  none (no 2'OH contacts in chain A)")


# ── step 3: sequence alignment + UniProt mapping ──────────────────

print("\n" + "=" * 60)
print("STEP 3: extract SEQRES and map to UniProt positions")
print("=" * 60)

from Bio.PDB import PDBParser
from Bio import pairwise2

parser = PDBParser(QUIET=True)

def extract_chain_sequence(pdb_path, pdb_id, chain_id):
    """extract amino acid sequence from a PDB chain (ATOM records only)."""
    from nasbench import AA_RES, AA_3TO1
    s = parser.get_structure(pdb_id, pdb_path)
    m = s[0]
    ch = m[chain_id]
    residues = []
    for r in ch:
        rn = r.get_resname().strip().upper()
        if rn in AA_RES and r.get_id()[0] in (" ", "A"):
            residues.append((r.get_id()[1], rn, AA_3TO1.get(rn, "?")))
    return residues


def fetch_uniprot_seq(acc):
    """fetch sequence from UniProt."""
    import urllib.request
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            lines = resp.read().decode().strip().split("\n")
            header = lines[0]
            seq = "".join(lines[1:])
            return header, seq
    except Exception as e:
        print(f"  WARNING: could not fetch {acc}: {e}")
        return None, None


# -- 3IAB chain A (Pop6, P53218) --
print("\n[3IAB chain A] Pop6 / P53218 (S. cerevisiae)")
rna_chain_res = extract_chain_sequence(path_rna, "3IAB", "A")
rna_pdb_seq = "".join(r[2] for r in rna_chain_res)
rna_pdb_resnums = [r[0] for r in rna_chain_res]
print(f"  PDB chain A: {len(rna_chain_res)} residues, range {rna_pdb_resnums[0]}-{rna_pdb_resnums[-1]}")

_, pop6_seq = fetch_uniprot_seq("P53218")
if pop6_seq:
    print(f"  UniProt P53218: {len(pop6_seq)} residues")

    # align PDB sequence to UniProt sequence
    alns = pairwise2.align.globalms(pop6_seq, rna_pdb_seq, 2, -1, -5, -0.5,
                                     one_alignment_only=True)
    if alns:
        aln = alns[0]
        # build mapping: PDB index -> UniProt position
        rna_pdb2uniprot = {}
        up_pos = 0  # counter for UniProt positions (1-based)
        pdb_idx = 0  # counter for PDB residues
        for a, b in zip(aln.seqA, aln.seqB):
            if a != "-":
                up_pos += 1
            if b != "-":
                pdb_resnum = rna_pdb_resnums[pdb_idx]
                if a != "-":
                    rna_pdb2uniprot[pdb_resnum] = up_pos
                pdb_idx += 1
        print(f"  alignment mapped {len(rna_pdb2uniprot)}/{len(rna_chain_res)} residues")

# -- 3U6Y chains A,C (Alba2, Q9YAX2) --
print("\n[3U6Y] Alba2 / Q9YAX2 (A. pernix)")
# use chain A as representative (chain C is identical copy)
dna_chain_res = extract_chain_sequence(path_dna, "3U6Y", "A")
dna_pdb_seq = "".join(r[2] for r in dna_chain_res)
dna_pdb_resnums = [r[0] for r in dna_chain_res]
print(f"  PDB chain A: {len(dna_chain_res)} residues, range {dna_pdb_resnums[0]}-{dna_pdb_resnums[-1]}")

_, alba2_seq = fetch_uniprot_seq("Q9YAX2")
if alba2_seq:
    print(f"  UniProt Q9YAX2: {len(alba2_seq)} residues")

    alns = pairwise2.align.globalms(alba2_seq, dna_pdb_seq, 2, -1, -5, -0.5,
                                     one_alignment_only=True)
    if alns:
        aln = alns[0]
        dna_pdb2uniprot = {}
        up_pos = 0
        pdb_idx = 0
        for a, b in zip(aln.seqA, aln.seqB):
            if a != "-":
                up_pos += 1
            if b != "-":
                pdb_resnum = dna_pdb_resnums[pdb_idx]
                if a != "-":
                    dna_pdb2uniprot[pdb_resnum] = up_pos
                pdb_idx += 1
        print(f"  alignment mapped {len(dna_pdb2uniprot)}/{len(dna_chain_res)} residues")


# ── save critical columns table ───────────────────────────────────

print("\n--- saving critical columns ---")
crit_out = os.path.join(OUT_DIR, "alba_critical_columns.tsv")
crit_fields = [
    "pdb_source", "pdb_chain", "pdb_resnum", "pdb_resname", "aa1",
    "uniprot_id", "uniprot_position",
    "n_bb", "n_base", "n_2oh", "n_sr", "n_total",
    "dominant", "classification",
]

with open(crit_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=crit_fields, delimiter="\t")
    w.writeheader()

    # RNA contacts (3IAB chain A)
    for r in rna_residues:
        up_pos = rna_pdb2uniprot.get(r["resnum"], "")
        w.writerow({
            "pdb_source": "3IAB",
            "pdb_chain": r["chain"],
            "pdb_resnum": r["resnum"],
            "pdb_resname": r["resname"],
            "aa1": r["aa1"],
            "uniprot_id": "P53218",
            "uniprot_position": up_pos,
            "n_bb": r["n_bb"],
            "n_base": r["n_base"],
            "n_2oh": r["n_2oh"],
            "n_sr": r["n_sr"],
            "n_total": r["n_total"],
            "dominant": r["dominant"],
            "classification": r["classification"],
        })

    # DNA contacts (3U6Y — all contacting chains)
    for r in dna_residues:
        # map using chain A numbering (chains A and C are identical copies)
        up_pos = dna_pdb2uniprot.get(r["resnum"], "")
        w.writerow({
            "pdb_source": "3U6Y",
            "pdb_chain": r["chain"],
            "pdb_resnum": r["resnum"],
            "pdb_resname": r["resname"],
            "aa1": r["aa1"],
            "uniprot_id": "Q9YAX2",
            "uniprot_position": up_pos,
            "n_bb": r["n_bb"],
            "n_base": r["n_base"],
            "n_2oh": r["n_2oh"],
            "n_sr": r["n_sr"],
            "n_total": r["n_total"],
            "dominant": r["dominant"],
            "classification": r["classification"],
        })

print(f"  saved: {crit_out}")


# ── step 4: cross-structure alignment and comparison ──────────────

print("\n" + "=" * 60)
print("STEP 4: cross-structure comparison")
print("=" * 60)

# these are DIFFERENT proteins (Pop6 vs Alba2/Ape10b2) sharing the Alba fold.
# align them to each other to find which positions are equivalent.
print("\n--- aligning Pop6 (3IAB) vs Alba2 (3U6Y) ---")
cross_alns = pairwise2.align.globalms(rna_pdb_seq, dna_pdb_seq, 2, -1, -5, -0.5,
                                       one_alignment_only=True)

if cross_alns:
    cross_aln = cross_alns[0]
    print(f"  alignment score: {cross_aln.score:.1f}")
    print(f"  Pop6 length: {len(rna_pdb_seq)}, Alba2 length: {len(dna_pdb_seq)}")

    # build position mapping: Pop6 resnum -> Alba2 resnum
    rna_idx, dna_idx = 0, 0
    rna2dna_map = {}  # Pop6 resnum -> Alba2 resnum
    for a, b in zip(cross_aln.seqA, cross_aln.seqB):
        rna_rn = rna_pdb_resnums[rna_idx] if (a != "-" and rna_idx < len(rna_pdb_resnums)) else None
        dna_rn = dna_pdb_resnums[dna_idx] if (b != "-" and dna_idx < len(dna_pdb_resnums)) else None
        if a != "-" and b != "-" and rna_rn is not None and dna_rn is not None:
            rna2dna_map[rna_rn] = dna_rn
        if a != "-":
            rna_idx += 1
        if b != "-":
            dna_idx += 1

    print(f"  aligned positions: {len(rna2dna_map)}")

    # show alignment snippet (first 80 chars)
    print(f"\n  alignment (first 80 positions):")
    print(f"  Pop6 : {cross_aln.seqA[:80]}")
    print(f"  Alba2: {cross_aln.seqB[:80]}")

# build sets of contacting resnums for each structure
rna_contact_resnums = set(r["resnum"] for r in rna_residues)
# for DNA, use chain A residues only for comparison (avoid double-counting dimer)
dna_contact_resnums_A = set(r["resnum"] for r in dna_residues if r["chain"] == "A")
dna_contact_resnums_C = set(r["resnum"] for r in dna_residues if r["chain"] == "C")
dna_contact_resnums_all = dna_contact_resnums_A | dna_contact_resnums_C

print(f"\n--- contact footprint summary ---")
print(f"  3IAB (RNA): {len(rna_contact_resnums)} residues contact RNA")
print(f"  3U6Y (DNA) chain A: {len(dna_contact_resnums_A)} residues contact DNA")
print(f"  3U6Y (DNA) chain C: {len(dna_contact_resnums_C)} residues contact DNA")

# find structurally equivalent positions that contact NA in both structures
# map RNA contact positions to DNA numbering via cross-alignment
rna_mapped_to_dna = set()
for rn in rna_contact_resnums:
    if rn in rna2dna_map:
        rna_mapped_to_dna.add(rna2dna_map[rn])

overlap = rna_mapped_to_dna & dna_contact_resnums_A
rna_only = rna_mapped_to_dna - dna_contact_resnums_A
dna_only = dna_contact_resnums_A - rna_mapped_to_dna

print(f"\n--- structural equivalence (via sequence alignment) ---")
print(f"  RNA contact positions mappable to Alba2: {len(rna_mapped_to_dna)}")
print(f"  overlap (contact NA in both structures): {len(overlap)}")
print(f"  RNA-only (contact RNA but not DNA at equivalent position): {len(rna_only)}")
print(f"  DNA-only (contact DNA but not RNA at equivalent position): {len(dna_only)}")

if overlap:
    print(f"\n  overlapping positions (Alba2 numbering):")
    for dna_rn in sorted(overlap):
        # find the corresponding residues
        dna_res = [r for r in dna_residues if r["chain"] == "A" and r["resnum"] == dna_rn]
        # find the Pop6 resnum that maps to this
        rna_rn = [k for k, v in rna2dna_map.items() if v == dna_rn]
        rna_res = [r for r in rna_residues if r["resnum"] in rna_rn]
        if dna_res and rna_res:
            d = dna_res[0]
            r = rna_res[0]
            print(f"    Pop6 {r['aa1']}{r['resnum']} <-> Alba2 {d['aa1']}{d['resnum']}"
                  f"  | RNA: bb={r['n_bb']} base={r['n_base']} 2oh={r['n_2oh']}"
                  f"  | DNA: bb={d['n_bb']} base={d['n_base']}")

# show 2'OH-specific residues and whether their equivalent position contacts DNA
print(f"\n--- 2'OH-contacting residues and structural equivalents ---")
for r in rna_residues:
    if r["n_2oh"] > 0:
        dna_equiv = rna2dna_map.get(r["resnum"])
        if dna_equiv:
            dna_match = [d for d in dna_residues
                         if d["chain"] == "A" and d["resnum"] == dna_equiv]
            if dna_match:
                d = dna_match[0]
                print(f"  Pop6 {r['aa1']}{r['resnum']} (2OH={r['n_2oh']}) <-> "
                      f"Alba2 {d['aa1']}{dna_equiv} (contacts DNA: "
                      f"bb={d['n_bb']} base={d['n_base']})")
            else:
                print(f"  Pop6 {r['aa1']}{r['resnum']} (2OH={r['n_2oh']}) <-> "
                      f"Alba2 position {dna_equiv}: NO DNA contact")
        else:
            print(f"  Pop6 {r['aa1']}{r['resnum']} (2OH={r['n_2oh']}) <-> "
                  f"no alignment to Alba2")


# ── save specificity map ──────────────────────────────────────────

print("\n--- saving specificity map ---")
spec_out = os.path.join(OUT_DIR, "alba_specificity_map.tsv")
spec_fields = [
    "alba2_resnum", "alba2_aa", "pop6_resnum", "pop6_aa",
    "contacts_rna", "contacts_dna",
    "rna_2oh", "rna_base", "rna_bb",
    "dna_base", "dna_bb",
    "specificity",
]

# build reverse mapping: dna_resnum -> rna_resnum
dna2rna_map = {v: k for k, v in rna2dna_map.items()}

with open(spec_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=spec_fields, delimiter="\t")
    w.writeheader()

    # iterate over all aligned positions using Alba2 (chain A) numbering
    all_dna_positions = set(r[0] for r in dna_chain_res)
    all_mapped_positions = set(rna2dna_map.values()) | all_dna_positions

    for dna_rn in sorted(all_mapped_positions):
        # find corresponding residue info
        dna_res_info = next((r for r in dna_chain_res if r[0] == dna_rn), None)
        rna_rn = dna2rna_map.get(dna_rn)
        rna_res_info = next((r for r in rna_chain_res if r[0] == rna_rn), None) if rna_rn else None

        # find contact data
        dna_contact = next((r for r in dna_residues if r["chain"] == "A" and r["resnum"] == dna_rn), None)
        rna_contact = next((r for r in rna_residues if r["resnum"] == rna_rn), None) if rna_rn else None

        contacts_rna = rna_contact is not None
        contacts_dna = dna_contact is not None

        # classify specificity
        if contacts_rna and contacts_dna:
            specificity = "generalist"
        elif contacts_rna and not contacts_dna:
            specificity = "RNA-only"
        elif contacts_dna and not contacts_rna:
            specificity = "DNA-only"
        else:
            specificity = "non-contact"

        # refine: if 2'OH contacts exist, mark as RNA-specific even if generalist
        if contacts_rna and rna_contact["n_2oh"] > 0:
            specificity = "RNA-specific-2OH" if not contacts_dna else "generalist-2OH"

        row = {
            "alba2_resnum": dna_rn if dna_res_info else "",
            "alba2_aa": dna_res_info[2] if dna_res_info else "",
            "pop6_resnum": rna_rn if rna_res_info else "",
            "pop6_aa": rna_res_info[2] if rna_res_info else "",
            "contacts_rna": contacts_rna,
            "contacts_dna": contacts_dna,
            "rna_2oh": rna_contact["n_2oh"] if rna_contact else 0,
            "rna_base": rna_contact["n_base"] if rna_contact else 0,
            "rna_bb": rna_contact["n_bb"] if rna_contact else 0,
            "dna_base": dna_contact["n_base"] if dna_contact else 0,
            "dna_bb": dna_contact["n_bb"] if dna_contact else 0,
            "specificity": specificity,
        }
        w.writerow(row)

print(f"  saved: {spec_out}")


# ── final summary ─────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SUMMARY: Alba domain (PF01918) NA contact analysis")
print("=" * 60)

print(f"\n  RNA structure: 3IAB (Pop6/Alba, S. cerevisiae)")
print(f"    chain A contacting residues: {len(rna_residues)}")
print(f"    chain A SI: {si_chain_a:.3f}")
print(f"    2'OH-contacting residues: {len(oh_residues)}")

print(f"\n  DNA structure: 3U6Y (Alba2, A. pernix, symmetry-expanded)")
print(f"    contacting residues (all chains): {len(dna_residues)}")
print(f"    chain A contacting residues: {len(dna_contact_resnums_A)}")
print(f"    aggregate SI: {dna_summary['SI']:.3f}")

print(f"\n  cross-structure comparison (aligned positions):")
print(f"    overlapping contact positions: {len(overlap)}")
print(f"    RNA-only contact positions: {len(rna_only)}")
print(f"    DNA-only contact positions: {len(dna_only)}")

total_contact = len(overlap) + len(rna_only) + len(dna_only)
if total_contact > 0:
    pct_shared = len(overlap) / total_contact * 100
    print(f"    shared binding surface: {pct_shared:.0f}% of contact positions overlap")
    print(f"    → {'same' if pct_shared > 50 else 'different'} binding surface "
          f"for RNA vs DNA")

print(f"\n  output files:")
print(f"    {rna_out}")
print(f"    {dna_out}")
print(f"    {crit_out}")
print(f"    {spec_out}")
print()
