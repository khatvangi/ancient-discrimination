#!/usr/bin/env python3
"""
alba convergence analysis: cross-reference ancestral states (Track 1)
with specificity-determining residue positions (Track 2).

key question: does the root/deep ancestor have RNA-like or DNA-like
residues at the critical binding positions?
"""

import sys
from pathlib import Path
from collections import defaultdict

ASR_DIR = Path("results/asr/alba")

# ── step 1: parse alignment to build domain_pos → alignment_col mapping ──

def parse_alignment(fasta_path):
    """return dict: seq_id -> aligned_sequence (with gaps)"""
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

def build_dompos_to_alncol(aligned_seq):
    """map ungapped domain position (1-based) -> alignment column (1-based)"""
    mapping = {}
    dom_pos = 0
    for col_idx, char in enumerate(aligned_seq):
        if char != "-":
            dom_pos += 1
            mapping[dom_pos] = col_idx + 1  # 1-based
    return mapping

# ── step 2: define contact positions within the pfam domain ──

# pop6 (P53218): pfam domain starts at UniProt position 45
# so domain_pos = uniprot_pos - 44
POP6_OFFSET = 44

# alba2 (Q9YAX2): pfam domain starts at UniProt position 7
# so domain_pos = uniprot_pos - 6
ALBA2_OFFSET = 6

# rna contacts from critical_columns (3IAB/P53218)
# (uniprot_pos, aa, classification)
RNA_CONTACTS = [
    (51,  "K", "base-only"),
    (52,  "N", "base+backbone"),
    (54,  "N", "backbone-only"),
    (55,  "I", "backbone-only"),
    (56,  "K", "2OH-contacting"),  # KEY 2'OH RESIDUE
    (59,  "V", "base-only"),
    (60,  "N", "base-only"),
    (61,  "K", "backbone-only"),
    (64,  "K", "backbone-only"),
    (89,  "Q", "base-only"),
    (90,  "K", "base-only"),
    (93,  "S", "base+backbone"),
    (96,  "E", "base-only"),
    (97,  "I", "base-only"),
    (100, "K", "backbone-only"),
]

# contacts outside domain (uniprot 20, 129, 131, 132, 133, 134) — can't trace
RNA_OUTSIDE_DOMAIN = [20, 129, 131, 132, 133, 134]

# dna contacts from critical_columns (3U6Y/Q9YAX2)
DNA_CONTACTS = [
    (10,  "R", "backbone-only"),
    (13,  "R", "base+backbone"),
    (40,  "R", "backbone-only"),
    (41,  "G", "backbone-only"),
    (42,  "R", "sugar-ring-only"),
    (43,  "N", "backbone-only"),
    (45,  "N", "backbone-only"),
    (46,  "R", "backbone-only"),  # also contacts RNA (generalist)
]

# contacts outside domain (uniprot 86, 88)
DNA_OUTSIDE_DOMAIN = [86, 88]

# ── step 3: parse IQ-TREE .state file ──

def parse_state_file(state_path):
    """return dict: node_id -> {site_int: (ml_state, pp, full_probs_dict)}"""
    aa_order = list("ARNDCQEGHILKMFPSTWYV")
    nodes = defaultdict(dict)
    with open(state_path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Node\t"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            node_id = parts[0]
            site = int(parts[1])
            state = parts[2]
            probs = {}
            for i, aa in enumerate(aa_order):
                probs[aa] = float(parts[3 + i])
            pp = probs.get(state, 0.0)
            nodes[node_id][site] = (state, pp, probs)
    return nodes

# ── step 4: parse ancestors fasta for ML sequences ──

def parse_ancestors(fasta_path):
    """return dict: node_id -> sequence"""
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

# ── main analysis ──

def main():
    # parse alignment
    aln = parse_alignment(ASR_DIR / "alba_aligned.fasta")

    # find the two reference sequences
    pop6_key = None
    alba2_key = None
    for k in aln:
        if "P53218" in k:
            pop6_key = k
        if "Q9YAX2" in k:
            alba2_key = k

    if not pop6_key or not alba2_key:
        print("ERROR: couldn't find reference sequences in alignment")
        sys.exit(1)

    # build domain position -> alignment column mappings
    pop6_map = build_dompos_to_alncol(aln[pop6_key])
    alba2_map = build_dompos_to_alncol(aln[alba2_key])

    aln_len = len(aln[pop6_key])
    print(f"alignment length: {aln_len} columns")
    print(f"pop6 domain: {len(pop6_map)} positions (UniProt starts at {POP6_OFFSET + 1})")
    print(f"alba2 domain: {len(alba2_map)} positions (UniProt starts at {ALBA2_OFFSET + 1})")

    # map contacts to alignment columns
    print("\n" + "=" * 80)
    print("RNA CONTACTS (Pop6/3IAB) mapped to alignment columns")
    print("=" * 80)
    rna_aln_cols = {}
    for uniprot_pos, aa, classification in RNA_CONTACTS:
        dom_pos = uniprot_pos - POP6_OFFSET
        aln_col = pop6_map.get(dom_pos)
        if aln_col:
            rna_aln_cols[aln_col] = (uniprot_pos, aa, classification, dom_pos)
            tag = " *** 2'OH ***" if "2OH" in classification else ""
            print(f"  UniProt {uniprot_pos} ({aa}) -> domain pos {dom_pos} -> aln col {aln_col} [{classification}]{tag}")
        else:
            print(f"  UniProt {uniprot_pos} ({aa}) -> domain pos {dom_pos} -> NOT IN ALIGNMENT")

    print(f"\n  NOTE: {len(RNA_OUTSIDE_DOMAIN)} RNA contacts outside Pfam domain: {RNA_OUTSIDE_DOMAIN}")
    print(f"  (includes 2'OH contact at UniProt 131)")

    print("\n" + "=" * 80)
    print("DNA CONTACTS (Alba2/3U6Y) mapped to alignment columns")
    print("=" * 80)
    dna_aln_cols = {}
    for uniprot_pos, aa, classification in DNA_CONTACTS:
        dom_pos = uniprot_pos - ALBA2_OFFSET
        aln_col = alba2_map.get(dom_pos)
        if aln_col:
            dna_aln_cols[aln_col] = (uniprot_pos, aa, classification, dom_pos)
            print(f"  UniProt {uniprot_pos} ({aa}) -> domain pos {dom_pos} -> aln col {aln_col} [{classification}]")
        else:
            print(f"  UniProt {uniprot_pos} ({aa}) -> domain pos {dom_pos} -> NOT IN ALIGNMENT")

    print(f"\n  NOTE: {len(DNA_OUTSIDE_DOMAIN)} DNA contacts outside Pfam domain: {DNA_OUTSIDE_DOMAIN}")

    # identify shared vs exclusive columns
    rna_only_cols = set(rna_aln_cols.keys()) - set(dna_aln_cols.keys())
    dna_only_cols = set(dna_aln_cols.keys()) - set(rna_aln_cols.keys())
    shared_cols = set(rna_aln_cols.keys()) & set(dna_aln_cols.keys())

    print(f"\n  RNA-only alignment columns: {sorted(rna_only_cols)}")
    print(f"  DNA-only alignment columns: {sorted(dna_only_cols)}")
    print(f"  Shared alignment columns:   {sorted(shared_cols)}")

    # ── parse ancestral states ──
    state_data = parse_state_file(ASR_DIR / "alba_asr.state")
    ancestors = parse_ancestors(ASR_DIR / "alba_ancestors.fasta")

    # key nodes to analyze (from ProNA2020 analysis)
    key_nodes = [
        ("Node1",  "root (universal ancestor)"),
        ("Node12", "archaeal ancestor (91% BS, best-reconstructed)"),
        ("Node23", "crenarchaeal ancestor (Sulfolobales group)"),
        ("Node18", "desulfurococcales ancestor"),
        ("Node19", "Staphylothermus+Ignicoccus ancestor"),
    ]

    # get the all critical columns (union)
    all_critical_cols = sorted(set(rna_aln_cols.keys()) | set(dna_aln_cols.keys()))

    print("\n" + "=" * 80)
    print("ANCESTRAL STATES AT CRITICAL BINDING POSITIONS")
    print("=" * 80)

    # header for the detailed table
    print(f"\n{'col':>4} {'rna_aa':>6} {'dna_aa':>6} {'type':>10}", end="")
    for node_id, desc in key_nodes:
        print(f" | {node_id:>8}", end="")
    print()
    print("-" * (30 + 12 * len(key_nodes)))

    for col in all_critical_cols:
        rna_info = rna_aln_cols.get(col, None)
        dna_info = dna_aln_cols.get(col, None)

        rna_aa = rna_info[1] if rna_info else "-"
        dna_aa = dna_info[1] if dna_info else "-"

        if col in shared_cols:
            col_type = "shared"
        elif col in rna_only_cols:
            col_type = "RNA"
        else:
            col_type = "DNA"

        # mark 2'OH
        if rna_info and "2OH" in rna_info[2]:
            col_type += "*"

        print(f"{col:>4} {rna_aa:>6} {dna_aa:>6} {col_type:>10}", end="")

        for node_id, desc in key_nodes:
            if node_id in state_data and col in state_data[node_id]:
                state, pp, probs = state_data[node_id][col]
                # check if ancestral state matches RNA or DNA modern residue
                match = ""
                if state == rna_aa and rna_info:
                    match = "=RNA"
                elif state == dna_aa and dna_info:
                    match = "=DNA"
                print(f" | {state}({pp:.2f}){match}", end="")
            else:
                print(f" | {'?':>8}", end="")
        print()

    # ── summary per node ──
    print("\n" + "=" * 80)
    print("SUMMARY: ANCESTRAL RESIDUE IDENTITY AT SPECIFICITY POSITIONS")
    print("=" * 80)

    for node_id, desc in key_nodes:
        if node_id not in state_data:
            continue

        rna_match = 0
        dna_match = 0
        neither = 0
        total = 0
        high_pp = 0  # sites with PP >= 0.7

        # check RNA-only columns
        for col in rna_only_cols:
            rna_aa = rna_aln_cols[col][1]
            if col in state_data[node_id]:
                state, pp, probs = state_data[node_id][col]
                total += 1
                if pp >= 0.7:
                    high_pp += 1
                if state == rna_aa:
                    rna_match += 1
                else:
                    neither += 1

        # check DNA-only columns
        for col in dna_only_cols:
            dna_aa = dna_aln_cols[col][1]
            if col in state_data[node_id]:
                state, pp, probs = state_data[node_id][col]
                total += 1
                if pp >= 0.7:
                    high_pp += 1
                if state == dna_aa:
                    dna_match += 1
                else:
                    neither += 1

        # check shared columns — compare to both
        for col in shared_cols:
            rna_aa = rna_aln_cols[col][1]
            dna_aa = dna_aln_cols[col][1]
            if col in state_data[node_id]:
                state, pp, probs = state_data[node_id][col]
                total += 1
                if pp >= 0.7:
                    high_pp += 1
                if state == rna_aa:
                    rna_match += 1
                elif state == dna_aa:
                    dna_match += 1
                else:
                    neither += 1

        print(f"\n{node_id} — {desc}")
        print(f"  matches RNA specialist:  {rna_match}/{total} ({100*rna_match/total:.0f}%)")
        print(f"  matches DNA specialist:  {dna_match}/{total} ({100*dna_match/total:.0f}%)")
        print(f"  matches neither:         {neither}/{total}")
        print(f"  high-confidence (PP>=0.7): {high_pp}/{total}")

    # ── the key 2'OH residue ──
    print("\n" + "=" * 80)
    print("KEY 2'OH RESIDUE ANALYSIS")
    print("=" * 80)

    # the 2'OH contact in Pop6 is K56 (UniProt) = domain pos 12 = alignment col 13
    oh2_col = None
    for col, info in rna_aln_cols.items():
        if "2OH" in info[2]:
            oh2_col = col
            break

    if oh2_col:
        rna_aa = rna_aln_cols[oh2_col][1]
        print(f"\nPop6 2'OH residue: {rna_aa} at UniProt {rna_aln_cols[oh2_col][0]}, alignment col {oh2_col}")

        # what's at this position in Alba2 (DNA binder)?
        alba2_state = aln[alba2_key][oh2_col - 1]  # 0-based indexing
        print(f"Alba2 at same column: {alba2_state}")

        print(f"\nAncestral states at this column:")
        for node_id, desc in key_nodes:
            if node_id in state_data and oh2_col in state_data[node_id]:
                state, pp, probs = state_data[node_id][oh2_col]
                # also show probability of the RNA specialist residue (K)
                p_rna = probs.get(rna_aa, 0)
                print(f"  {node_id:>8} ({desc}): {state} (PP={pp:.3f}), P({rna_aa})={p_rna:.3f}")

    # ── also check the out-of-domain 2'OH contact ──
    print(f"\nNOTE: Pop6 also has a 2'OH contact at UniProt 131 (E)")
    print(f"  This position is OUTSIDE the Pfam domain and cannot be traced phylogenetically.")
    print(f"  It maps to the specificity_map as alba2 pos 74 (I) -> RNA-specific-2OH")

    # ── write summary table ──
    outpath = ASR_DIR / "alba_convergence_summary.tsv"
    with open(outpath, "w") as f:
        header = ["aln_col", "rna_uniprot", "rna_aa", "dna_uniprot", "dna_aa",
                   "col_type", "rna_classification", "dna_classification"]
        for node_id, _ in key_nodes:
            header.extend([f"{node_id}_state", f"{node_id}_pp"])
        f.write("\t".join(header) + "\n")

        for col in all_critical_cols:
            rna_info = rna_aln_cols.get(col, None)
            dna_info = dna_aln_cols.get(col, None)

            rna_up = str(rna_info[0]) if rna_info else ""
            rna_aa = rna_info[1] if rna_info else ""
            rna_cls = rna_info[2] if rna_info else ""
            dna_up = str(dna_info[0]) if dna_info else ""
            dna_aa = dna_info[1] if dna_info else ""
            dna_cls = dna_info[2] if dna_info else ""

            if col in shared_cols:
                col_type = "shared"
            elif col in rna_only_cols:
                col_type = "RNA-only"
            else:
                col_type = "DNA-only"

            row = [str(col), rna_up, rna_aa, dna_up, dna_aa, col_type, rna_cls, dna_cls]

            for node_id, _ in key_nodes:
                if node_id in state_data and col in state_data[node_id]:
                    state, pp, _ = state_data[node_id][col]
                    row.extend([state, f"{pp:.4f}"])
                else:
                    row.extend(["", ""])

            f.write("\t".join(row) + "\n")

    print(f"\nconvergence table saved to {outpath}")

if __name__ == "__main__":
    main()
