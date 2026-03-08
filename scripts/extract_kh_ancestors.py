#!/usr/bin/env python3
"""
extract ancestral sequences from IQ-TREE 3 .state file.
parse per-site marginal posterior probabilities, build ML ancestor
sequences, and identify low-confidence sites.

also runs ProNA2020 on all tip + ancestor sequences if available.
"""

import os
import sys
from collections import defaultdict

STATE_FILE = "results/asr/kh/kh_asr.state"
TREE_FILE = "results/asr/kh/kh_asr.treefile"
ALN_FILE = "results/asr/kh/kh_aligned.fasta"
OUTDIR = "results/asr/kh"

AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")


def parse_state_file(path):
    """parse kh_asr.state → {node_id: [(site, ml_state, pp_dict), ...]}"""
    nodes = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Node\t"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 23:
                continue
            node = parts[0]
            site = int(parts[1])
            state = parts[2]
            probs = {AA_ORDER[i]: float(parts[3 + i]) for i in range(20)}
            # handle gap states — IQ-TREE reconstructs gaps at alignment positions
            if state == "-":
                pp = 1.0  # gap is certain at this position
            else:
                pp = probs.get(state, 0.0)
            nodes[node].append((site, state, pp, probs))
    return dict(nodes)


def build_ancestor_sequences(nodes):
    """build ML sequences and PP vectors for each internal node."""
    ancestors = {}
    for node, sites in nodes.items():
        seq = ""
        pps = []
        for site, state, pp, probs in sorted(sites, key=lambda x: x[0]):
            seq += state
            pps.append(pp)
        ancestors[node] = {"seq": seq, "pps": pps, "n_sites": len(seq)}
    return ancestors


def parse_alignment(path):
    """parse FASTA alignment → {header: seq}."""
    seqs = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:]
                seqs[current] = ""
            elif current:
                seqs[current] += line
    return seqs


def identify_root_and_deep_nodes(tree_str):
    """identify the root node and find which internal nodes are deepest.

    the tree is rooted at Node1 (outermost parentheses).
    from the newick, we can identify:
    - Node1 = root
    - Node2 = bacteria+archaea split (bootstrap 100)
    - nodes with the deepest branching points
    """
    # extract all node names from the tree
    import re
    node_names = sorted(set(re.findall(r'(Node\d+)', tree_str)))
    return node_names


def main():
    os.chdir("/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination")

    # --- parse ancestral states ---
    print("parsing ancestral state file...")
    nodes = parse_state_file(STATE_FILE)
    print(f"  found {len(nodes)} internal nodes")

    # --- build ML ancestor sequences ---
    print("building ML ancestor sequences...")
    ancestors = build_ancestor_sequences(nodes)

    # --- identify nodes ---
    with open(TREE_FILE) as f:
        tree_str = f.read().strip()
    node_names = identify_root_and_deep_nodes(tree_str)
    print(f"  internal nodes: {', '.join(node_names)}")

    # --- write ancestor FASTA ---
    anc_fasta = os.path.join(OUTDIR, "kh_ancestors.fasta")
    with open(anc_fasta, "w") as f:
        for node in sorted(ancestors.keys(), key=lambda x: int(x.replace("Node", ""))):
            a = ancestors[node]
            # remove gap columns (positions where ML state is '-' shouldn't happen
            # but IQ-TREE sometimes outputs them)
            seq = a["seq"].replace("-", "")
            f.write(f">{node}\n{seq}\n")
    print(f"  wrote {len(ancestors)} ancestors to {anc_fasta}")

    # --- write PP table ---
    pp_file = os.path.join(OUTDIR, "kh_ancestors_pp.tsv")
    with open(pp_file, "w") as f:
        f.write("node\tsite\taa\tpp\tconfident\n")
        for node in sorted(ancestors.keys(), key=lambda x: int(x.replace("Node", ""))):
            node_data = nodes[node]
            for site, state, pp, probs in sorted(node_data, key=lambda x: x[0]):
                confident = "yes" if pp >= 0.7 else "NO"
                f.write(f"{node}\t{site}\t{state}\t{pp:.5f}\t{confident}\n")
    print(f"  wrote PP table to {pp_file}")

    # --- summary of reconstruction quality per node ---
    print("\n=== RECONSTRUCTION QUALITY ===")
    print(f"{'Node':<10} {'Sites':>5} {'PP≥0.7':>6} {'PP≥0.9':>6} {'PP<0.5':>6} {'Mean PP':>8}")
    print("-" * 48)

    # focus on key nodes
    for node in sorted(ancestors.keys(), key=lambda x: int(x.replace("Node", ""))):
        a = ancestors[node]
        pps = a["pps"]
        n = len(pps)
        high = sum(1 for p in pps if p >= 0.7)
        vhigh = sum(1 for p in pps if p >= 0.9)
        low = sum(1 for p in pps if p < 0.5)
        mean_pp = sum(pps) / n if n > 0 else 0
        print(f"{node:<10} {n:>5} {high:>6} {vhigh:>6} {low:>6} {mean_pp:>8.4f}")

    # --- identify the root node (Node1) and key deep ancestors ---
    print("\n=== ROOT NODE SEQUENCE (Node1) ===")
    root = ancestors.get("Node1")
    if root:
        print(f"  sequence ({root['n_sites']} aa): {root['seq']}")
        low_pp_sites = [(i+1, root["pps"][i]) for i in range(len(root["pps"])) if root["pps"][i] < 0.7]
        if low_pp_sites:
            print(f"  ambiguous sites (PP<0.7): {len(low_pp_sites)}")
            for site, pp in low_pp_sites:
                node_data = nodes["Node1"]
                site_data = [d for d in node_data if d[0] == site][0]
                top3 = sorted(site_data[3].items(), key=lambda x: -x[1])[:3]
                top3_str = ", ".join(f"{aa}={p:.3f}" for aa, p in top3)
                print(f"    site {site}: ML={site_data[1]} PP={pp:.3f} | top3: {top3_str}")

    # --- now load tip sequences for combined FASTA ---
    print("\n=== BUILDING COMBINED FASTA FOR ProNA2020 ===")
    tip_seqs = parse_alignment(ALN_FILE)
    combined_fasta = os.path.join(OUTDIR, "kh_all_for_prona.fasta")
    with open(combined_fasta, "w") as f:
        # tips first (ungapped)
        for header, seq in tip_seqs.items():
            ungapped = seq.replace("-", "")
            f.write(f">tip|{header}\n{ungapped}\n")
        # then ancestors
        for node in sorted(ancestors.keys(), key=lambda x: int(x.replace("Node", ""))):
            seq = ancestors[node]["seq"].replace("-", "")
            f.write(f">ancestor|{node}\n{seq}\n")
    n_total = len(tip_seqs) + len(ancestors)
    print(f"  wrote {n_total} sequences ({len(tip_seqs)} tips + {len(ancestors)} ancestors)")
    print(f"  file: {combined_fasta}")

    # --- check the critical alignment columns from Track 2 ---
    crit_file = os.path.join(OUTDIR, "kh_critical_columns.tsv")
    if os.path.exists(crit_file):
        print("\n=== ANCESTRAL STATES AT CRITICAL SPECIFICITY COLUMNS ===")
        print("(from Track 2 NAS-Bench per-residue analysis)")

        # read critical columns
        import csv
        critical = []
        with open(crit_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                critical.append(row)

        # the specificity map gives alignment columns between the two PDB sequences
        # we need to map those to the phylogenetic alignment columns
        # for now, report the ancestral state at the root for each critical position
        spec_file = os.path.join(OUTDIR, "kh_specificity_map.tsv")
        if os.path.exists(spec_file):
            print("\nspecificity residue analysis from kh_specificity_map.tsv:")
            with open(spec_file) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    cls = row.get("classification", "")
                    if cls in ("rna_specific", "dna_specific", "generalist"):
                        rna_rn = row.get("rna_resname", "-")
                        dna_rn = row.get("dna_resname", "-")
                        rna_prof = row.get("rna_contact_profile", "")
                        dna_prof = row.get("dna_contact_profile", "")
                        col = row.get("aln_col", "?")
                        print(f"  col {col:>3}: {cls:<15} RNA={rna_rn}({rna_prof}) DNA={dna_rn}({dna_prof})")

    print("\n=== DONE ===")
    print(f"all output files in {OUTDIR}/")
    print("next: run ProNA2020 on kh_all_for_prona.fasta")


if __name__ == "__main__":
    main()
