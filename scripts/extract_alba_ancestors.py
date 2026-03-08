#!/usr/bin/env python3
"""
extract ancestral sequences from IQ-TREE 3 .state file for Alba domain.
parse per-site marginal posterior probabilities, build ML ancestor
sequences, and identify low-confidence sites.

also builds combined FASTA for ProNA2020.
"""

import os
import sys
import re
from collections import defaultdict

STATE_FILE = "results/asr/alba/alba_asr.state"
TREE_FILE = "results/asr/alba/alba_asr.treefile"
ALN_FILE = "results/asr/alba/alba_aligned.fasta"
OUTDIR = "results/asr/alba"

AA_ORDER = list("ARNDCQEGHILKMFPSTWYV")


def parseStateFile(path):
    """parse alba_asr.state -> {node_id: [(site, ml_state, pp, probs), ...]}"""
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


def buildAncestorSequences(nodes):
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


def parseAlignment(path):
    """parse FASTA alignment -> {header: seq}."""
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


def classifyTreeNodes(tree_str):
    """extract information about which tips descend from which internal nodes.
    returns:
      node_names: list of all Node names
      node_tips: dict {node: [list of tip names under this node]}
    """
    node_names = sorted(set(re.findall(r'(Node\d+)', tree_str)))
    return node_names


def findNodeTaxonomy(tree_str, tip_seqs):
    """try to identify which nodes correspond to key clades.
    parse the newick to find which tips are under each internal node.

    for Alba, the key clades are:
    - root (Node1)
    - base of Archaea
    - base of Crenarchaeota (DNA-binding lineage)
    - base of Eukaryota
    """
    # identify which tips belong to which domain of life
    archaea_tips = []
    eukaryota_tips = []
    crenarchaeota_tips = []

    for tip in tip_seqs:
        parts = tip.split('|')
        if len(parts) >= 3:
            dol = parts[2]
            if dol == 'Archaea':
                archaea_tips.append(tip)
                # identify crenarchaeota (Thermoproteota in NCBI taxonomy)
                # these include Sulfolobus, Aeropyrum, Pyrobaculum, Thermofilum, etc.
                org = parts[0].lower()
                if any(x in org for x in ['aeropyrum', 'sulfolobus', 'saccharolobus',
                                           'sulfo', 'metallosphaera', 'ignicoccus',
                                           'hyperthermus', 'staphylothermus',
                                           'thermofilum', 'caldivirga', 'pyrobaculum',
                                           'nanoarchaeum', 'cenarchaeum', 'nitrosopumilus']):
                    # note: some of these are Thaumarchaeota/DPANN, not Crenarchaeota sensu stricto
                    pass
            elif dol == 'Eukaryota':
                eukaryota_tips.append(tip)

    print(f"  archaea tips: {len(archaea_tips)}")
    print(f"  eukaryota tips: {len(eukaryota_tips)}")

    return archaea_tips, eukaryota_tips


def main():
    os.chdir("/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination")

    # --- parse ancestral states ---
    print("parsing ancestral state file...")
    nodes = parseStateFile(STATE_FILE)
    print(f"  found {len(nodes)} internal nodes")

    # --- build ML ancestor sequences ---
    print("building ML ancestor sequences...")
    ancestors = buildAncestorSequences(nodes)

    # --- identify nodes ---
    with open(TREE_FILE) as f:
        tree_str = f.read().strip()
    node_names = classifyTreeNodes(tree_str)
    print(f"  internal nodes: {', '.join(node_names)}")

    # --- parse tip alignment for taxonomy info ---
    tip_seqs = parseAlignment(ALN_FILE)
    archaea_tips, eukaryota_tips = findNodeTaxonomy(tree_str, tip_seqs)

    # --- write ancestor FASTA (ungapped) ---
    anc_fasta = os.path.join(OUTDIR, "alba_ancestors.fasta")
    with open(anc_fasta, "w") as f:
        for node in sorted(ancestors.keys(), key=lambda x: int(x.replace("Node", ""))):
            a = ancestors[node]
            seq = a["seq"].replace("-", "")
            f.write(f">{node}\n{seq}\n")
    print(f"  wrote {len(ancestors)} ancestors to {anc_fasta}")

    # --- write PP table ---
    pp_file = os.path.join(OUTDIR, "alba_ancestors_pp.tsv")
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
    print(f"{'Node':<10} {'Sites':>5} {'PP>=0.7':>7} {'PP>=0.9':>7} {'PP<0.5':>7} {'Mean PP':>8}")
    print("-" * 50)

    for node in sorted(ancestors.keys(), key=lambda x: int(x.replace("Node", ""))):
        a = ancestors[node]
        pps = a["pps"]
        n = len(pps)
        high = sum(1 for p in pps if p >= 0.7)
        vhigh = sum(1 for p in pps if p >= 0.9)
        low = sum(1 for p in pps if p < 0.5)
        mean_pp = sum(pps) / n if n > 0 else 0
        print(f"{node:<10} {n:>5} {high:>7} {vhigh:>7} {low:>7} {mean_pp:>8.4f}")

    # --- root node analysis ---
    print("\n=== ROOT NODE SEQUENCE (Node1) ===")
    root = ancestors.get("Node1")
    if root:
        print(f"  full sequence ({root['n_sites']} columns): {root['seq']}")
        ungapped = root['seq'].replace('-', '')
        print(f"  ungapped ({len(ungapped)} aa): {ungapped}")

        low_pp_sites = [(i+1, root["pps"][i]) for i in range(len(root["pps"])) if root["pps"][i] < 0.7]
        high_pp = sum(1 for p in root["pps"] if p >= 0.7)
        print(f"  confident sites (PP>=0.7): {high_pp}/{root['n_sites']} ({100*high_pp/root['n_sites']:.1f}%)")

        if low_pp_sites:
            print(f"  ambiguous sites (PP<0.7): {len(low_pp_sites)}")
            for site, pp in low_pp_sites:
                node_data = nodes["Node1"]
                site_data = [d for d in node_data if d[0] == site][0]
                top3 = sorted(site_data[3].items(), key=lambda x: -x[1])[:3]
                top3_str = ", ".join(f"{aa}={p:.3f}" for aa, p in top3)
                print(f"    site {site}: ML={site_data[1]} PP={pp:.3f} | top3: {top3_str}")

    # --- identify key nodes from tree topology ---
    # Node1 is the root (outermost split)
    # We need to figure out which node is the MRCA of archaea
    # from the tree: the archaea clade is under Node12 (91% bootstrap)
    print("\n=== KEY INTERNAL NODES (from tree topology) ===")
    print("  reading tree structure to identify key clades...")

    # the tree structure from iqtree shows:
    # Node1 = root (split: yeast/fungi vs everything else)
    # Node12 = archaeal clade (91% bootstrap) that also includes some eukaryotic sequences
    # this is common — the unrooted tree can place the root anywhere
    # let's look at node-by-node to understand the topology

    # print key nodes and their reconstructed sequences
    key_nodes = ["Node1", "Node2", "Node3", "Node12", "Node13"]
    print(f"\n  key node reconstructions:")
    for nname in key_nodes:
        if nname in ancestors:
            a = ancestors[nname]
            ungapped = a["seq"].replace("-", "")
            mean_pp = sum(a["pps"]) / len(a["pps"]) if a["pps"] else 0
            high_pp = sum(1 for p in a["pps"] if p >= 0.7)
            pct = 100 * high_pp / len(a["pps"]) if a["pps"] else 0
            print(f"    {nname}: {len(ungapped)} aa, mean PP={mean_pp:.3f}, "
                  f"sites PP>=0.7: {high_pp}/{len(a['pps'])} ({pct:.0f}%)")
            print(f"      seq: {ungapped[:60]}{'...' if len(ungapped) > 60 else ''}")

    # --- build combined FASTA for ProNA2020 ---
    print("\n=== BUILDING COMBINED FASTA FOR ProNA2020 ===")
    combined_fasta = os.path.join(OUTDIR, "alba_all_for_prona.fasta")
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

    print("\n=== DONE ===")
    print(f"all output files in {OUTDIR}/")
    print("next: run ProNA2020 on alba_all_for_prona.fasta")


if __name__ == "__main__":
    main()
