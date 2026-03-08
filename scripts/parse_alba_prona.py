#!/usr/bin/env python3
"""
parse ProNA2020 raw output for Alba domain and produce:
  1. alba_prona_all_nodes.tsv — clean table for all tips + ancestors
  2. summary report identifying key evolutionary transitions
"""

import os
import re

PRONA_RAW = "results/asr/alba/alba_prona_raw_output.txt"
TREE_FILE = "results/asr/alba/alba_asr.treefile"
OUTDIR = "results/asr/alba"


def parsePronaOutput(path):
    """parse ProNA2020 raw output -> list of dicts."""
    results = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            # skip header, separator, warnings
            if not line or line.startswith('-') or line.startswith('Protein') or \
               line.startswith('/') or line.startswith('loading') or 'Warning' in line or \
               'Deprecated' in line or 'DeprecationWarning' in line or 'UserWarning' in line:
                continue

            # expected format:
            # name   len  dna_prob dna_bind  rna_prob rna_bind  pro_prob pro_bind  nuc_prob
            parts = line.split()
            if len(parts) < 9:
                continue

            # the protein name may contain spaces — but our format doesn't
            # try to parse from the end
            try:
                nuc_prob = float(parts[-1])
                pro_bind = parts[-2]
                pro_prob = float(parts[-3])
                rna_bind = parts[-4]
                rna_prob = float(parts[-5])
                dna_bind = parts[-6]
                dna_prob = float(parts[-7])
                seq_len = int(parts[-8])
                name = ' '.join(parts[:-8])
            except (ValueError, IndexError):
                continue

            results.append({
                'name': name,
                'seq_len': seq_len,
                'dna_prob': dna_prob,
                'dna_bind': dna_bind,
                'rna_prob': rna_prob,
                'rna_bind': rna_bind,
                'pro_prob': pro_prob,
                'pro_bind': pro_bind,
                'nuc_prob': nuc_prob,
            })

    return results


def classifyNode(name):
    """classify a result as tip or ancestor, and extract metadata."""
    if name.startswith('tip|'):
        rest = name[4:]
        parts = rest.split('|')
        organism = parts[0].replace('_', ' ') if len(parts) >= 1 else 'Unknown'
        accession = parts[1] if len(parts) >= 2 else ''
        domain_of_life = parts[2] if len(parts) >= 3 else 'Unknown'
        has_pdb = parts[3] if len(parts) >= 4 else 'no'
        return 'tip', organism, accession, domain_of_life, has_pdb
    elif name.startswith('ancestor|'):
        node = name.replace('ancestor|', '')
        return 'ancestor', node, '', '', ''
    else:
        return 'unknown', name, '', '', ''


def main():
    os.chdir("/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination")

    # --- parse ProNA2020 output ---
    print("parsing ProNA2020 raw output...")
    results = parsePronaOutput(PRONA_RAW)
    print(f"  parsed {len(results)} results")

    tips = [r for r in results if r['name'].startswith('tip|')]
    ancestors = [r for r in results if r['name'].startswith('ancestor|')]
    print(f"  tips: {len(tips)}, ancestors: {len(ancestors)}")

    # --- write clean TSV ---
    tsv_path = os.path.join(OUTDIR, "alba_prona_all_nodes.tsv")
    with open(tsv_path, "w") as f:
        f.write("type\tname\torganism\taccession\tdomain_of_life\thas_pdb\t"
                "seq_len\tdna_prob\tdna_bind\trna_prob\trna_bind\t"
                "pro_prob\tpro_bind\tnuc_prob\tDI\n")

        for r in results:
            node_type, label, acc, dol, pdb = classifyNode(r['name'])

            # compute discrimination index (DI) = RNA_prob - DNA_prob
            di = r['rna_prob'] - r['dna_prob']

            if node_type == 'tip':
                f.write(f"tip\t{label}\t{label}\t{acc}\t{dol}\t{pdb}\t"
                        f"{r['seq_len']}\t{r['dna_prob']:.4f}\t{r['dna_bind']}\t"
                        f"{r['rna_prob']:.4f}\t{r['rna_bind']}\t"
                        f"{r['pro_prob']:.4f}\t{r['pro_bind']}\t"
                        f"{r['nuc_prob']:.4f}\t{di:.4f}\n")
            elif node_type == 'ancestor':
                f.write(f"ancestor\t{label}\t{label}\t\t\t\t"
                        f"{r['seq_len']}\t{r['dna_prob']:.4f}\t{r['dna_bind']}\t"
                        f"{r['rna_prob']:.4f}\t{r['rna_bind']}\t"
                        f"{r['pro_prob']:.4f}\t{r['pro_bind']}\t"
                        f"{r['nuc_prob']:.4f}\t{di:.4f}\n")

    print(f"  wrote {tsv_path}")

    # ================================================================
    # SUMMARY REPORT
    # ================================================================
    print()
    print("=" * 80)
    print("ALBA DOMAIN — ProNA2020 SUMMARY")
    print("=" * 80)

    # --- tips by domain of life ---
    print("\n--- TIP SEQUENCES ---")
    print(f"{'Organism':<50s} {'DoL':<10s} {'PDB':<6s} "
          f"{'P_DNA':>7s} {'P_RNA':>7s} {'DI':>7s} {'P_nuc':>7s}")
    print("-" * 100)

    archaea_tips = []
    eukaryota_tips = []

    for r in tips:
        _, label, acc, dol, pdb = classifyNode(r['name'])
        di = r['rna_prob'] - r['dna_prob']
        short_org = label[:50]
        print(f"{short_org:<50s} {dol:<10s} {pdb:<6s} "
              f"{r['dna_prob']:>7.4f} {r['rna_prob']:>7.4f} {di:>+7.4f} {r['nuc_prob']:>7.4f}")

        if dol == 'Archaea':
            archaea_tips.append(r)
        elif dol == 'Eukaryota':
            eukaryota_tips.append(r)

    # --- means by domain ---
    print(f"\n  mean DNA_prob / RNA_prob / DI:")
    for label, group in [("Archaea", archaea_tips), ("Eukaryota", eukaryota_tips)]:
        if group:
            mean_dna = sum(r['dna_prob'] for r in group) / len(group)
            mean_rna = sum(r['rna_prob'] for r in group) / len(group)
            mean_di = mean_rna - mean_dna
            print(f"    {label} (n={len(group)}): DNA={mean_dna:.4f}, RNA={mean_rna:.4f}, DI={mean_di:+.4f}")

    # --- key ancestral nodes ---
    print("\n--- KEY ANCESTRAL NODES ---")
    print(f"{'Node':<10s} {'P_DNA':>7s} {'P_RNA':>7s} {'DI':>7s} {'P_nuc':>7s} {'Interpretation'}")
    print("-" * 80)

    # read tree to identify what each node represents
    with open(TREE_FILE) as f:
        tree_str = f.read().strip()

    # node annotations based on tree topology analysis
    # from the tree:
    # Node1 = root (fungi vs everything)
    # Node2 = everything except fungi
    # Node3 = metazoa/nematostella + big clade
    # Node12 = archaea clade (91% bootstrap) — THIS IS THE KEY ARCHAEAL ANCESTOR
    # Node13 = deep archaea (Methanopyrus + rest)
    # Node19 = Aeropyrum (the DNA-binder!)
    # Node40 = trypanosomatids + plants
    # Node43 = fungi (Sclerotinia + Chaetomium)

    key_interpretations = {
        'Node1': 'ROOT — LUCA Alba ancestor',
        'Node2': 'non-fungal clade',
        'Node3': 'metazoa + archaea + protists',
        'Node12': 'ARCHAEAL ANCESTOR (91% BS)',
        'Node13': 'deep archaea (without M. kandleri Q8TXF9)',
        'Node14': 'core archaea',
        'Node17': 'Thermoproteota + euryarchaeota',
        'Node18': 'Thermoproteota ancestor',
        'Node19': 'Aeropyrum clade (DNA-binding!)',
        'Node40': 'trypanosomatids + plants',
        'Node41': 'trypanosomatid ancestor',
        'Node43': 'fungi (Sclerotinia + Chaetomium)',
        'Node4': 'metazoa ancestor',
        'Node5': 'mouse + C. elegans',
    }

    for r in ancestors:
        node = r['name'].replace('ancestor|', '')
        if node in key_interpretations:
            di = r['rna_prob'] - r['dna_prob']
            interp = key_interpretations[node]
            print(f"{node:<10s} {r['dna_prob']:>7.4f} {r['rna_prob']:>7.4f} "
                  f"{di:>+7.4f} {r['nuc_prob']:>7.4f} {interp}")

    # --- all ancestors sorted by node number ---
    print("\n--- ALL ANCESTRAL NODES ---")
    print(f"{'Node':<10s} {'P_DNA':>7s} {'P_RNA':>7s} {'DI':>7s} {'P_nuc':>7s}")
    print("-" * 50)

    for r in sorted(ancestors, key=lambda x: int(x['name'].replace('ancestor|Node', ''))):
        node = r['name'].replace('ancestor|', '')
        di = r['rna_prob'] - r['dna_prob']
        print(f"{node:<10s} {r['dna_prob']:>7.4f} {r['rna_prob']:>7.4f} "
              f"{di:>+7.4f} {r['nuc_prob']:>7.4f}")

    # --- the money question: does DI change from root to Aeropyrum? ---
    print("\n" + "=" * 80)
    print("KEY EVOLUTIONARY TRANSITIONS")
    print("=" * 80)

    # find the path from root (Node1) to Aeropyrum DNA-binders (Node19)
    print("\npath from root to Aeropyrum DNA-binding clade:")
    path_nodes = ['Node1', 'Node3', 'Node12', 'Node13', 'Node14',
                  'Node17', 'Node18', 'Node19']

    node_data = {}
    for r in ancestors:
        node = r['name'].replace('ancestor|', '')
        node_data[node] = r

    print(f"  {'Node':<10s} {'P_DNA':>7s} {'P_RNA':>7s} {'DI':>7s} {'Meaning'}")
    print(f"  {'-'*60}")
    for node in path_nodes:
        if node in node_data:
            r = node_data[node]
            di = r['rna_prob'] - r['dna_prob']
            meaning = key_interpretations.get(node, '')
            print(f"  {node:<10s} {r['dna_prob']:>7.4f} {r['rna_prob']:>7.4f} "
                  f"{di:>+7.4f} {meaning}")

    # compare actual Aeropyrum sequences
    print("\nactual Aeropyrum tip sequences:")
    for r in tips:
        if 'Aeropyrum' in r['name']:
            _, label, acc, dol, pdb = classifyNode(r['name'])
            di = r['rna_prob'] - r['dna_prob']
            print(f"  {acc} (PDB:{pdb}) P_DNA={r['dna_prob']:.4f} P_RNA={r['rna_prob']:.4f} DI={di:+.4f}")

    # compare with known RNA-binding tips
    print("\nknown RNA-binding tip sequences (Pop6, Rpp25):")
    for r in tips:
        name = r['name']
        if 'P53218' in name or 'Q91WE3' in name:
            _, label, acc, dol, pdb = classifyNode(name)
            di = r['rna_prob'] - r['dna_prob']
            print(f"  {acc} (PDB:{pdb}) P_DNA={r['dna_prob']:.4f} P_RNA={r['rna_prob']:.4f} DI={di:+.4f}")

    print("\n" + "=" * 80)
    print("INTERPRETATION")
    print("=" * 80)

    # compute root DI
    root = node_data.get('Node1')
    aeropyrum_node = node_data.get('Node19')
    arch_ancestor = node_data.get('Node12')

    if root:
        root_di = root['rna_prob'] - root['dna_prob']
        print(f"\n  root DI:              {root_di:+.4f} (P_RNA={root['rna_prob']:.4f}, P_DNA={root['dna_prob']:.4f})")
    if arch_ancestor:
        arch_di = arch_ancestor['rna_prob'] - arch_ancestor['dna_prob']
        print(f"  archaeal ancestor DI: {arch_di:+.4f} (P_RNA={arch_ancestor['rna_prob']:.4f}, P_DNA={arch_ancestor['dna_prob']:.4f})")
    if aeropyrum_node:
        aero_di = aeropyrum_node['rna_prob'] - aeropyrum_node['dna_prob']
        print(f"  Aeropyrum clade DI:   {aero_di:+.4f} (P_RNA={aeropyrum_node['rna_prob']:.4f}, P_DNA={aeropyrum_node['dna_prob']:.4f})")

    print()
    if root and aeropyrum_node:
        root_di = root['rna_prob'] - root['dna_prob']
        aero_di = aeropyrum_node['rna_prob'] - aeropyrum_node['dna_prob']
        if root_di > 0 and aero_di < root_di:
            print("  --> the DI decreases from root toward crenarchaeota,")
            print("      consistent with ancestral RNA-preference shifting toward DNA-binding")
        elif root_di > 0 and aero_di >= root_di:
            print("  --> the DI does NOT decrease toward crenarchaeota")
            print("      the RNA preference persists even in the DNA-binding lineage")
        else:
            print("  --> root DI is negative or zero; need careful interpretation")

    print("\n=== DONE ===")


if __name__ == "__main__":
    main()
