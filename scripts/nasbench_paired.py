#!/usr/bin/env python3
"""
NAS-Bench paired comparison — compute SI for the same LUCA domain family
bound to RNA vs DNA, measuring ΔSI as a specificity index.

reads: data/pdb_na_cocrystals.tsv (families with example PDB IDs)
uses: nasbench.compute_si() for each structure
outputs: results/nasbench_paired_si.tsv

logic:
  for each family with both RNA and DNA co-crystals,
  try up to 3 example PDBs per type. pick the one with
  most contacts (n_total) as the representative.
  report ΔSI = SI(RNA) - SI(DNA) for the family.

interpretation:
  ΔSI near 0 → domain uses same binding mode for both substrates (generalist)
  ΔSI > 0 → more specific contacts when binding RNA
  ΔSI < 0 → more specific contacts when binding DNA
"""

import os, sys, gc

# add scripts dir to path so we can import nasbench
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nasbench import compute_si, fetch

def load_dual_families(tsv_path):
    """load families that have both RNA and DNA co-crystal examples."""
    families = []
    with open(tsv_path) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            pfam_id = parts[0]
            pfam_name = parts[1]
            n_rna = int(parts[2])
            n_dna = int(parts[3])
            rna_pdbs = parts[5].split(',') if parts[5] != 'none' else []
            dna_pdbs = parts[6].split(',') if parts[6] != 'none' else []

            # only families with both RNA and DNA examples
            if rna_pdbs and dna_pdbs:
                families.append({
                    'pfam_id': pfam_id,
                    'pfam_name': pfam_name,
                    'n_rna': n_rna,
                    'n_dna': n_dna,
                    'rna_pdbs': rna_pdbs,
                    'dna_pdbs': dna_pdbs,
                })
    return families

def best_result(pdb_list, max_try=3):
    """try up to max_try PDBs, return the one with most contacts.
    some PDBs may fail (no protein-NA interface, download error, etc).
    """
    best = None
    tried = 0
    for pdb_id in pdb_list:
        if tried >= max_try:
            break
        tried += 1
        try:
            path = fetch(pdb_id)
            r = compute_si(path, pdb_id)
            # pick the one with most contacts for robust SI
            if best is None or r['n_total'] > best['n_total']:
                best = r
            gc.collect()
        except Exception as e:
            print(f"  warning: {pdb_id} failed — {str(e)[:60]}")
    return best

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tsv_path = os.path.join(base_dir, 'data', 'pdb_na_cocrystals.tsv')
    out_path = os.path.join(base_dir, 'results', 'nasbench_paired_si.tsv')

    families = load_dual_families(tsv_path)
    print(f"found {len(families)} families with dual RNA/DNA co-crystal evidence\n")

    # header
    header = "pfam_id\tpfam_name\trna_pdb\trna_SI\trna_bb%\trna_base%\trna_2oh%\trna_contacts"
    header += "\tdna_pdb\tdna_SI\tdna_bb%\tdna_base%\tdna_2oh%\tdna_contacts"
    header += "\tdelta_SI\tinterpretation"

    results = []

    # print table header
    print(f"{'Family':<20} {'RNA_PDB':>7} {'RNA_SI':>7} {'DNA_PDB':>7} {'DNA_SI':>7} {'ΔSI':>7} {'Interp':<15}")
    print("-" * 85)

    for fam in families:
        pfam_id = fam['pfam_id']
        pfam_name = fam['pfam_name']
        print(f"\n--- {pfam_id} ({pfam_name}) ---")

        # compute best RNA result
        print(f"  RNA: trying {fam['rna_pdbs'][:3]}...")
        rna_r = best_result(fam['rna_pdbs'], max_try=3)

        # compute best DNA result
        print(f"  DNA: trying {fam['dna_pdbs'][:3]}...")
        dna_r = best_result(fam['dna_pdbs'], max_try=3)

        if rna_r is None or dna_r is None:
            print(f"  SKIP: {'no RNA result' if rna_r is None else 'no DNA result'}")
            continue

        delta_si = rna_r['SI'] - dna_r['SI']

        # interpret ΔSI
        if abs(delta_si) < 0.05:
            interp = "GENERALIST"     # same binding mode
        elif delta_si > 0.15:
            interp = "RNA-specific"   # more specific contacts with RNA
        elif delta_si < -0.15:
            interp = "DNA-specific"   # more specific contacts with DNA
        elif delta_si > 0:
            interp = "slight-RNA"
        else:
            interp = "slight-DNA"

        print(f"  {pfam_name:<20} {rna_r['pdb_id']:>7} {rna_r['SI']:>7.3f} {dna_r['pdb_id']:>7} {dna_r['SI']:>7.3f} {delta_si:>+7.3f} {interp:<15}")

        results.append({
            'pfam_id': pfam_id,
            'pfam_name': pfam_name,
            'rna_r': rna_r,
            'dna_r': dna_r,
            'delta_si': delta_si,
            'interp': interp,
        })

    # summary table
    print(f"\n\n{'='*85}")
    print(f"SUMMARY: Paired SI comparison for {len(results)} LUCA families")
    print(f"{'='*85}")
    print(f"{'Family':<20} {'RNA_PDB':>7} {'RNA_SI':>7} {'DNA_PDB':>7} {'DNA_SI':>7} {'ΔSI':>7} {'Interp':<15}")
    print("-" * 85)
    for res in sorted(results, key=lambda x: x['delta_si'], reverse=True):
        r = res['rna_r']
        d = res['dna_r']
        print(f"{res['pfam_name']:<20} {r['pdb_id']:>7} {r['SI']:>7.3f} {d['pdb_id']:>7} {d['SI']:>7.3f} {res['delta_si']:>+7.3f} {res['interp']:<15}")

    # statistics
    generalists = [r for r in results if r['interp'] == 'GENERALIST']
    rna_spec = [r for r in results if 'RNA' in r['interp']]
    dna_spec = [r for r in results if 'DNA' in r['interp']]
    print(f"\ngeneralist (|ΔSI| < 0.05): {len(generalists)}/{len(results)}")
    print(f"rna-leaning (ΔSI > 0):     {len(rna_spec)}/{len(results)}")
    print(f"dna-leaning (ΔSI < 0):     {len(dna_spec)}/{len(results)}")

    if results:
        delta_vals = [r['delta_si'] for r in results]
        import numpy as np
        print(f"mean ΔSI: {np.mean(delta_vals):+.3f}")
        print(f"median ΔSI: {np.median(delta_vals):+.3f}")

    # write TSV output
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        f.write(header + '\n')
        for res in sorted(results, key=lambda x: x['delta_si'], reverse=True):
            r = res['rna_r']
            d = res['dna_r']
            row = f"{res['pfam_id']}\t{res['pfam_name']}"
            row += f"\t{r['pdb_id']}\t{r['SI']:.4f}\t{r['bb_frac']:.4f}\t{r['base_frac']:.4f}\t{r['oh_frac']:.4f}\t{r['n_total']}"
            row += f"\t{d['pdb_id']}\t{d['SI']:.4f}\t{d['bb_frac']:.4f}\t{d['base_frac']:.4f}\t{d['oh_frac']:.4f}\t{d['n_total']}"
            row += f"\t{res['delta_si']:.4f}\t{res['interp']}"
            f.write(row + '\n')
    print(f"\nresults written to {out_path}")


if __name__ == '__main__':
    main()
