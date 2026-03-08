#!/usr/bin/env python3
"""
NAS-Bench paired comparison v2 — handles mixed RNA+DNA structures.

for structures with both RNA and DNA chains, uses compute_si_split
to measure contacts to each NA type separately. this fixes the confound
where mixed structures (like RNAP elongation complexes) gave identical
or misleading SI values.

reads: data/pdb_na_cocrystals.tsv
uses: nasbench.compute_si(), nasbench.compute_si_split()
outputs: results/nasbench_paired_si_v2.tsv
"""

import os, sys, gc
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nasbench import compute_si, compute_si_split, fetch


def load_dual_families(tsv_path):
    """load families that have both RNA and DNA co-crystal examples."""
    families = []
    with open(tsv_path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            pfam_id = parts[0]
            pfam_name = parts[1]
            n_rna = int(parts[2])
            n_dna = int(parts[3])
            n_both = int(parts[4])
            rna_pdbs = parts[5].split(',') if parts[5] != 'none' else []
            dna_pdbs = parts[6].split(',') if parts[6] != 'none' else []
            if rna_pdbs and dna_pdbs:
                families.append({
                    'pfam_id': pfam_id,
                    'pfam_name': pfam_name,
                    'n_rna': n_rna,
                    'n_dna': n_dna,
                    'n_both': n_both,
                    'rna_pdbs': rna_pdbs,
                    'dna_pdbs': dna_pdbs,
                })
    return families


def smart_si(pdb_list, target_na_type, max_try=3):
    """compute SI for the target NA type from a list of PDB candidates.

    strategy:
    1. try each PDB in order
    2. if the structure has ONLY the target NA type → use compute_si directly
    3. if the structure is mixed → use compute_si_split and take the target type
    4. pick the result with the most contacts for robustness
    """
    best = None
    best_method = None
    tried = 0

    for pdb_id in pdb_list:
        if tried >= max_try:
            break
        tried += 1

        try:
            path = fetch(pdb_id)

            # first check if it's mixed using compute_si_split
            split_result = compute_si_split(path, pdb_id)

            has_rna = split_result['rna'] is not None
            has_dna = split_result['dna'] is not None

            if target_na_type == "RNA":
                if has_rna and not has_dna:
                    # pure RNA structure — use original compute_si for full accuracy
                    r = compute_si(path, pdb_id)
                    method = "pure"
                elif has_rna:
                    # mixed structure — use RNA portion from split
                    r = split_result['rna']
                    method = f"split(mixed:{'+'.join(split_result['rna_chains'])}RNA,{'+'.join(split_result['dna_chains'])}DNA)"
                else:
                    print(f"    {pdb_id}: no RNA chains found, skip")
                    continue
            else:  # DNA
                if has_dna and not has_rna:
                    r = compute_si(path, pdb_id)
                    method = "pure"
                elif has_dna:
                    r = split_result['dna']
                    method = f"split(mixed:{'+'.join(split_result['rna_chains'])}RNA,{'+'.join(split_result['dna_chains'])}DNA)"
                else:
                    print(f"    {pdb_id}: no DNA chains found, skip")
                    continue

            # pick result with most contacts
            if best is None or r['n_total'] > best['n_total']:
                best = r
                best_method = method
                best['_method'] = method
                best['_pdb_id'] = pdb_id

            gc.collect()

        except Exception as e:
            print(f"    {pdb_id}: error — {str(e)[:60]}")

    return best


def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    tsv_path = os.path.join(base_dir, 'data', 'pdb_na_cocrystals.tsv')
    out_path = os.path.join(base_dir, 'results', 'nasbench_paired_si_v2.tsv')

    families = load_dual_families(tsv_path)
    print(f"found {len(families)} families with dual RNA/DNA co-crystal evidence\n")

    results = []

    for fam in families:
        pfam_id = fam['pfam_id']
        pfam_name = fam['pfam_name']
        print(f"\n{'='*60}")
        print(f"{pfam_id} ({pfam_name}) — n_rna={fam['n_rna']}, n_dna={fam['n_dna']}, n_both={fam['n_both']}")

        # compute RNA SI
        print(f"  RNA targets: {fam['rna_pdbs'][:3]}")
        rna_r = smart_si(fam['rna_pdbs'], "RNA", max_try=3)
        if rna_r:
            print(f"  → RNA: {rna_r['_pdb_id']} SI={rna_r['SI']:.3f} n={rna_r['n_total']} ({rna_r['_method']})")
        else:
            print(f"  → RNA: FAILED (no valid structure)")

        # compute DNA SI
        print(f"  DNA targets: {fam['dna_pdbs'][:3]}")
        dna_r = smart_si(fam['dna_pdbs'], "DNA", max_try=3)
        if dna_r:
            print(f"  → DNA: {dna_r['_pdb_id']} SI={dna_r['SI']:.3f} n={dna_r['n_total']} ({dna_r['_method']})")
        else:
            print(f"  → DNA: FAILED (no valid structure)")

        if rna_r is None or dna_r is None:
            print(f"  SKIP — missing {'RNA' if rna_r is None else 'DNA'} data")
            continue

        delta_si = rna_r['SI'] - dna_r['SI']

        # interpret with quality flags
        min_contacts = min(rna_r['n_total'], dna_r['n_total'])
        quality = "HIGH" if min_contacts >= 100 else ("MED" if min_contacts >= 30 else "LOW")

        if abs(delta_si) < 0.05:
            interp = "GENERALIST"
        elif delta_si > 0.15:
            interp = "RNA-specific"
        elif delta_si < -0.15:
            interp = "DNA-specific"
        elif delta_si > 0:
            interp = "slight-RNA"
        else:
            interp = "slight-DNA"

        print(f"  ΔSI = {delta_si:+.3f} → {interp} (quality: {quality})")

        results.append({
            'pfam_id': pfam_id,
            'pfam_name': pfam_name,
            'rna_pdb': rna_r['_pdb_id'],
            'rna_method': rna_r['_method'],
            'rna_si': rna_r['SI'],
            'rna_bb': rna_r['bb_frac'],
            'rna_base': rna_r['base_frac'],
            'rna_2oh': rna_r['oh_frac'],
            'rna_n': rna_r['n_total'],
            'dna_pdb': dna_r['_pdb_id'],
            'dna_method': dna_r['_method'],
            'dna_si': dna_r['SI'],
            'dna_bb': dna_r['bb_frac'],
            'dna_base': dna_r['base_frac'],
            'dna_2oh': dna_r['oh_frac'],
            'dna_n': dna_r['n_total'],
            'delta_si': delta_si,
            'quality': quality,
            'interp': interp,
        })

    # summary
    print(f"\n\n{'='*90}")
    print(f"SUMMARY: Paired NAS-Bench v2 for {len(results)} LUCA families")
    print(f"{'='*90}")
    print(f"{'Family':<18} {'RNA_PDB':>7} {'RNA_SI':>7} {'DNA_PDB':>7} {'DNA_SI':>7} {'ΔSI':>7} {'Q':>4} {'Interp':<15} {'Method'}")
    print("-" * 110)
    for r in sorted(results, key=lambda x: x['delta_si'], reverse=True):
        rna_meth = "pure" if r['rna_method'] == "pure" else "split"
        dna_meth = "pure" if r['dna_method'] == "pure" else "split"
        methods = f"RNA:{rna_meth}, DNA:{dna_meth}"
        print(f"{r['pfam_name']:<18} {r['rna_pdb']:>7} {r['rna_si']:>7.3f} {r['dna_pdb']:>7} {r['dna_si']:>7.3f} {r['delta_si']:>+7.3f} {r['quality']:>4} {r['interp']:<15} {methods}")

    # statistics
    if results:
        deltas = [r['delta_si'] for r in results]
        high_q = [r for r in results if r['quality'] == 'HIGH']
        generalists = [r for r in results if r['interp'] == 'GENERALIST']

        print(f"\nall families: mean ΔSI = {np.mean(deltas):+.3f}, median = {np.median(deltas):+.3f}")
        if high_q:
            hq_deltas = [r['delta_si'] for r in high_q]
            print(f"high-quality only ({len(high_q)}): mean ΔSI = {np.mean(hq_deltas):+.3f}, median = {np.median(hq_deltas):+.3f}")
        print(f"\ngeneralist (|ΔSI| < 0.05): {len(generalists)}/{len(results)}")
        print(f"rna-leaning (ΔSI > 0):     {len([r for r in results if r['delta_si'] > 0.05])}/{len(results)}")
        print(f"dna-leaning (ΔSI < 0):     {len([r for r in results if r['delta_si'] < -0.05])}/{len(results)}")

    # write TSV
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        cols = ['pfam_id', 'pfam_name', 'rna_pdb', 'rna_method', 'rna_si', 'rna_bb', 'rna_base', 'rna_2oh', 'rna_n',
                'dna_pdb', 'dna_method', 'dna_si', 'dna_bb', 'dna_base', 'dna_2oh', 'dna_n',
                'delta_si', 'quality', 'interp']
        f.write('\t'.join(cols) + '\n')
        for r in sorted(results, key=lambda x: x['delta_si'], reverse=True):
            vals = [str(r[c]) if isinstance(r[c], (str, int)) else f"{r[c]:.4f}" for c in cols]
            f.write('\t'.join(vals) + '\n')
    print(f"\nresults written to {out_path}")


if __name__ == '__main__':
    main()
