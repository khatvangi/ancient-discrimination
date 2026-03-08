#!/usr/bin/env python3
"""compute intersections between LUCA domains, RNA-binding domains,
DNA-binding domains, and our 51 GREEN families.

uses three RNA-binding tiers:
  - go_only: 487 pfam domains (GO:0003723 descendants, conservative lower bound)
  - rbpworld: 916 pfam domains (RBPWorld/EuRBPDB curated, primary estimate)
  - union: 968 pfam domains (both combined, upper bound)

key intersections:
  A = LUCA × RNA-binding         → headline number
  B = LUCA × DNA-binding
  C = LUCA × dual-binding (RNA + DNA)
  D = pre-LUCA × RNA-binding     → strongest generalism candidates
  E = GREEN-51 × LUCA           → project viability gate
  F = GREEN-51 × LUCA × RNA     → our families that are LUCA + RNA-binding

note on Wehbi LUCA count: the paper reports 969 LUCA pfams, but the robust
classification (ClassifiedPFAMs.csv from their GitHub) gives 871, with ~98
reclassified as "unclassifiable." we use the 871 robust classification.

outputs:
  results/phase1_intersections.json — all counts and lists
  results/phase1_census.tsv — per-pfam table with all annotations
"""

import os
import csv
import json
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')


def loadWehbi(filepath):
    """load Wehbi age classification → dict of pfam_id → age_class."""
    pfam_age = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pfam_age[row['pfam_id']] = row['age_class']
    return pfam_age


def loadRnaBindingCombined(filepath):
    """load combined RNA-binding reference → dict of pfam_id → tier."""
    rna_tiers = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            rna_tiers[row['pfam_id']] = row['tier']
    return rna_tiers


def loadDnaBinding(filepath):
    """load DNA-binding pfam list → set of pfam_ids."""
    dna_pfams = set()
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            dna_pfams.add(row['pfam_id'])
    return dna_pfams


def loadIprToPfam(filepath, source_filter=None):
    """load IPR→Pfam mapping → set of pfam_ids.
    if source_filter given, only include rows matching that source.
    """
    pfams = set()
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if source_filter is None or source_filter in row.get('source', ''):
                pfams.add(row['pfam_id'])
    return pfams


def main():
    # load all datasets
    print("loading datasets...")

    wehbi_path = os.path.join(PROJECT_DIR, 'data', 'wehbi_luca', 'luca_pfams.tsv')
    rna_combined_path = os.path.join(PROJECT_DIR, 'data', 'rbpworld', 'rna_binding_combined.tsv')
    dna_path = os.path.join(PROJECT_DIR, 'data', 'rbpworld', 'dna_binding_pfams.tsv')
    ipr_map_path = os.path.join(PROJECT_DIR, 'data', 'ipr_to_pfam_map.tsv')

    pfam_age = loadWehbi(wehbi_path)
    rna_tiers = loadRnaBindingCombined(rna_combined_path)
    dna_pfams = loadDnaBinding(dna_path)
    green_pfams = loadIprToPfam(ipr_map_path, source_filter='green51')
    universal_pfams = loadIprToPfam(ipr_map_path)

    # classify RNA-binding into tiers
    rna_go_only = {p for p, t in rna_tiers.items() if t in ('go_only', 'both')}
    rna_rbpworld = {p for p, t in rna_tiers.items() if t in ('rbpworld_only', 'both')}
    rna_union = set(rna_tiers.keys())

    # extract LUCA and pre-LUCA sets
    luca_pfams = {p for p, a in pfam_age.items() if a == 'LUCA'}
    pre_luca_pfams = {p for p, a in pfam_age.items() if a == 'preLUCA'}
    ancient_pfams = luca_pfams | pre_luca_pfams  # LUCA-age or older

    print(f"  wehbi pfams (all ages):   {len(pfam_age)}")
    print(f"  LUCA pfams:               {len(luca_pfams)}")
    print(f"  pre-LUCA pfams:           {len(pre_luca_pfams)}")
    print(f"  RNA-binding (GO only):    {len(rna_go_only)}")
    print(f"  RNA-binding (RBPWorld):    {len(rna_rbpworld)}")
    print(f"  RNA-binding (union):       {len(rna_union)}")
    print(f"  DNA-binding:              {len(dna_pfams)}")
    print(f"  GREEN-51 pfams:           {len(green_pfams)}")
    print(f"  universal-415 pfams:      {len(universal_pfams)}")

    # ================================================================
    # compute intersections for each RNA-binding tier
    # ================================================================
    results = {}

    for tier_name, rna_set in [
        ('go_lower_bound', rna_go_only),
        ('rbpworld_primary', rna_rbpworld),
        ('union_upper_bound', rna_union),
    ]:
        tier = {}

        # A: LUCA × RNA-binding
        luca_rna = luca_pfams & rna_set
        tier['luca_rna_binding_count'] = len(luca_rna)
        tier['luca_rna_binding_fraction'] = round(len(luca_rna) / len(luca_pfams) * 100, 1)
        tier['luca_rna_binding_ids'] = sorted(luca_rna)

        # B: LUCA × DNA-binding
        luca_dna = luca_pfams & dna_pfams
        tier['luca_dna_binding_count'] = len(luca_dna)
        tier['luca_dna_binding_fraction'] = round(len(luca_dna) / len(luca_pfams) * 100, 1)

        # C: LUCA × dual-binding
        luca_dual = luca_pfams & rna_set & dna_pfams
        tier['luca_dual_binding_count'] = len(luca_dual)
        tier['luca_dual_binding_ids'] = sorted(luca_dual)

        # D: pre-LUCA × RNA-binding
        pre_luca_rna = pre_luca_pfams & rna_set
        tier['pre_luca_rna_binding_count'] = len(pre_luca_rna)
        tier['pre_luca_rna_binding_fraction'] = round(
            len(pre_luca_rna) / len(pre_luca_pfams) * 100, 1
        ) if pre_luca_pfams else 0
        tier['pre_luca_rna_binding_ids'] = sorted(pre_luca_rna)

        # E: GREEN-51 × LUCA (same regardless of RNA tier)
        green_luca = green_pfams & luca_pfams
        tier['green51_luca_overlap'] = len(green_luca)
        tier['green51_luca_ids'] = sorted(green_luca)

        # F: GREEN-51 × LUCA × RNA-binding
        green_luca_rna = green_pfams & luca_pfams & rna_set
        tier['green51_luca_rna_count'] = len(green_luca_rna)
        tier['green51_luca_rna_ids'] = sorted(green_luca_rna)

        # ancient (LUCA + pre-LUCA) × RNA-binding
        ancient_rna = ancient_pfams & rna_set
        tier['ancient_rna_binding_count'] = len(ancient_rna)
        tier['ancient_rna_binding_fraction'] = round(
            len(ancient_rna) / len(ancient_pfams) * 100, 1
        )

        results[tier_name] = tier

    # also compute tier-independent values
    luca_dna = luca_pfams & dna_pfams
    pre_luca_dna = pre_luca_pfams & dna_pfams
    green_luca = green_pfams & luca_pfams

    results['metadata'] = {
        'wehbi_total_pfams': len(pfam_age),
        'wehbi_luca_count': len(luca_pfams),
        'wehbi_pre_luca_count': len(pre_luca_pfams),
        'wehbi_luca_note': (
            'Wehbi et al. 2024 report 969 LUCA-age domains; we use their 871 '
            'robustly classified domains (excluding 98 reclassified as unclassifiable '
            'in their supplementary data on GitHub)'
        ),
        'rna_go_only_count': len(rna_go_only),
        'rna_rbpworld_count': len(rna_rbpworld),
        'rna_union_count': len(rna_union),
        'dna_binding_count': len(dna_pfams),
        'green51_pfam_count': len(green_pfams),
        'green51_total_ipr': 51,
        'green51_unmapped_note': (
            '26 of 51 GREEN InterPro entries have no Pfam member database; '
            'these cannot be intersected with Wehbi Pfam-based age classification'
        ),
        'luca_dna_binding_count': len(luca_dna),
        'luca_dna_binding_fraction': round(len(luca_dna) / len(luca_pfams) * 100, 1),
        'pre_luca_dna_binding_count': len(pre_luca_dna),
        'green51_luca_overlap': len(green_luca),
        'green51_luca_ids': sorted(green_luca),
    }

    # ================================================================
    # write phase1_intersections.json
    # ================================================================
    json_path = os.path.join(PROJECT_DIR, 'results', 'phase1_intersections.json')
    os.makedirs(os.path.dirname(json_path), exist_ok=True)
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nwrote: {json_path}")

    # ================================================================
    # write phase1_census.tsv — per-pfam table
    # ================================================================
    census_path = os.path.join(PROJECT_DIR, 'results', 'phase1_census.tsv')
    with open(census_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'pfam_id', 'age_class', 'rna_go', 'rna_rbpworld', 'rna_union',
            'dna_binding', 'in_green51', 'in_universal415'
        ])
        # iterate over all pfams that appear in any dataset
        all_pfams = set(pfam_age.keys()) | rna_union | dna_pfams | green_pfams
        for pfam_id in sorted(all_pfams):
            writer.writerow([
                pfam_id,
                pfam_age.get(pfam_id, ''),
                pfam_id in rna_go_only,
                pfam_id in rna_rbpworld,
                pfam_id in rna_union,
                pfam_id in dna_pfams,
                pfam_id in green_pfams,
                pfam_id in universal_pfams,
            ])
    print(f"wrote: {census_path} ({len(all_pfams)} rows)")

    # ================================================================
    # print headline numbers
    # ================================================================
    print("\n" + "=" * 70)
    print("HEADLINE RESULTS: LUCA × RNA-BINDING INTERSECTION")
    print("=" * 70)

    for tier_name, label in [
        ('go_lower_bound', 'GO-derived (lower bound)'),
        ('rbpworld_primary', 'RBPWorld/EuRBPDB (primary)'),
        ('union_upper_bound', 'Union (upper bound)'),
    ]:
        t = results[tier_name]
        print(f"\n  {label}:")
        print(f"    LUCA × RNA:     {t['luca_rna_binding_count']:3d} / {len(luca_pfams)} = "
              f"{t['luca_rna_binding_fraction']}%")
        print(f"    pre-LUCA × RNA: {t['pre_luca_rna_binding_count']:3d} / {len(pre_luca_pfams)} = "
              f"{t['pre_luca_rna_binding_fraction']}%")
        print(f"    ancient × RNA:  {t['ancient_rna_binding_count']:3d} / {len(ancient_pfams)} = "
              f"{t['ancient_rna_binding_fraction']}%")
        print(f"    LUCA × dual:    {t['luca_dual_binding_count']:3d}")

    print(f"\n  DNA-binding (all tiers use same set):")
    print(f"    LUCA × DNA:     {results['metadata']['luca_dna_binding_count']:3d} / {len(luca_pfams)} = "
          f"{results['metadata']['luca_dna_binding_fraction']}%")

    print(f"\n  GREEN-51 (of 25 Pfam-mappable):")
    print(f"    GREEN × LUCA:   {results['metadata']['green51_luca_overlap']}")
    print(f"    GREEN × LUCA Pfam IDs: {results['metadata']['green51_luca_ids']}")

    # ================================================================
    # checkpoint and kill switch evaluation
    # ================================================================
    print("\n" + "=" * 70)
    print("CHECKPOINT — CRITICAL GATES")
    print("=" * 70)

    primary = results['rbpworld_primary']

    if primary['luca_rna_binding_count'] > 30:
        print(f"  [PASS] LUCA × RNA count ({primary['luca_rna_binding_count']}) > 30 → proceed")
    else:
        print(f"  [FAIL] LUCA × RNA count ({primary['luca_rna_binding_count']}) <= 30 → INVESTIGATE")

    green_overlap = results['metadata']['green51_luca_overlap']
    if green_overlap >= 10:
        print(f"  [PASS] GREEN-51 × LUCA overlap ({green_overlap}) >= 10 → use our families")
    elif green_overlap >= 5:
        print(f"  [WARN] GREEN-51 × LUCA overlap ({green_overlap}) is 5-9 → marginal, consider Wehbi directly")
    else:
        print(f"  [FAIL] GREEN-51 × LUCA overlap ({green_overlap}) < 5 → PIVOT to Wehbi LUCA domains")

    # cross-validation against literature
    print(f"\n  cross-validation:")
    print(f"    Anantharaman 2002 estimates ~40-45 RNA metabolism domains in LUCA")
    print(f"    our primary estimate: {primary['luca_rna_binding_count']} LUCA RNA-binding domains")
    print(f"    Crapitto 2022 consensus LUCA: 111 GO terms, translation-dominated")

    return 0


if __name__ == '__main__':
    sys.exit(main())
