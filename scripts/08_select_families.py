#!/usr/bin/env python3
"""select domain families for AF3/Boltz-2 structural prediction.

selection criteria:
  1. LUCA or pre-LUCA age (Wehbi classification)
  2. annotated as NA-binding (RNA, DNA, or both in PDB)
  3. has structural representatives in PDB (protein-NA complexes)
  4. has representatives across domains of life (Bac + Arc minimum)

priority tiers:
  tier 1: experimentally dual (RNA+DNA PDB structures) — direct thesis test
  tier 2: RNA-only or DNA-only with high structure coverage — controls
  tier 3: thesis-critical domains (HhH, ancient vocabulary) even if sparse

the key experimental design:
  for each selected domain, AF3/Boltz-2 will predict:
    - protein + poly-U RNA (10-mer)
    - protein + poly-dT DNA (10-mer)
  comparing ipTM/interface scores between RNA and DNA predictions
  tests whether "RNA-annotated" domains also bind DNA (and vice versa)

outputs:
  results/phase2_selected_families.tsv
"""

import os
import csv
import json
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')

# ============================================================
# family selection with rationale
# based on PDB survey, functional classification, and thesis relevance
# ============================================================

# format: (pfam_id, tier, selection_rationale)
SELECTED_FAMILIES = [
    # === TIER 1: EXPERIMENTALLY DUAL (RNA + DNA in PDB) ===
    # these are the core thesis test — domains that already show dual binding
    # experimentally but are annotated as single-substrate specialists

    ('PF00575', 1, 'S1/OB-fold: 308 RNA + 257 DNA structures; '
     'universal ancient fold; strongest dual-binding evidence'),
    ('PF00270', 1, 'DEAD-box helicase: 284 RNA + 56 DNA; '
     'universal RNA processing; Alva fragment 3 (P-loop)'),
    ('PF02171', 1, 'Piwi/Argonaute: 75 RNA + 11 DNA; '
     'RNA-guided DNA cleavage; direct test of dual mechanism'),
    ('PF00013', 1, 'KH_1: 23 RNA + 10 DNA; '
     'classic RNA-binding domain with unexpected DNA structures'),
    ('PF01336', 1, 'tRNA_bind/OB-fold: 14 RNA + 11 DNA; '
     'nearly equal binding; pre-LUCA age; in aaRS'),
    ('PF00398', 1, 'RrnaAD/KsgA: 27 RNA + 15 DNA; '
     'rRNA methylase with DNA structures; Alva fragment 2 (HhH)'),
    ('PF01479', 1, 'S4: 1689 RNA + 68 DNA; '
     'ribosomal protein S4; Alva fragment 10; ancient vocabulary'),
    ('PF01131', 1, 'Topoisom_bac: 8 RNA + 27 DNA; '
     'topology enzyme; primarily DNA but RNA structures exist'),

    # === TIER 2: STRONG RNA OR DNA SPECIALISTS (controls) ===
    # these serve as positive controls — if AF3 also predicts dual binding
    # for known specialists, that weakens the signal

    # RNA specialists (expect RNA >> DNA in AF3)
    ('PF01472', 2, 'PUA: 29 RNA + 0 DNA; '
     'RNA modification scaffold; clean RNA-only control'),
    ('PF01509', 2, 'TruB_N: 30 RNA + 0 DNA; '
     'pseudouridine synthase; clean RNA-only control'),
    ('PF00448', 2, 'SRP54: 17 RNA + 0 DNA; '
     'signal recognition particle; translation targeting'),
    ('PF00588', 2, 'SpoU_methylase: 7 RNA + 0 DNA; '
     'rRNA 2-O-ribose methyltransferase; RNA-only'),

    # DNA specialists (expect DNA >> RNA in AF3)
    ('PF00580', 2, 'UvrD-helicase: 0 RNA + 32 DNA; '
     'DNA repair helicase; clean DNA-only control'),
    ('PF00589', 2, 'Phage_integrase: 0 RNA + 33 DNA; '
     'site-specific recombination; DNA-only control'),
    ('PF00136', 2, 'DNA_pol_B: 9 RNA + 196 DNA; '
     'DNA polymerase B; primarily DNA; RNA from primers'),

    # === TIER 3: THESIS-CRITICAL DESPITE SPARSE PDB ===
    # these are the most important for the ancient generalism hypothesis

    ('PF00633', 3, 'HHH/helix-hairpin-helix: 0 RNA + 8 DNA; '
     'Weil-Ktorza 2025 "ambidextrous" motif; Alva fragment 2; '
     'DNA-only in PDB but predicted generalist — key thesis test'),
    ('PF00009', 3, 'GTP_EFTU: 390 RNA + 14 DNA; '
     'pre-LUCA EF-Tu; core translation; Alva fragment 3 (P-loop)'),
    ('PF00488', 3, 'MutS_V: 1 RNA + 43 DNA; '
     'mismatch repair with 1 RNA structure; interesting anomaly'),
    ('PF00133', 3, 'tRNA-synt_1: 16 RNA + 0 DNA; '
     'pre-LUCA class I aaRS; Rossmann fold (Alva fragment 7)'),

    # === TIER 3b: ADDITIONAL TRANSLATION CONTROLS ===
    ('PF00152', 3, 'tRNA-synt_2: 12 RNA + 0 DNA; '
     'class II aaRS core; translation control'),
    ('PF00347', 3, 'Ribosomal_L6: 1663 RNA + 53 DNA; '
     'ribosomal protein; mostly ribosome co-crystals'),
]


def main():
    # load functional classification
    func_path = os.path.join(PROJECT_DIR, 'results',
                             'phase1_functional_composition.tsv')
    func = {}
    with open(func_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            func[row['pfam_id']] = row

    # load PDB survey
    pdb_path = os.path.join(PROJECT_DIR, 'results', 'pdb_na_complex_survey.tsv')
    pdb_data = {}
    if os.path.exists(pdb_path):
        with open(pdb_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pdb_data[row['pfam_id']] = row

    # load Alva mapping
    alva_path = os.path.join(PROJECT_DIR, 'results',
                             'phase1_alva_luca_mapping.tsv')
    alva_pfams = set()
    if os.path.exists(alva_path):
        with open(alva_path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                # extract Pfam IDs from mapped_domains field
                for part in row.get('mapped_domains', '').split('; '):
                    if part.startswith('PF'):
                        alva_pfams.add(part.split('(')[0])

    # write selection table
    out_path = os.path.join(PROJECT_DIR, 'results',
                            'phase2_selected_families.tsv')
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow([
            'pfam_id', 'pfam_name', 'age_class', 'binding_type',
            'functional_category', 'tier', 'pdb_rna_count', 'pdb_dna_count',
            'has_alva_fragment', 'selection_rationale',
        ])

        for pfam_id, tier, rationale in SELECTED_FAMILIES:
            finfo = func.get(pfam_id, {})
            pinfo = pdb_data.get(pfam_id, {})

            writer.writerow([
                pfam_id,
                finfo.get('pfam_name', '?'),
                finfo.get('age_class', '?'),
                finfo.get('binding_type', '?'),
                finfo.get('functional_category', '?'),
                tier,
                pinfo.get('rna_count', '?'),
                pinfo.get('dna_count', '?'),
                pfam_id in alva_pfams,
                rationale,
            ])

    print(f"wrote {len(SELECTED_FAMILIES)} selected families to {out_path}")

    # print summary
    print(f"\n{'=' * 70}")
    print(f"PHASE 2: SELECTED FAMILIES FOR STRUCTURAL PREDICTION")
    print(f"{'=' * 70}")

    from collections import Counter
    tier_counts = Counter(t for _, t, _ in SELECTED_FAMILIES)
    for tier in sorted(tier_counts):
        print(f"\n  tier {tier}: {tier_counts[tier]} families")
        for pfam_id, t, rationale in SELECTED_FAMILIES:
            if t == tier:
                finfo = func.get(pfam_id, {})
                name = finfo.get('pfam_name', '?')
                binding = finfo.get('binding_type', '?')
                age = finfo.get('age_class', '?')
                has_alva = '*' if pfam_id in alva_pfams else ' '
                print(f"    {has_alva} {pfam_id} {name:20s} {age:8s} "
                      f"{binding:4s} {rationale[:60]}")

    # experimental design summary
    n = len(SELECTED_FAMILIES)
    print(f"\n{'=' * 70}")
    print(f"EXPERIMENTAL DESIGN")
    print(f"{'=' * 70}")
    print(f"  families selected: {n}")
    print(f"  predictions per family: ~9 species × 2 substrates (RNA, DNA) = ~18")
    print(f"  total predictions: ~{n * 18}")
    print(f"  plus ~{n * 3} known-substrate controls = ~{n * 18 + n * 3} total")

    # binding type composition
    rna_sel = sum(1 for pid, _, _ in SELECTED_FAMILIES
                  if func.get(pid, {}).get('binding_type') == 'RNA')
    dna_sel = sum(1 for pid, _, _ in SELECTED_FAMILIES
                  if func.get(pid, {}).get('binding_type') == 'DNA')
    print(f"\n  annotated RNA-binding: {rna_sel}")
    print(f"  annotated DNA-binding: {dna_sel}")
    print(f"  ratio: {rna_sel}:{dna_sel} (annotated RNA:DNA)")

    dual_in_pdb = sum(1 for pid, _, _ in SELECTED_FAMILIES
                      if pdb_data.get(pid, {}).get('is_dual') == 'DUAL')
    print(f"  experimentally dual (RNA+DNA in PDB): {dual_in_pdb}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
