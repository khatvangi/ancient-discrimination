#!/usr/bin/env python3
"""map Alva et al. 2015 ancient peptide fragments onto LUCA NA-binding domains.

the 40 fragments are subfold-level motifs (30-70 aa) that appear across
multiple SCOP folds. our LUCA domains are Pfam families. the bridge is:
  fragment → SCOP fold → Pfam structural homology → LUCA age classification

three mapping strategies:
  1. ribosomal protein name → known Pfam domain
  2. SCOP fold → Pfam domains via InterPro superfamily cross-reference
  3. NA-binding annotation comparison

outputs:
  results/phase1_alva_luca_mapping.tsv — per-fragment mapping table
  results/phase1_alva_summary.json — summary statistics
"""

import os
import csv
import json
import sys
from collections import Counter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')

# ============================================================
# manual SCOP fold → LUCA Pfam mapping
# based on InterPro superfamily cross-reference (SSF → Pfam)
# each SCOP fold maps to the LUCA NA-binding Pfams that share
# that structural superfamily
# ============================================================
SCOP_TO_LUCA_PFAMS = {
    # c.37 = P-loop NTPase (SSF52540 / IPR027417)
    # the largest superfamily — GTPases, helicases, kinases
    'c.37': [
        'PF00009',  # GTP_EFTU — EF-Tu GTPase (preLUCA)
        'PF00270',  # DEAD — DEAD-box RNA helicase
        'PF00448',  # SRP54 — SRP GTPase
        'PF00488',  # MutS_V — mismatch repair ATPase
        'PF00580',  # UvrD-helicase — helicase SF-I
        'PF01926',  # MMR_HSR1 — ribosome-associated GTPase
        'PF04851',  # ResIII — restriction enzyme (preLUCA)
        'PF05191',  # ADK_lid — adenylate kinase
        'PF06733',  # DEAD_2 — DEAD-box variant
    ],

    # a.60 = SAM-dependent methyltransferase (SSF53335 / IPR029063)
    'a.60': [
        'PF00398',  # RrnaAD — rRNA adenine dimethylase
        'PF01170',  # UPF0020 — SAM-dependent MTase
        'PF01189',  # Methyltr_RsmB-F — 16S rRNA MTase
        'PF01555',  # N6_N4_Mtase — DNA methyltransferase
        'PF02475',  # Met_10 — tRNA MTase
        'PF13847',  # Methyltransf_31 — MTase (preLUCA)
    ],

    # c.2 = Rossmann fold (SSF51735 / IPR036291)
    # class I aminoacyl-tRNA synthetase catalytic domains
    'c.2': [
        'PF00133',  # tRNA-synt_1 — class I aaRS core (preLUCA)
        'PF00579',  # tRNA-synt_1b
        'PF00749',  # tRNA-synt_1c
        'PF00750',  # tRNA-synt_1d
        'PF09334',  # tRNA-synt_1g
    ],

    # b.40 = OB-fold (SSF50199 / IPR012340)
    # includes S1 domain, tRNA-binding, SNase
    'b.40': [
        'PF00565',  # SNase — staphylococcal nuclease fold
        'PF00575',  # S1 — S1 RNA-binding domain
        'PF01336',  # tRNA_anti-codon — anticodon binding (preLUCA)
        'PF01588',  # tRNA_bind — tRNA-binding domain
    ],

    # winged helix / HTH variants
    'a.4': [
        'PF01325',  # Fe_dep_repress — iron-dependent repressor HTH
        'PF09339',  # HTH_IclR — IclR transcription factor
        'PF13404',  # HTH_AsnC-type — AsnC transcription factor
    ],

    # a.24 = ferritin-like 4-helical bundle (SSF47240)
    # primarily iron-storage/redox, few NA-binders
    'a.24': [],

    # d.58 = ferredoxin-like (SSF54862 / IPR036010)
    # primarily FeS-binding; some RNA-binding via RNase H-like
    'd.58': [
        'PF02171',  # Piwi — Argonaute/PIWI (RNase H-like fold)
    ],

    # d.37 = flavodoxin-like (SSF52218 / IPR029039)
    # primarily cofactor binding; no direct NA-binding LUCA Pfams
    'd.37': [],

    # c.1 = TIM barrel (SSF51351)
    'c.1': [
        'PF01207',  # Dus — dihydrouridine synthase
    ],

    # c.66 = S-adenosylmethionine synthetase (often grouped with c.2)
    'c.66': [],

    # b.84 = meander beta (from fragment 12, L27)
    'b.84': [],

    # d.66 = S4/S5-like SH3 barrel
    'd.66': [
        'PF01479',  # S4 — S4 RNA-binding domain
    ],
}

# ============================================================
# ribosomal protein → Pfam mapping
# maps Alva's ribosomal protein annotations to specific Pfam IDs
# ============================================================
RIBOSOMAL_TO_PFAM = {
    'S19e': ['PF01090'],    # Ribosomal_S19e
    'S13':  ['PF00416'],    # Ribosomal_S13
    'L29':  ['PF00831'],    # Ribosomal_L29
    'S5':   ['PF00333'],    # Ribosomal_S5
    'S3':   ['PF00189'],    # Ribosomal_S3_C
    'S4':   ['PF01479'],    # S4 RNA-binding domain — IN LUCA LIST
    'TL5':  ['PF00673'],    # Ribosomal_L5
    'L27':  ['PF01016'],    # Ribosomal_L27
    'L7/12': ['PF00542'],   # Ribosomal_L7Ae
    # these ribosomal proteins have Pfam domains IN our LUCA NA-binding list:
    # (not in Alva fragments, but present in our LUCA set for context)
    # L4 → PF00573 (Ribosomal_L4) — LUCA
    # L6 → PF00347 (Ribosomal_L6) — LUCA
    # L10 → PF00466 (Ribosomal_L10) — LUCA
}


def loadWehbiAges(filepath):
    """load Wehbi age classification for all pfams."""
    ages = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ages[row['pfam_id']] = row['age_class']
    return ages


def loadFunctionalClassification(filepath):
    """load our functional classification of LUCA NA-binding domains."""
    classes = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            classes[row['pfam_id']] = {
                'name': row['pfam_name'],
                'age': row['age_class'],
                'binding': row['binding_type'],
                'function': row['functional_category'],
            }
    return classes


def loadAlvaFragments(filepath):
    """load Alva 2015 fragments."""
    fragments = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            fragments.append(row)
    return fragments


def main():
    wehbi_path = os.path.join(PROJECT_DIR, 'data', 'wehbi_luca', 'luca_pfams.tsv')
    func_path = os.path.join(PROJECT_DIR, 'results', 'phase1_functional_composition.tsv')
    alva_path = os.path.join(PROJECT_DIR, 'data', 'alva_fragments', 'fragments.tsv')

    pfam_ages = loadWehbiAges(wehbi_path)
    func_classes = loadFunctionalClassification(func_path)
    fragments = loadAlvaFragments(alva_path)

    # set of all LUCA+preLUCA NA-binding Pfam IDs
    luca_na_pfams = set(func_classes.keys())

    print(f"loaded {len(fragments)} Alva fragments")
    print(f"loaded {len(func_classes)} classified LUCA NA-binding domains")
    print(f"loaded {len(pfam_ages)} Wehbi age classifications")

    # ================================================================
    # map each fragment to LUCA domains
    # ================================================================
    results = []

    for frag in fragments:
        fid = frag['fragment_id']
        na_bind = frag.get('na_binding', '').strip()
        ribosomal = frag.get('ribosomal', '').strip()
        scop_folds = [s.strip() for s in frag.get('scop_folds', '').split(',')
                      if s.strip() and s.strip() != 'all_folds']

        # collect all mapped LUCA Pfam IDs for this fragment
        mapped_luca_pfams = set()
        mapping_sources = []

        # strategy 1: ribosomal protein → Pfam
        if ribosomal:
            rpfams = RIBOSOMAL_TO_PFAM.get(ribosomal, [])
            for pid in rpfams:
                if pid in luca_na_pfams:
                    mapped_luca_pfams.add(pid)
                    mapping_sources.append(f"ribosomal:{ribosomal}→{pid}")
                else:
                    # check if it's at least LUCA-age (but not in NA-binding list)
                    age = pfam_ages.get(pid, 'unknown')
                    if age in ('LUCA', 'preLUCA'):
                        mapping_sources.append(
                            f"ribosomal:{ribosomal}→{pid}(LUCA,not_NA-annotated)")
                    else:
                        mapping_sources.append(
                            f"ribosomal:{ribosomal}→{pid}(age={age})")

        # strategy 2: SCOP fold → Pfam via structural homology
        for fold in scop_folds:
            # try exact match, then parent fold
            fold_key = fold
            if fold_key not in SCOP_TO_LUCA_PFAMS:
                # try parent fold (e.g., c.37.1.20 → c.37)
                parts = fold_key.split('.')
                if len(parts) >= 2:
                    fold_key = '.'.join(parts[:2])
                if fold_key not in SCOP_TO_LUCA_PFAMS:
                    fold_key = parts[0] + '.' + parts[1] if len(parts) >= 2 else fold_key

            luca_pfams_for_fold = SCOP_TO_LUCA_PFAMS.get(fold_key, [])
            for pid in luca_pfams_for_fold:
                mapped_luca_pfams.add(pid)
                mapping_sources.append(f"SCOP:{fold}→{pid}")

        # determine binding classification of mapped domains
        rna_count = 0
        dna_count = 0
        mapped_functions = []
        for pid in sorted(mapped_luca_pfams):
            info = func_classes.get(pid, {})
            bt = info.get('binding', '?')
            fn = info.get('function', '?')
            if bt == 'RNA':
                rna_count += 1
            elif bt == 'DNA':
                dna_count += 1
            mapped_functions.append(f"{pid}({info.get('name', '?')},{bt},{fn})")

        results.append({
            'fragment_id': fid,
            'num_folds': frag['num_folds'],
            'num_superfamilies': frag['num_superfamilies'],
            'ribosomal': ribosomal,
            'alva_na_binding': na_bind,
            'scop_folds': frag.get('scop_folds', ''),
            'mapped_luca_pfams': len(mapped_luca_pfams),
            'mapped_rna': rna_count,
            'mapped_dna': dna_count,
            'mapped_domains': '; '.join(mapped_functions),
            'mapping_sources': '; '.join(mapping_sources),
        })

    # ================================================================
    # write output TSV
    # ================================================================
    out_path = os.path.join(PROJECT_DIR, 'results', 'phase1_alva_luca_mapping.tsv')
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=[
            'fragment_id', 'num_folds', 'num_superfamilies', 'ribosomal',
            'alva_na_binding', 'scop_folds', 'mapped_luca_pfams',
            'mapped_rna', 'mapped_dna', 'mapped_domains', 'mapping_sources',
        ])
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    print(f"\nwrote {len(results)} fragment mappings to {out_path}")

    # ================================================================
    # summary statistics
    # ================================================================
    mapped_any = [r for r in results if r['mapped_luca_pfams'] > 0]
    mapped_rna = [r for r in results if r['mapped_rna'] > 0]
    mapped_dna = [r for r in results if r['mapped_dna'] > 0]
    na_annotated = [r for r in results if r['alva_na_binding']]

    print(f"\n{'=' * 70}")
    print("ALVA FRAGMENT → LUCA DOMAIN MAPPING SUMMARY")
    print("=" * 70)
    print(f"  total Alva fragments:              {len(results)}")
    print(f"  fragments with LUCA domain match:  {len(mapped_any)} ({len(mapped_any)/len(results)*100:.0f}%)")
    print(f"  → mapped to LUCA RNA-binding:      {len(mapped_rna)}")
    print(f"  → mapped to LUCA DNA-binding:      {len(mapped_dna)}")
    print(f"  fragments with NA annotation:      {len(na_annotated)} ({len(na_annotated)/len(results)*100:.0f}%)")

    # cross-tabulate: Alva's NA annotation vs our LUCA mapping
    print(f"\n{'=' * 70}")
    print("CROSS-TABULATION: ALVA NA ANNOTATION vs LUCA MAPPING")
    print("=" * 70)
    print(f"  {'Alva annotation':20s} {'maps→LUCA RNA':>14s} {'maps→LUCA DNA':>14s} {'no LUCA map':>12s}")
    print(f"  {'-'*20} {'-'*14} {'-'*14} {'-'*12}")

    for annot_label, annot_values in [
        ('RNA', ['RNA']),
        ('DNA', ['DNA']),
        ('DNA+RNA', ['DNA+RNA']),
        ('none', ['']),
    ]:
        subset = [r for r in results if r['alva_na_binding'] in annot_values]
        has_rna = sum(1 for r in subset if r['mapped_rna'] > 0)
        has_dna = sum(1 for r in subset if r['mapped_dna'] > 0)
        no_map = sum(1 for r in subset if r['mapped_luca_pfams'] == 0)
        print(f"  {annot_label:20s} {has_rna:14d} {has_dna:14d} {no_map:12d}")

    # key fragments table
    print(f"\n{'=' * 70}")
    print("KEY FRAGMENTS: NA-BINDING WITH LUCA DOMAIN MATCHES")
    print("=" * 70)
    for r in results:
        if r['mapped_luca_pfams'] > 0 and r['alva_na_binding']:
            print(f"\n  fragment {r['fragment_id']}:")
            print(f"    Alva: {r['alva_na_binding']}-binding, "
                  f"ribosomal={r['ribosomal'] or 'none'}, "
                  f"{r['num_folds']} folds / {r['num_superfamilies']} superfamilies")
            print(f"    LUCA matches ({r['mapped_luca_pfams']}): "
                  f"{r['mapped_rna']} RNA + {r['mapped_dna']} DNA")
            print(f"    domains: {r['mapped_domains']}")

    # how many unique LUCA Pfams are reached by any fragment?
    all_reached = set()
    for r in results:
        for part in r['mapped_domains'].split('; '):
            if part and part.startswith('PF'):
                pid = part.split('(')[0]
                all_reached.add(pid)

    print(f"\n{'=' * 70}")
    print("COVERAGE: HOW MUCH OF LUCA'S NA-BINDING PROTEOME IS 'ANCIENT VOCABULARY'?")
    print("=" * 70)
    print(f"  unique LUCA NA-binding Pfams reached by Alva fragments: "
          f"{len(all_reached)} / {len(func_classes)} ({len(all_reached)/len(func_classes)*100:.0f}%)")

    # breakdown by functional category
    reached_by_function = Counter()
    total_by_function = Counter()
    for pid, info in func_classes.items():
        fn = info['function']
        total_by_function[fn] += 1
        if pid in all_reached:
            reached_by_function[fn] += 1

    print(f"\n  {'function':25s} {'reached':>8s} {'total':>6s} {'coverage':>9s}")
    print(f"  {'-'*25} {'-'*8} {'-'*6} {'-'*9}")
    for fn in sorted(total_by_function.keys()):
        r = reached_by_function.get(fn, 0)
        t = total_by_function[fn]
        pct = r / t * 100 if t > 0 else 0
        print(f"  {fn:25s} {r:8d} {t:6d} {pct:8.0f}%")

    # save summary JSON
    summary = {
        'total_fragments': len(results),
        'fragments_mapped_to_luca': len(mapped_any),
        'fragments_mapped_to_rna': len(mapped_rna),
        'fragments_mapped_to_dna': len(mapped_dna),
        'fragments_na_annotated': len(na_annotated),
        'unique_luca_pfams_reached': len(all_reached),
        'total_luca_na_pfams': len(func_classes),
        'coverage_fraction': round(len(all_reached) / len(func_classes), 3),
        'key_finding': (
            f"{len(all_reached)} of {len(func_classes)} LUCA NA-binding Pfam domains "
            f"({len(all_reached)/len(func_classes)*100:.0f}%) contain structural "
            f"elements from Alva's ancient peptide vocabulary"
        ),
    }
    json_path = os.path.join(PROJECT_DIR, 'results', 'phase1_alva_summary.json')
    with open(json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\nwrote summary to {json_path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
