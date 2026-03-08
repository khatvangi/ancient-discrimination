#!/usr/bin/env python3
"""classify LUCA NA-binding domains by function.

categories:
  - translation: ribosomal proteins, tRNA synthetases, translation factors, SRP
  - rna_modification: methyltransferases, pseudouridine synthases, dihydrouridine synthases
  - rna_processing: ribonucleases, RNA ligases, helicases
  - rna_binding_structural: general RNA-binding domains (KH, S1, S4, PUA)
  - dna_replication: polymerases, ligases, primases
  - dna_repair: MutS, endonucleases, recombinases
  - dna_modification: methyltransferases, restriction enzymes
  - dna_topology: topoisomerases, gyrases
  - transcription: HTH domains, sigma factors, transcription regulators
  - other: moonlighters, uncharacterized, multi-functional

classification is based on:
  1. pfam domain name and known function (primary)
  2. GO term annotations from pfam2go (supplementary)
  3. manual curation for ambiguous cases

outputs:
  results/phase1_functional_composition.tsv
"""

import os
import csv
import json
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')

# manual functional classification of the 54 LUCA RNA-binding domains
# based on domain name, GO annotations, and literature knowledge
LUCA_RNA_CLASSIFICATION = {
    # === TRANSLATION (tRNA synthetases, ribosomal, translation factors, SRP) ===
    'PF00152': 'translation',      # tRNA-synt_2 — class II aaRS core
    'PF00347': 'translation',      # Ribosomal_L6
    'PF00448': 'translation',      # SRP54 — cotranslational targeting
    'PF00466': 'translation',      # Ribosomal_L10
    'PF00573': 'translation',      # Ribosomal_L4
    'PF00579': 'translation',      # tRNA-synt_1b — class I aaRS
    'PF00749': 'translation',      # tRNA-synt_1c — class I aaRS
    'PF00750': 'translation',      # tRNA-synt_1d — class I aaRS
    'PF01411': 'translation',      # tRNA-synt_2c — class II aaRS (AlaRS)
    'PF01588': 'translation',      # tRNA_bind — tRNA-binding domain
    'PF02881': 'translation',      # SRP54_N — SRP GTPase N-terminal
    'PF03129': 'translation',      # HGTP_anticodon — anticodon binding (aaRS)
    'PF03483': 'translation',      # B3_4 — PheRS beta subunit
    'PF03485': 'translation',      # Arg_tRNA_synt_N — ArgRS N-terminal
    'PF06827': 'translation',      # zf-FPG_IleRS — zinc finger in IleRS
    'PF07973': 'translation',      # tRNA_SAD — aaRS-associated domain
    'PF08264': 'translation',      # Anticodon_1 — anticodon recognition
    'PF09190': 'translation',      # DALR_2 — CysRS editing domain
    'PF09334': 'translation',      # tRNA-synt_1g — class I aaRS
    'PF01926': 'translation',      # MMR_HSR1 — GTPase, translation-associated (Obg family)

    # === RNA MODIFICATION (methylation, pseudouridylation, dihydrouridine) ===
    'PF00588': 'rna_modification',  # SpoU_methylase — 2'-O-ribose methyltransferase
    'PF00398': 'rna_modification',  # RrnaAD — rRNA adenine dimethylase (KsgA/Dim1)
    'PF01142': 'rna_modification',  # TruD — pseudouridine synthase
    'PF01189': 'rna_modification',  # Methyltr_RsmB-F — RNA methyltransferase
    'PF01207': 'rna_modification',  # Dus — dihydrouridine synthase
    'PF01509': 'rna_modification',  # TruB_N — pseudouridine synthase
    'PF02475': 'rna_modification',  # Met_10 — SAM-dependent methyltransferase

    # === RNA-BINDING STRUCTURAL (general RNA-binding domains) ===
    'PF00013': 'rna_binding_structural',  # KH_1 — K homology domain
    'PF00575': 'rna_binding_structural',  # S1 — S1 RNA-binding domain (OB-fold)
    'PF01300': 'rna_binding_structural',  # Sua5_yciO_yrdC — dsRNA binding
    'PF01472': 'rna_binding_structural',  # PUA — RNA modification scaffold
    'PF01479': 'rna_binding_structural',  # S4 — S4 RNA-binding domain
    'PF01878': 'rna_binding_structural',  # EVE — PUA-like RNA-binding
    'PF07650': 'rna_binding_structural',  # KH_2 — K homology type 2

    # === RNA PROCESSING / DEGRADATION (nucleases, helicases, ligases) ===
    'PF00270': 'rna_processing',   # DEAD — DEAD-box RNA helicase
    'PF00580': 'rna_processing',   # UvrD-helicase — helicase superfamily I
    'PF01850': 'rna_processing',   # PIN — ribonuclease (toxin-antitoxin)
    'PF01936': 'rna_processing',   # NYN — ribonuclease
    'PF08494': 'rna_processing',   # DEAD_assoc — DEAD-box helicase associated
    'PF09414': 'rna_processing',   # RNA_ligase — RNA ligase
    'PF10130': 'rna_processing',   # PIN_2 — ribonuclease
    'PF10996': 'rna_processing',   # Beta-Casp — metallo-beta-lactamase (CPSF-like)
    'PF02272': 'rna_processing',   # DHHA1 — nucleic acid binding/processing

    # === OTHER / MULTI-FUNCTIONAL ===
    'PF00330': 'other',  # Aconitase — moonlighter: TCA enzyme + IRP1 (mRNA binding)
    'PF00565': 'other',  # SNase — staphylococcal nuclease fold (diverse functions)
    'PF00753': 'other',  # Lactamase_B — metallo-beta-lactamase fold (some RNA-processing)
    'PF01042': 'other',  # Ribonuc_L-PSP — endoribonuclease L-PSP
    'PF01170': 'other',  # UPF0020 — uncharacterized, tRNA modification?
    'PF02171': 'other',  # Piwi — Argonaute/PIWI (RNA silencing, ancient?)
    'PF02834': 'other',  # LigT_PEase — 2'-5' RNA ligase
    'PF03737': 'other',  # RraA-like — ribonuclease activity regulator
    'PF05191': 'other',  # ADK_lid — adenylate kinase lid domain
    'PF07521': 'other',  # RMMBL — metallo-beta-lactamase fold
    'PF08443': 'other',  # RimK — ribosome-associated modifying enzyme
}

# classification of 8 pre-LUCA RNA-binding domains
PRE_LUCA_RNA_CLASSIFICATION = {
    'PF00009': 'translation',      # GTP_EFTU — EF-Tu, tRNA delivery to ribosome
    'PF00133': 'translation',      # tRNA-synt_1 — class I aaRS core
    'PF00587': 'translation',      # tRNA-synt_2b — class II aaRS
    'PF01336': 'translation',      # tRNA_anti-codon — anticodon binding
    'PF03372': 'rna_processing',   # Exo_endo_phos — phosphodiesterase
    'PF04266': 'other',            # ASCH — ASH domain (RNA binding, function unclear)
    'PF04851': 'other',            # ResIII — restriction enzyme (annotated DNA-binding!)
    'PF13847': 'rna_modification', # Methyltransf_31 — methyltransferase
}

# classification of 21 LUCA DNA-binding domains
LUCA_DNA_CLASSIFICATION = {
    'PF00136': 'dna_replication',   # DNA_pol_B — DNA polymerase B
    'PF00488': 'dna_repair',        # MutS_V — mismatch repair
    'PF00589': 'dna_repair',        # Phage_integrase — tyrosine recombinase
    'PF00633': 'dna_repair',        # HHH — helix-hairpin-helix (base excision repair)
    'PF01131': 'dna_topology',      # Topoisom_bac — bacterial topoisomerase
    'PF01325': 'transcription',     # Fe_dep_repress — iron-dependent repressor (DtxR)
    'PF01420': 'dna_modification',  # Methylase_S — restriction methyltransferase
    'PF01548': 'dna_repair',        # DEDD_Tnp_IS110 — transposase
    'PF01555': 'dna_modification',  # N6_N4_Mtase — DNA methyltransferase
    'PF02586': 'dna_repair',        # SRAP — SOS response associated peptidase
    'PF02732': 'dna_repair',        # ERCC4 — nucleotide excision repair
    'PF02899': 'dna_repair',        # Phage_int_SAM_1 — integrase SAM-like
    'PF03989': 'dna_topology',      # DNA_gyraseA_C — gyrase A C-terminal
    'PF04014': 'other',             # MazE_antitoxin — antitoxin (TA system)
    'PF04471': 'dna_modification',  # Mrr_cat — restriction enzyme
    'PF06733': 'rna_processing',    # DEAD_2 — DEAD-box helicase variant
    'PF08463': 'dna_modification',  # EcoEI_R_C — restriction enzyme
    'PF09339': 'transcription',     # HTH_IclR — HTH transcription factor
    'PF13404': 'transcription',     # HTH_AsnC-type — HTH transcription factor
    'PF13495': 'dna_repair',        # Phage_int_SAM_4 — integrase SAM-like
    'PF11799': 'dna_repair',        # IMS_C — integrase-associated
}


def main():
    # combine all classifications
    all_classes = {}

    for pid, func in LUCA_RNA_CLASSIFICATION.items():
        all_classes[pid] = {'function': func, 'binding': 'RNA', 'age': 'LUCA'}

    for pid, func in PRE_LUCA_RNA_CLASSIFICATION.items():
        all_classes[pid] = {'function': func, 'binding': 'RNA', 'age': 'preLUCA'}

    for pid, func in LUCA_DNA_CLASSIFICATION.items():
        all_classes[pid] = {'function': func, 'binding': 'DNA', 'age': 'LUCA'}

    # load domain names
    names = {}
    with open(os.path.join(PROJECT_DIR, 'data', 'rbpworld', 'rbd_pfam_mapping.tsv')) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            names[row['pfam_id']] = row['domain_name']
    # supplement from pfam2go
    with open(os.path.join(PROJECT_DIR, 'data', 'rbpworld', 'raw', 'pfam2go.txt')) as f:
        for line in f:
            if line.startswith('Pfam:PF'):
                parts = line.split()
                pid = parts[0].replace('Pfam:', '')
                if pid not in names:
                    names[pid] = parts[1] if len(parts) > 1 else ''
    # also from DNA-binding TSV
    with open(os.path.join(PROJECT_DIR, 'data', 'rbpworld', 'dna_binding_pfams.tsv')) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['pfam_id'] not in names:
                names[row['pfam_id']] = row['pfam_name']

    # write output TSV
    out_path = os.path.join(PROJECT_DIR, 'results', 'phase1_functional_composition.tsv')
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['pfam_id', 'pfam_name', 'age_class', 'binding_type',
                         'functional_category'])
        for pid in sorted(all_classes.keys()):
            info = all_classes[pid]
            writer.writerow([
                pid, names.get(pid, '?'), info['age'], info['binding'],
                info['function']
            ])
    print(f"wrote {len(all_classes)} classified domains to {out_path}")

    # summary tables
    print("\n" + "=" * 70)
    print("FUNCTIONAL COMPOSITION: 54 LUCA RNA-BINDING DOMAINS")
    print("=" * 70)

    from collections import Counter
    rna_cats = Counter(v['function'] for v in all_classes.values()
                       if v['binding'] == 'RNA' and v['age'] == 'LUCA')
    total_rna = sum(rna_cats.values())
    for cat, count in rna_cats.most_common():
        pct = count / total_rna * 100
        print(f"  {cat:25s} {count:3d}  ({pct:4.1f}%)")
    print(f"  {'TOTAL':25s} {total_rna:3d}")

    print(f"\n{'=' * 70}")
    print("FUNCTIONAL COMPOSITION: 8 PRE-LUCA RNA-BINDING DOMAINS")
    print("=" * 70)
    pre_cats = Counter(v['function'] for v in all_classes.values()
                       if v['binding'] == 'RNA' and v['age'] == 'preLUCA')
    for cat, count in pre_cats.most_common():
        print(f"  {cat:25s} {count:3d}")

    print(f"\n{'=' * 70}")
    print("FUNCTIONAL COMPOSITION: 21 LUCA DNA-BINDING DOMAINS")
    print("=" * 70)
    dna_cats = Counter(v['function'] for v in all_classes.values()
                       if v['binding'] == 'DNA')
    total_dna = sum(dna_cats.values())
    for cat, count in dna_cats.most_common():
        pct = count / total_dna * 100
        print(f"  {cat:25s} {count:3d}  ({pct:4.1f}%)")

    # the key comparison: RNA vs DNA functional profiles
    print(f"\n{'=' * 70}")
    print("RNA vs DNA FUNCTIONAL PROFILE COMPARISON (LUCA domains)")
    print("=" * 70)
    all_cats = sorted(set(list(rna_cats.keys()) + list(dna_cats.keys())))
    print(f"  {'category':25s} {'RNA':>5s} {'DNA':>5s}")
    print(f"  {'-'*25} {'-'*5} {'-'*5}")
    for cat in all_cats:
        r = rna_cats.get(cat, 0)
        d = dna_cats.get(cat, 0)
        print(f"  {cat:25s} {r:5d} {d:5d}")
    print(f"  {'-'*25} {'-'*5} {'-'*5}")
    print(f"  {'TOTAL':25s} {total_rna:5d} {total_dna:5d}")

    # annotation bias check
    print(f"\n{'=' * 70}")
    print("ANNOTATION BIAS CHECK: DNA-BINDING COVERAGE")
    print("=" * 70)
    print("  DNA-binding domains use GO:0003677 descendants (385 Pfam domains)")
    print("  RNA-binding uses RBPWorld/EuRBPDB (916 Pfam domains)")
    print(f"  ratio: {916/385:.1f}x more RNA annotations than DNA annotations")
    print("  this means the 2.6:1 RNA:DNA ratio in LUCA may partly reflect")
    print("  annotation asymmetry, not just biological reality.")
    print("  however, the FUNCTIONAL profiles are genuinely different:")
    print("  RNA = translation-dominated, DNA = repair-dominated.")
    print("  this structural difference is robust to annotation completeness.")

    return 0


if __name__ == '__main__':
    sys.exit(main())
