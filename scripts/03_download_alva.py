#!/usr/bin/env python3
"""parse Alva et al. 2015 eLife ancient peptide fragment data.

source: Alva, Söding & Lupas (2015) "A vocabulary of ancient peptides at
the origin of folded proteins." eLife 4:e09410.
doi: 10.7554/eLife.09410
PMC: PMC4739770

the 40 primordial fragments are encoded from Table 1 of the paper.
also downloads the supplementary MSA/accession file (fig3-data1) for
future cross-referencing.

outputs:
  - data/alva_fragments/raw/ (downloaded supplement)
  - data/alva_fragments/fragments.tsv (structured table of 40 fragments)
"""

import os
import csv
import sys
import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')
RAW_DIR = os.path.join(PROJECT_DIR, 'data', 'alva_fragments', 'raw')
OUT_TSV = os.path.join(PROJECT_DIR, 'data', 'alva_fragments', 'fragments.tsv')

# supplementary data URL (MSA details for the 40 fragments)
SUPP_URL = (
    'https://elifesciences.org/download/'
    'aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMDk0MTAv'
    'ZWxpZmUtMDk0MTAtZmlnMy1kYXRhMS12Mi5kb2N4/'
    'elife-09410-fig3-data1-v2.docx'
)

# table 1 from Alva et al. 2015 — 40 primordial fragments
# columns: fragment_id, num_folds, num_superfamilies, ribosomal, na_binding,
#          metal_fes, ligands, ancient_folds
# "—" in the paper means no annotation
FRAGMENTS_TABLE1 = [
    (1,  14, 20, 'S19e', 'DNA', '', '', 'a.4'),
    (2,  8,  15, 'S13',  'DNA+RNA', '', '', 'a.60'),
    (3,  3,  3,  '',     'DNA', '', '', 'c.37.1.20'),
    (4,  2,  2,  'L29',  'RNA', '', '', 'a.2'),
    (5,  2,  2,  'S5',   'RNA', '', '', ''),
    (6,  2,  2,  'S3',   'DNA+RNA', '', '', ''),
    (7,  6,  8,  '',     'DNA', '', 'CTP,SAM,FMN', 'c.2,c.66'),
    (8,  10, 10, '',     '',    '', 'FAD,NAD,cofactors', 'c.2,c.66'),
    (9,  2,  2,  '',     'DNA', '', '', ''),
    (10, 3,  4,  'S4',   'RNA', '', '', ''),
    (11, 2,  2,  'TL5',  'RNA', '', '', ''),
    (12, 4,  8,  'L27',  'RNA', '', '', 'b.84.1'),
    (13, 5,  30, '',     'DNA', '', '', 'all_folds'),
    (14, 2,  2,  '',     'DNA', 'ZN', '', ''),
    (15, 5,  7,  '',     '',    '', 'FMN,AMP,FAD', 'd.37.1'),
    (16, 3,  3,  '',     '',    '', 'ATP,GTP,ADP', 'c.37'),
    (17, 3,  3,  '',     '',    'HEM', '', 'a.24'),
    (18, 2,  2,  '',     '',    'SF4', '', 'd.58.1,d.58'),
    (19, 3,  3,  '',     '',    'ZN', '', ''),
    (20, 2,  2,  '',     '',    'FES', '', ''),
    (21, 2,  2,  '',     '',    'CA', '', ''),
    (22, 2,  2,  '',     '',    'FEC,FE', '', ''),
    (23, 2,  2,  '',     '',    '', 'COA', ''),
    (24, 2,  2,  'L7/12', '',   '', '', ''),
    (25, 3,  3,  '',     '',    '', '', ''),
    (26, 2,  2,  '',     '',    '', '', ''),
    (27, 7,  12, '',     '',    '', '', 'a.24,a.118'),
    (28, 2,  2,  '',     '',    '', '', ''),
    (29, 2,  2,  '',     '',    '', '', ''),
    (30, 2,  2,  '',     '',    '', '', ''),
    (31, 2,  2,  '',     '',    '', '', ''),
    (32, 2,  2,  '',     '',    '', '', ''),
    (33, 2,  2,  '',     '',    '', '', ''),
    (34, 3,  4,  '',     '',    '', '', ''),
    (35, 2,  2,  '',     '',    '', '', 'd.58'),
    (36, 2,  2,  '',     '',    '', '', ''),
    (37, 2,  2,  '',     '',    '', '', ''),
    (38, 2,  2,  '',     '',    '', '', ''),
    (39, 2,  2,  '',     '',    '', '', ''),
    (40, 2,  2,  '',     '',    '', '', ''),
]

COLUMNS = [
    'fragment_id', 'num_folds', 'num_superfamilies', 'ribosomal',
    'na_binding', 'metal_fes', 'ligands', 'scop_folds'
]


def main():
    os.makedirs(RAW_DIR, exist_ok=True)

    # step 1: try to download supplementary docx
    print("step 1: downloading supplementary MSA data...")
    supp_path = os.path.join(RAW_DIR, 'elife-09410-fig3-data1-v2.docx')
    if os.path.exists(supp_path):
        print(f"  already exists: {supp_path}")
    else:
        try:
            resp = requests.get(SUPP_URL, timeout=60)
            resp.raise_for_status()
            with open(supp_path, 'wb') as f:
                f.write(resp.content)
            print(f"  saved: {supp_path} ({len(resp.content)} bytes)")
        except Exception as e:
            print(f"  WARNING: could not download supplement: {e}")
            print(f"  (this is non-critical — table data is encoded directly)")

    # step 2: write structured TSV from Table 1
    print("\nstep 2: writing fragments TSV...")
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(COLUMNS)
        for row in FRAGMENTS_TABLE1:
            writer.writerow(row)

    print(f"  wrote {len(FRAGMENTS_TABLE1)} fragments to {OUT_TSV}")

    # step 3: summary statistics
    print("\nstep 3: summary statistics...")
    na_binding = [r for r in FRAGMENTS_TABLE1 if r[4]]  # column 4 = na_binding
    rna_binding = [r for r in FRAGMENTS_TABLE1 if 'RNA' in r[4]]
    dna_binding = [r for r in FRAGMENTS_TABLE1 if 'DNA' in r[4]]
    dual_binding = [r for r in FRAGMENTS_TABLE1 if 'RNA' in r[4] and 'DNA' in r[4]]
    ribosomal = [r for r in FRAGMENTS_TABLE1 if r[3]]
    metal_fes = [r for r in FRAGMENTS_TABLE1 if r[5]]

    print(f"  total fragments: {len(FRAGMENTS_TABLE1)}")
    print(f"  NA-binding:      {len(na_binding)} ({100*len(na_binding)/len(FRAGMENTS_TABLE1):.1f}%)")
    print(f"    RNA-binding:   {len(rna_binding)}")
    print(f"    DNA-binding:   {len(dna_binding)}")
    print(f"    dual (RNA+DNA): {len(dual_binding)}")
    print(f"  ribosomal:       {len(ribosomal)} ({100*len(ribosomal)/len(FRAGMENTS_TABLE1):.1f}%)")
    print(f"  metal/FeS:       {len(metal_fes)} ({100*len(metal_fes)/len(FRAGMENTS_TABLE1):.1f}%)")

    # list NA-binding fragments
    print("\n  NA-binding fragments:")
    for r in na_binding:
        print(f"    fragment {r[0]:2d}: {r[4]:10s} (folds={r[1]}, ribosomal={r[3] or '-'})")

    # checkpoint
    print("\n--- CHECKPOINT ---")
    if len(FRAGMENTS_TABLE1) != 40:
        print(f"  ERROR: expected 40 fragments, got {len(FRAGMENTS_TABLE1)}")
    else:
        print(f"  OK: 40 fragments")

    expected_na = 13  # from thesis: "33% bind nucleic acids"
    if len(na_binding) != expected_na:
        print(f"  CHECK: expected {expected_na} NA-binding, got {len(na_binding)}")
    else:
        print(f"  OK: {len(na_binding)} NA-binding fragments ({len(na_binding)/40*100:.0f}%)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
