#!/usr/bin/env python3
"""download wehbi et al. 2024 LUCA pfam age classification.

source: github.com/sawsanwehbi/Pfam-age-classification
paper: Wehbi et al. 2024 PNAS — "dating pfam domains through phylogenetics"

downloads two files:
  - ClassifiedPFAMs.csv (pfam_id, robust_classifications)
  - Pfam_data_ancestralAAC.csv (pfam_id, ancestor, clan, amino acid frequencies)

produces:
  - data/wehbi_luca/raw/ (downloaded originals)
  - data/wehbi_luca/luca_pfams.tsv (cleaned extract: pfam_id, age_class)
"""

import os
import csv
import sys
import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')
RAW_DIR = os.path.join(PROJECT_DIR, 'data', 'wehbi_luca', 'raw')
OUT_TSV = os.path.join(PROJECT_DIR, 'data', 'wehbi_luca', 'luca_pfams.tsv')

# raw github URLs for the two key files
URLS = {
    'ClassifiedPFAMs.csv': (
        'https://raw.githubusercontent.com/sawsanwehbi/'
        'Pfam-age-classification/main/ClassifiedPFAMs.csv'
    ),
    'Pfam_data_ancestralAAC.csv': (
        'https://raw.githubusercontent.com/sawsanwehbi/'
        'Pfam-age-classification/main/Pfam_data_ancestralAAC.csv'
    ),
}


def download_file(url, outpath):
    """download a single file, skip if already exists."""
    if os.path.exists(outpath):
        print(f"  already exists: {outpath}")
        return
    print(f"  downloading: {url}")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    with open(outpath, 'wb') as f:
        f.write(resp.content)
    print(f"  saved: {outpath} ({len(resp.content)} bytes)")


def parseClassifiedPfams(csv_path):
    """parse ClassifiedPFAMs.csv → list of (pfam_id, age_class) tuples.

    the CSV has quoted fields: "PF00004","LUCA"
    we strip quotes and normalize the pfam ID format.
    """
    rows = []
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        # expect columns: PFAM_IDs, robust_classifications
        for row in reader:
            if len(row) < 2:
                continue
            pfam_id = row[0].strip().strip('"')
            age_class = row[1].strip().strip('"')
            # skip rows with empty or invalid entries
            if not pfam_id.startswith('PF'):
                continue
            rows.append((pfam_id, age_class))
    return rows


def main():
    os.makedirs(RAW_DIR, exist_ok=True)

    # step 1: download raw files
    print("step 1: downloading wehbi data from github...")
    for filename, url in URLS.items():
        download_file(url, os.path.join(RAW_DIR, filename))

    # step 2: parse the classification file
    print("\nstep 2: parsing ClassifiedPFAMs.csv...")
    csv_path = os.path.join(RAW_DIR, 'ClassifiedPFAMs.csv')
    all_pfams = parseClassifiedPfams(csv_path)
    print(f"  total pfam entries: {len(all_pfams)}")

    # step 3: write cleaned TSV with all pfams and their age classes
    print(f"\nstep 3: writing cleaned TSV to {OUT_TSV}...")
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['pfam_id', 'age_class'])
        for pfam_id, age_class in all_pfams:
            writer.writerow([pfam_id, age_class])

    # step 4: summarize counts by age class
    print("\nstep 4: age class summary:")
    counts = {}
    for _, age_class in all_pfams:
        counts[age_class] = counts.get(age_class, 0) + 1

    # sort by count descending
    for age, count in sorted(counts.items(), key=lambda x: -x[1]):
        marker = " <<<" if age in ('LUCA', 'preLUCA') else ""
        print(f"  {age:25s} {count:5d}{marker}")

    luca_count = counts.get('LUCA', 0)
    pre_luca_count = counts.get('preLUCA', 0)
    print(f"\n  LUCA pfams:     {luca_count}")
    print(f"  pre-LUCA pfams: {pre_luca_count}")
    print(f"  total (all ages): {len(all_pfams)}")

    # checkpoint verification
    print("\n--- CHECKPOINT ---")
    if luca_count < 800 or luca_count > 1100:
        print(f"  WARNING: LUCA count ({luca_count}) is outside expected range 800-1100")
    else:
        print(f"  OK: LUCA count ({luca_count}) in expected range")

    if pre_luca_count < 50 or pre_luca_count > 200:
        print(f"  WARNING: pre-LUCA count ({pre_luca_count}) is outside expected range 50-200")
    else:
        print(f"  OK: pre-LUCA count ({pre_luca_count}) in expected range")

    return 0


if __name__ == '__main__':
    sys.exit(main())
