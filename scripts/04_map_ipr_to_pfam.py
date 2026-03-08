#!/usr/bin/env python3
"""map InterPro IDs from green_families.txt to Pfam IDs.

uses InterPro REST API: /api/entry/interpro/IPRXXXXXX
extracts member_databases.pfam entries.

important: one IPR can map to multiple Pfam domains, and one Pfam domain
can appear in multiple IPR entries. downstream intersections must join
on Pfam accession, not InterPro ID.

also maps the 415 universal families for broader context.

outputs:
  - data/ipr_to_pfam_map.tsv (ipr_id, pfam_id, pfam_name, ipr_name, source)
"""

import os
import csv
import sys
import json
import time
import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')
PARENT_DIR = os.path.join(PROJECT_DIR, '..')  # Vocabulary-Proteins root

GREEN_FILE = os.path.join(PARENT_DIR, 'green_families.txt')
UNIVERSAL_FILE = os.path.join(PARENT_DIR, 'canon_union_universal_na_binding_families.txt')
OUT_TSV = os.path.join(PROJECT_DIR, 'data', 'ipr_to_pfam_map.tsv')

INTERPRO_API = 'https://www.ebi.ac.uk/interpro/api/entry/interpro/{}'


def loadIprIds(filepath, label):
    """load InterPro IDs from a file (one per line)."""
    ids = []
    with open(filepath) as f:
        for line in f:
            ipr = line.strip()
            if ipr.startswith('IPR'):
                ids.append(ipr)
    print(f"  loaded {len(ids)} InterPro IDs from {label}")
    return ids


def queryInterproPfamMapping(ipr_id, max_retries=3):
    """query InterPro API for a single entry, extract Pfam cross-references.

    returns: list of (pfam_id, pfam_name) tuples, plus ipr_name
    """
    url = INTERPRO_API.format(ipr_id)
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, headers={'Accept': 'application/json'}, timeout=30)
            if resp.status_code == 404:
                return [], f"NOT_FOUND({ipr_id})"
            resp.raise_for_status()
            data = resp.json()

            meta = data.get('metadata', {})
            ipr_name = meta.get('name', {})
            if isinstance(ipr_name, dict):
                ipr_name = ipr_name.get('name', '') or ipr_name.get('short', '')

            # extract pfam member databases
            member_dbs = meta.get('member_databases', {})
            pfam_entries = member_dbs.get('pfam', {})

            results = []
            for pfam_id, pfam_info in pfam_entries.items():
                pfam_name = pfam_info if isinstance(pfam_info, str) else str(pfam_info)
                results.append((pfam_id, pfam_name))

            return results, ipr_name

        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                return [], f"ERROR({e})"

    return [], "MAX_RETRIES"


def main():
    # step 1: load IPR IDs
    print("step 1: loading InterPro IDs...")
    green_iprs = loadIprIds(GREEN_FILE, "GREEN-51")
    universal_iprs = loadIprIds(UNIVERSAL_FILE, "universal-415")

    # combine, tracking source
    all_iprs = {}
    for ipr in green_iprs:
        all_iprs[ipr] = 'green51'
    for ipr in universal_iprs:
        if ipr in all_iprs:
            all_iprs[ipr] = 'green51+universal'
        else:
            all_iprs[ipr] = 'universal_only'

    print(f"  total unique IPR IDs: {len(all_iprs)}")
    print(f"  green51-only: {sum(1 for v in all_iprs.values() if v == 'green51')}")
    print(f"  in both: {sum(1 for v in all_iprs.values() if v == 'green51+universal')}")
    print(f"  universal-only: {sum(1 for v in all_iprs.values() if v == 'universal_only')}")

    # step 2: query InterPro API for each IPR ID
    print("\nstep 2: querying InterPro API for Pfam mappings...")
    rows = []
    unmapped = []
    errors = []

    for i, (ipr_id, source) in enumerate(sorted(all_iprs.items())):
        pfam_list, ipr_name = queryInterproPfamMapping(ipr_id)

        if not pfam_list:
            if 'NOT_FOUND' in str(ipr_name) or 'ERROR' in str(ipr_name):
                errors.append((ipr_id, ipr_name))
                print(f"  [{i+1}/{len(all_iprs)}] {ipr_id}: {ipr_name}")
            else:
                unmapped.append((ipr_id, ipr_name))
                print(f"  [{i+1}/{len(all_iprs)}] {ipr_id}: no Pfam mapping ({ipr_name})")
        else:
            for pfam_id, pfam_name in pfam_list:
                rows.append((ipr_id, pfam_id, pfam_name, ipr_name, source))
            if len(pfam_list) > 1:
                print(f"  [{i+1}/{len(all_iprs)}] {ipr_id}: {len(pfam_list)} Pfam domains ({ipr_name})")

        # rate limit: 0.15s between requests
        if (i + 1) % 10 == 0:
            print(f"  ... {i+1}/{len(all_iprs)} queried")
        time.sleep(0.15)

    # step 3: write output
    print(f"\nstep 3: writing {len(rows)} mappings to {OUT_TSV}...")
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)
    with open(OUT_TSV, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['ipr_id', 'pfam_id', 'pfam_name', 'ipr_name', 'source'])
        for row in rows:
            writer.writerow(row)

    # step 4: summary
    print("\n--- SUMMARY ---")
    unique_pfams = set(r[1] for r in rows)
    green_pfams = set(r[1] for r in rows if 'green51' in r[4])
    universal_pfams = set(r[1] for r in rows if 'universal' in r[4])

    print(f"  total IPR→Pfam mappings: {len(rows)}")
    print(f"  unique Pfam IDs: {len(unique_pfams)}")
    print(f"  Pfam IDs from GREEN-51: {len(green_pfams)}")
    print(f"  Pfam IDs from universal-415: {len(universal_pfams)}")
    print(f"  IPRs with no Pfam mapping: {len(unmapped)}")
    print(f"  IPRs with errors: {len(errors)}")

    if unmapped:
        print(f"\n  unmapped IPRs:")
        for ipr_id, name in unmapped:
            print(f"    {ipr_id}: {name}")

    if errors:
        print(f"\n  errored IPRs:")
        for ipr_id, err in errors:
            print(f"    {ipr_id}: {err}")

    # check directionality: any pfam in multiple IPRs?
    pfam_to_iprs = {}
    for row in rows:
        pfam_to_iprs.setdefault(row[1], []).append(row[0])
    multi_ipr_pfams = {p: iprs for p, iprs in pfam_to_iprs.items() if len(iprs) > 1}
    if multi_ipr_pfams:
        print(f"\n  Pfam IDs appearing in multiple IPR entries: {len(multi_ipr_pfams)}")
        for pfam_id, iprs in sorted(multi_ipr_pfams.items())[:5]:
            print(f"    {pfam_id} → {', '.join(iprs)}")

    # checkpoint
    print("\n--- CHECKPOINT ---")
    if len(green_pfams) < 20:
        print(f"  WARNING: only {len(green_pfams)} unique Pfam IDs from GREEN-51 (expected ~40-60)")
    else:
        print(f"  OK: {len(green_pfams)} unique Pfam IDs from GREEN-51")

    return 0


if __name__ == '__main__':
    sys.exit(main())
