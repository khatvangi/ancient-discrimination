#!/usr/bin/env python3
"""
query_pfam_crossdomain.py
for each priority Pfam domain, look up:
1. InterPro ID mapping
2. bacteria and archaea representation from universality audit
3. or query InterPro API for taxonomy distribution
"""

import requests
import json
import time
import sys
import csv

# the dual-binder and key domains we found have PDB structures
DOMAINS_WITH_STRUCTURES = [
    "PF00013", "PF00575", "PF01479", "PF01472",
    "PF00633", "PF02171",
    "PF00009", "PF00133", "PF00152",
    "PF00347", "PF00573", "PF00466",
    "PF00749", "PF00750",
    "PF00270", "PF00448", "PF00580",
    "PF00587", "PF00488", "PF00136",
    "PF01555",
    "PF00398", "PF00579", "PF00588", "PF00589",
    "PF01131", "PF01336", "PF01509", "PF01926",
    "PF02272", "PF03372",
]


def findInterProMapping(pfam_id):
    """look up the InterPro ID for a given Pfam domain"""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}"
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        meta = data.get("metadata", {})
        name_info = meta.get("name", {})
        short_name = name_info.get("short", "unknown") if isinstance(name_info, dict) else str(name_info)
        full_name = name_info.get("name", "unknown") if isinstance(name_info, dict) else str(name_info)
        ipr = meta.get("integrated", None)
        return {
            "pfam_id": pfam_id,
            "short_name": short_name,
            "full_name": full_name,
            "interpro_id": ipr
        }
    except Exception as e:
        print(f"  error looking up {pfam_id}: {e}", file=sys.stderr)
        return {
            "pfam_id": pfam_id,
            "short_name": "unknown",
            "full_name": "unknown",
            "interpro_id": None
        }


def main():
    # load the universality audit for cross-referencing
    audit_path = "/storage/kiran-stuff/Vocabulary-Proteins/canon_union_universality_audit.tsv"
    audit = {}
    try:
        with open(audit_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                ipr = row["IPR"]
                audit[ipr] = {
                    "bacteria": int(row["Bacteria"]),
                    "archaea": int(row["Archaea"]),
                    "eukaryota": int(row["Eukaryota"]),
                }
    except Exception as e:
        print(f"warning: could not load audit file: {e}", file=sys.stderr)

    results = []
    for pfam_id in DOMAINS_WITH_STRUCTURES:
        print(f"looking up {pfam_id}...", file=sys.stderr)
        info = findInterProMapping(pfam_id)
        time.sleep(0.3)

        # check audit data
        ipr = info["interpro_id"]
        if ipr and ipr in audit:
            info["bacteria"] = audit[ipr]["bacteria"]
            info["archaea"] = audit[ipr]["archaea"]
            info["eukaryota"] = audit[ipr]["eukaryota"]
            info["in_audit"] = "yes"
        else:
            info["bacteria"] = "?"
            info["archaea"] = "?"
            info["eukaryota"] = "?"
            info["in_audit"] = "no"

        results.append(info)

    # output TSV
    header = ["pfam_id", "short_name", "full_name", "interpro_id",
              "in_audit", "bacteria", "archaea", "eukaryota"]
    print("\t".join(header))
    for r in results:
        print("\t".join(str(r[h]) for h in header))


if __name__ == "__main__":
    main()
