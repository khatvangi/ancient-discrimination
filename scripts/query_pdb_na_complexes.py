#!/usr/bin/env python3
"""
query_pdb_na_complexes.py
query RCSB PDB for protein-nucleic acid co-crystal structures for LUCA-age Pfam domains.
reads families from phase2_selected_families.tsv, queries the RCSB Search API for
structures with each Pfam domain + RNA, DNA, or both.

outputs: data/pdb_na_cocrystals.tsv

usage:
    python scripts/query_pdb_na_complexes.py
"""

import requests
import json
import time
import sys
import csv
import os

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"

# input/output paths (relative to project root)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_FILE = os.path.join(PROJECT_ROOT, "results", "phase2_selected_families.tsv")
OUTPUT_FILE = os.path.join(PROJECT_ROOT, "data", "pdb_na_cocrystals.tsv")


def buildQuery(pfam_id, na_type="RNA"):
    """
    build RCSB Search API query for structures containing a Pfam domain
    AND nucleic acid polymer(s).

    na_type: "RNA", "DNA", or "both" (RNA AND DNA in same structure)
    returns the JSON query body.
    """
    # pfam domain filter
    pfam_node = {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_polymer_entity_annotation.annotation_id",
            "operator": "exact_match",
            "value": pfam_id
        }
    }

    if na_type == "both":
        # require both RNA and DNA polymer entities
        rna_node = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_RNA",
                "operator": "greater",
                "value": 0
            }
        }
        dna_node = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_DNA",
                "operator": "greater",
                "value": 0
            }
        }
        nodes = [pfam_node, rna_node, dna_node]
    else:
        # single NA type
        attr = ("rcsb_entry_info.polymer_entity_count_RNA" if na_type == "RNA"
                else "rcsb_entry_info.polymer_entity_count_DNA")
        na_node = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": attr,
                "operator": "greater",
                "value": 0
            }
        }
        nodes = [pfam_node, na_node]

    query_body = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": False,
            "paginate": {
                "start": 0,
                "rows": 10  # fetch up to 10 example IDs
            }
        }
    }
    return query_body


def queryPDB(pfam_id, na_type="RNA"):
    """
    execute RCSB Search API query for a Pfam domain + nucleic acid type.
    returns (total_count, list_of_pdb_ids).
    """
    query_body = buildQuery(pfam_id, na_type)

    try:
        resp = requests.post(SEARCH_URL, json=query_body, timeout=30)
        if resp.status_code == 204:
            # 204 = no results
            return (0, [])
        resp.raise_for_status()
        data = resp.json()
        total = data.get("total_count", 0)
        ids = [r["identifier"] for r in data.get("result_set", [])]
        return (total, ids)
    except requests.exceptions.HTTPError as e:
        print(f"  HTTP error for {pfam_id} ({na_type}): {e}", file=sys.stderr)
        return (0, [])
    except Exception as e:
        print(f"  error for {pfam_id} ({na_type}): {e}", file=sys.stderr)
        return (0, [])


def loadFamilies(filepath):
    """
    read phase2_selected_families.tsv and return list of (pfam_id, pfam_name) tuples.
    """
    families = []
    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pid = row.get("pfam_id", "").strip()
            pname = row.get("pfam_name", "").strip()
            if pid:
                families.append((pid, pname))
    return families


def main():
    # load families
    if not os.path.exists(INPUT_FILE):
        print(f"ERROR: input file not found: {INPUT_FILE}", file=sys.stderr)
        sys.exit(1)

    families = loadFamilies(INPUT_FILE)
    print(f"loaded {len(families)} families from {INPUT_FILE}", file=sys.stderr)

    # ensure output directory exists
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)

    results = []

    for pfam_id, pfam_name in families:
        print(f"querying {pfam_id} ({pfam_name})...", file=sys.stderr)

        # query RNA complexes
        rna_count, rna_ids = queryPDB(pfam_id, "RNA")
        time.sleep(0.4)  # be polite to the API

        # query DNA complexes
        dna_count, dna_ids = queryPDB(pfam_id, "DNA")
        time.sleep(0.4)

        # query structures with BOTH RNA and DNA
        both_count, both_ids = queryPDB(pfam_id, "both")
        time.sleep(0.4)

        # take up to 5 example PDB IDs per category
        rna_examples = rna_ids[:5]
        dna_examples = dna_ids[:5]

        results.append({
            "pfam_id": pfam_id,
            "pfam_name": pfam_name,
            "n_rna": rna_count,
            "n_dna": dna_count,
            "n_both": both_count,
            "example_rna_pdbs": ",".join(rna_examples) if rna_examples else "none",
            "example_dna_pdbs": ",".join(dna_examples) if dna_examples else "none",
        })

        print(f"  RNA: {rna_count}, DNA: {dna_count}, both: {both_count}", file=sys.stderr)

    # write output TSV
    header = ["pfam_id", "pfam_name", "n_rna", "n_dna", "n_both",
              "example_rna_pdbs", "example_dna_pdbs"]
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"\nwrote {len(results)} rows to {OUTPUT_FILE}", file=sys.stderr)

    # also print summary to stdout
    print(f"\n{'='*80}")
    print(f"PDB nucleic acid co-crystal census for {len(results)} LUCA Pfam domains")
    print(f"{'='*80}")
    print(f"{'Pfam ID':<12} {'Name':<20} {'RNA':>6} {'DNA':>6} {'Both':>6}")
    print(f"{'-'*12} {'-'*20} {'-'*6} {'-'*6} {'-'*6}")
    for r in results:
        print(f"{r['pfam_id']:<12} {r['pfam_name']:<20} {r['n_rna']:>6} {r['n_dna']:>6} {r['n_both']:>6}")

    # totals
    total_rna = sum(r['n_rna'] for r in results)
    total_dna = sum(r['n_dna'] for r in results)
    total_both = sum(r['n_both'] for r in results)
    n_with_rna = sum(1 for r in results if r['n_rna'] > 0)
    n_with_dna = sum(1 for r in results if r['n_dna'] > 0)
    n_with_both = sum(1 for r in results if r['n_both'] > 0)
    n_dual = sum(1 for r in results if r['n_rna'] > 0 and r['n_dna'] > 0)

    print(f"\n{'='*80}")
    print(f"families with RNA structures:  {n_with_rna}/{len(results)}")
    print(f"families with DNA structures:  {n_with_dna}/{len(results)}")
    print(f"families with both RNA+DNA:    {n_with_both}/{len(results)}")
    print(f"families with dual evidence:   {n_dual}/{len(results)} (have both RNA-only and DNA-only structures)")
    print(f"total RNA structures:          {total_rna}")
    print(f"total DNA structures:          {total_dna}")
    print(f"total both structures:         {total_both}")


if __name__ == "__main__":
    main()
