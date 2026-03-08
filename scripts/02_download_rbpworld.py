#!/usr/bin/env python3
"""build RNA-binding and DNA-binding pfam domain lists using GO annotations.

strategy:
  1. download pfam2go mapping from gene ontology consortium
  2. use quickgo API to get all descendant GO terms of:
     - GO:0003723 (RNA binding) → captures tRNA, rRNA, mRNA, snRNA binding etc.
     - GO:0003677 (DNA binding) → captures ssDNA, dsDNA, chromatin binding etc.
     - GO:0003676 (nucleic acid binding) → general NA binders
  3. filter pfam2go for pfam domains annotated with any of these terms
  4. also extract EuRBPDB 791 RBD names from their HMM archive for cross-validation

sources:
  - pfam2go: current.geneontology.org/ontology/external2go/pfam2go
  - quickgo: ebi.ac.uk/QuickGO/services/ontology/go/terms/
  - eurbpdb: eurbpdb.gzsys.org.cn (Liao et al. 2020/2025)

outputs:
  - data/rbpworld/raw/pfam2go.txt (downloaded original)
  - data/rbpworld/rna_binding_pfams.tsv (pfam_id, pfam_name, go_terms, source)
  - data/rbpworld/dna_binding_pfams.tsv
  - data/rbpworld/na_binding_pfams.tsv (union: RNA or DNA or general NA)
"""

import os
import csv
import sys
import json
import time
import zipfile
import tempfile
import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')
RAW_DIR = os.path.join(PROJECT_DIR, 'data', 'rbpworld', 'raw')
OUT_DIR = os.path.join(PROJECT_DIR, 'data', 'rbpworld')

PFAM2GO_URL = 'https://current.geneontology.org/ontology/external2go/pfam2go'
QUICKGO_URL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}/descendants?relations=is_a'
EURBPDB_URL = 'http://eurbpdb.gzsys.org.cn/data/download/791_RBDs.PFam.gz'

# root GO terms for binding categories
RNA_BINDING_GO = 'GO:0003723'
DNA_BINDING_GO = 'GO:0003677'
NA_BINDING_GO = 'GO:0003676'   # general nucleic acid binding


def download_file(url, outpath, description="file"):
    """download a single file, skip if already exists."""
    if os.path.exists(outpath):
        print(f"  already exists: {outpath}")
        return True
    print(f"  downloading {description}: {url}")
    try:
        resp = requests.get(url, timeout=120)
        resp.raise_for_status()
        with open(outpath, 'wb') as f:
            f.write(resp.content)
        print(f"  saved: {outpath} ({len(resp.content)} bytes)")
        return True
    except Exception as e:
        print(f"  FAILED: {e}")
        return False


def fetchGoDescendants(go_term):
    """fetch all descendant GO terms (is_a relationship) from quickgo API."""
    url = QUICKGO_URL.format(go_term)
    print(f"  querying quickgo for descendants of {go_term}...")
    try:
        resp = requests.get(url, headers={'Accept': 'application/json'}, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        terms = set()
        for result in data.get('results', []):
            for d in result.get('descendants', []):
                terms.add(d)
        # always include the root term itself
        terms.add(go_term)
        print(f"  found {len(terms)} terms (including {go_term} itself)")
        return terms
    except Exception as e:
        print(f"  WARNING: quickgo query failed ({e}), using root term only")
        return {go_term}


def parsePfam2go(filepath):
    """parse pfam2go file → dict of pfam_id → [(go_term, go_name, pfam_name), ...]"""
    pfam_go = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('!') or not line:
                continue
            # format: Pfam:PF00013 KH_1 > GO:RNA binding ; GO:0003723
            try:
                left, right = line.split(' > ')
                pfam_part = left.split()
                pfam_id = pfam_part[0].replace('Pfam:', '')
                pfam_name = ' '.join(pfam_part[1:])

                # parse GO term
                go_name_part, go_id = right.rsplit(' ; ', 1)
                go_name = go_name_part.replace('GO:', '')
                go_id = go_id.strip()

                if pfam_id not in pfam_go:
                    pfam_go[pfam_id] = {'name': pfam_name, 'go_terms': []}
                pfam_go[pfam_id]['go_terms'].append((go_id, go_name))
            except (ValueError, IndexError):
                continue
    return pfam_go


def extractEurbpdbDomains(zip_path):
    """extract domain names from EuRBPDB HMM zip archive."""
    names = set()
    try:
        with zipfile.ZipFile(zip_path, 'r') as zf:
            for name in zf.namelist():
                # format: 791DomainPfaFile/DOMAIN_NAME.hmm
                if name.endswith('.hmm'):
                    domain = os.path.basename(name).replace('.hmm', '')
                    names.add(domain)
    except Exception as e:
        print(f"  WARNING: could not read EuRBPDB archive: {e}")
    return names


def findPfamsByGoTerms(pfam_go, target_go_terms):
    """find all pfam domains annotated with any of the target GO terms."""
    matches = {}
    for pfam_id, info in pfam_go.items():
        matching_gos = []
        for go_id, go_name in info['go_terms']:
            if go_id in target_go_terms:
                matching_gos.append(f"{go_id}:{go_name}")
        if matching_gos:
            matches[pfam_id] = {
                'name': info['name'],
                'go_terms': matching_gos
            }
    return matches


def writeTsv(filepath, pfam_dict, source_label):
    """write pfam domain list to TSV file."""
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['pfam_id', 'pfam_name', 'go_terms', 'source'])
        for pfam_id in sorted(pfam_dict.keys()):
            info = pfam_dict[pfam_id]
            go_str = '; '.join(info['go_terms'])
            writer.writerow([pfam_id, info['name'], go_str, source_label])
    print(f"  wrote {len(pfam_dict)} entries to {filepath}")


def main():
    os.makedirs(RAW_DIR, exist_ok=True)

    # step 1: download pfam2go
    print("step 1: downloading pfam2go mapping...")
    pfam2go_path = os.path.join(RAW_DIR, 'pfam2go.txt')
    if not download_file(PFAM2GO_URL, pfam2go_path, "pfam2go"):
        print("FATAL: cannot download pfam2go")
        return 1

    # step 2: download EuRBPDB HMM archive (for cross-validation)
    print("\nstep 2: downloading EuRBPDB 791 RBD archive...")
    eurbpdb_path = os.path.join(RAW_DIR, 'eurbpdb_791_rbds.zip')
    eurbpdb_ok = download_file(EURBPDB_URL, eurbpdb_path, "EuRBPDB RBDs")

    # step 3: get GO term descendants from quickgo
    print("\nstep 3: fetching GO term hierarchies from QuickGO...")
    rna_go_terms = fetchGoDescendants(RNA_BINDING_GO)
    time.sleep(1)  # be polite to API
    dna_go_terms = fetchGoDescendants(DNA_BINDING_GO)
    time.sleep(1)
    na_go_terms = fetchGoDescendants(NA_BINDING_GO)

    # step 4: parse pfam2go and find matches
    print("\nstep 4: parsing pfam2go and finding matches...")
    pfam_go = parsePfam2go(pfam2go_path)
    print(f"  total pfam entries with GO annotations: {len(pfam_go)}")

    rna_pfams = findPfamsByGoTerms(pfam_go, rna_go_terms)
    dna_pfams = findPfamsByGoTerms(pfam_go, dna_go_terms)
    na_pfams = findPfamsByGoTerms(pfam_go, na_go_terms)

    # combine: any pfam that binds RNA, DNA, or general NA
    all_na_pfams = {}
    for pfam_id in set(list(rna_pfams.keys()) + list(dna_pfams.keys()) + list(na_pfams.keys())):
        binding_types = []
        if pfam_id in rna_pfams:
            binding_types.append('RNA')
        if pfam_id in dna_pfams:
            binding_types.append('DNA')
        if pfam_id in na_pfams and pfam_id not in rna_pfams and pfam_id not in dna_pfams:
            binding_types.append('NA_general')

        # get name and go terms from whichever source has them
        name = (rna_pfams.get(pfam_id, {}).get('name') or
                dna_pfams.get(pfam_id, {}).get('name') or
                na_pfams.get(pfam_id, {}).get('name', ''))
        go_terms = []
        for d in [rna_pfams, dna_pfams, na_pfams]:
            if pfam_id in d:
                go_terms.extend(d[pfam_id]['go_terms'])
        go_terms = list(dict.fromkeys(go_terms))  # deduplicate preserving order

        all_na_pfams[pfam_id] = {
            'name': name,
            'go_terms': go_terms,
            'binding_type': '+'.join(binding_types)
        }

    # step 5: cross-validate with EuRBPDB
    eurbpdb_domains = set()
    if eurbpdb_ok:
        print("\nstep 5: cross-validating with EuRBPDB domains...")
        eurbpdb_domains = extractEurbpdbDomains(eurbpdb_path)
        print(f"  EuRBPDB domain names extracted: {len(eurbpdb_domains)}")

        # try to match EuRBPDB domain names to pfam names
        pfam_name_to_id = {}
        for pfam_id, info in pfam_go.items():
            pfam_name_to_id[info['name'].lower()] = pfam_id

        eurbpdb_matched = 0
        eurbpdb_new = 0
        for domain in eurbpdb_domains:
            # eurbpdb uses pfam short names
            if domain.lower() in pfam_name_to_id:
                pid = pfam_name_to_id[domain.lower()]
                if pid not in rna_pfams:
                    eurbpdb_new += 1
                    rna_pfams[pid] = {
                        'name': domain,
                        'go_terms': ['EuRBPDB:RNA-binding domain (EuRBPDB 791)']
                    }
                    if pid not in all_na_pfams:
                        all_na_pfams[pid] = {
                            'name': domain,
                            'go_terms': ['EuRBPDB:RNA-binding domain'],
                            'binding_type': 'RNA'
                        }
                else:
                    eurbpdb_matched += 1

        print(f"  EuRBPDB domains already in GO-based RNA set: {eurbpdb_matched}")
        print(f"  EuRBPDB domains added (new): {eurbpdb_new}")

    # step 6: write output TSVs
    print("\nstep 6: writing output files...")
    writeTsv(os.path.join(OUT_DIR, 'rna_binding_pfams.tsv'), rna_pfams, 'GO+EuRBPDB')
    writeTsv(os.path.join(OUT_DIR, 'dna_binding_pfams.tsv'), dna_pfams, 'GO')

    # write combined NA-binding TSV with binding_type column
    combined_path = os.path.join(OUT_DIR, 'na_binding_pfams.tsv')
    with open(combined_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['pfam_id', 'pfam_name', 'binding_type', 'go_terms', 'source'])
        for pfam_id in sorted(all_na_pfams.keys()):
            info = all_na_pfams[pfam_id]
            go_str = '; '.join(info['go_terms'])
            writer.writerow([pfam_id, info['name'], info['binding_type'], go_str, 'GO+EuRBPDB'])
    print(f"  wrote {len(all_na_pfams)} entries to {combined_path}")

    # step 7: summary and verification
    print("\n--- SUMMARY ---")
    print(f"  RNA-binding pfam domains:     {len(rna_pfams)}")
    print(f"  DNA-binding pfam domains:     {len(dna_pfams)}")
    print(f"  General NA-binding pfam domains: {len(na_pfams)}")
    print(f"  Total unique NA-binding:      {len(all_na_pfams)}")

    # identify dual binders
    dual_binders = set(rna_pfams.keys()) & set(dna_pfams.keys())
    rna_only = set(rna_pfams.keys()) - set(dna_pfams.keys())
    dna_only = set(dna_pfams.keys()) - set(rna_pfams.keys())
    print(f"\n  RNA-only:   {len(rna_only)}")
    print(f"  DNA-only:   {len(dna_only)}")
    print(f"  dual (RNA+DNA): {len(dual_binders)}")

    # spot checks
    print("\n--- SPOT CHECKS ---")
    spot_checks = {
        'PF00076': ('RRM_1', 'RNA'),
        'PF00046': ('Homeobox', 'DNA'),
        'PF00013': ('KH_1', 'RNA'),
        'PF00271': ('Helicase_C', 'RNA or DNA'),
        'PF00009': ('GTP_EFTU', 'neither (translation factor)'),
    }
    for pfam_id, (expected_name, expected_binding) in spot_checks.items():
        in_rna = pfam_id in rna_pfams
        in_dna = pfam_id in dna_pfams
        actual = 'RNA+DNA' if (in_rna and in_dna) else 'RNA' if in_rna else 'DNA' if in_dna else 'none'
        status = "OK" if (
            (expected_binding == 'RNA' and in_rna) or
            (expected_binding == 'DNA' and in_dna) or
            (expected_binding == 'neither' and not in_rna and not in_dna)
        ) else "CHECK"
        print(f"  {pfam_id} ({expected_name}): expected={expected_binding}, actual={actual} [{status}]")

    # checkpoint
    print("\n--- CHECKPOINT ---")
    if len(rna_pfams) < 100:
        print(f"  WARNING: RNA-binding count ({len(rna_pfams)}) seems low (<100)")
        print(f"  consider supplementing with literature curation")
    else:
        print(f"  OK: RNA-binding count ({len(rna_pfams)}) is reasonable")

    if len(dna_pfams) < 100:
        print(f"  WARNING: DNA-binding count ({len(dna_pfams)}) seems low (<100)")
    else:
        print(f"  OK: DNA-binding count ({len(dna_pfams)}) is reasonable")

    return 0


if __name__ == '__main__':
    sys.exit(main())
