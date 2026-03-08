#!/usr/bin/env python3
"""for each selected family, identify one ortholog per species.

species panel (3 per domain of life):
  bacteria:  eco (E. coli K-12), bsu (B. subtilis), tth (T. thermophilus)
  archaea:   mja (M. jannaschii), hvo (H. volcanii), sso (S. solfataricus)
  eukarya:   sce (S. cerevisiae), hsa (H. sapiens), ath (A. thaliana)

method:
  - query UniProt REST API for reviewed (Swiss-Prot) entries with target Pfam + species
  - prefer entries with: 3D structure, experimental evidence
  - fall back to any reviewed entry if no structure exists

outputs:
  results/phase2_orthologs.tsv — one row per ortholog
"""

import os
import csv
import json
import sys
import time
import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')

# species panel: taxon_id, short_name, domain_of_life, species_name
SPECIES_PANEL = [
    # bacteria
    (83333,  'eco', 'Bacteria', 'Escherichia coli K-12'),
    (224308, 'bsu', 'Bacteria', 'Bacillus subtilis 168'),
    (274,    'tth', 'Bacteria', 'Thermus thermophilus'),
    # archaea
    (243232, 'mja', 'Archaea', 'Methanocaldococcus jannaschii DSM 2661'),
    (2246,   'hvo', 'Archaea', 'Haloferax volcanii'),
    (273057, 'sso', 'Archaea', 'Sulfolobus solfataricus P2'),
    # eukarya
    (559292, 'sce', 'Eukarya', 'Saccharomyces cerevisiae S288C'),
    (9606,   'hsa', 'Eukarya', 'Homo sapiens'),
    (3702,   'ath', 'Eukarya', 'Arabidopsis thaliana'),
]

# UniProt REST API base
UNIPROT_API = 'https://rest.uniprot.org/uniprotkb/search'


def queryUniProt(pfam_id, taxon_id, max_results=5):
    """query UniProt for reviewed entries with a given Pfam domain in a given species.

    returns list of dicts with: accession, entry_name, protein_name, length,
    has_structure (bool), gene_name
    """
    # use UniProt's query syntax
    query = f'(xref:pfam-{pfam_id}) AND (taxonomy_id:{taxon_id}) AND (reviewed:true)'

    params = {
        'query': query,
        'format': 'json',
        'size': max_results,
        'fields': 'accession,id,protein_name,length,gene_names,xref_pdb,organism_name',
    }

    try:
        resp = requests.get(UNIPROT_API, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as e:
        # try alternate query syntax if first fails
        query = f'({pfam_id}) AND (taxonomy_id:{taxon_id}) AND (reviewed:true)'
        params['query'] = query
        try:
            resp = requests.get(UNIPROT_API, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e2:
            return []

    results = []
    for entry in data.get('results', []):
        acc = entry.get('primaryAccession', '')
        entry_name = entry.get('uniProtkbId', '')

        # extract protein name
        pname = ''
        if 'proteinDescription' in entry:
            rec = entry['proteinDescription'].get('recommendedName', {})
            if rec:
                pname = rec.get('fullName', {}).get('value', '')
            if not pname:
                sub = entry['proteinDescription'].get('submissionNames', [])
                if sub:
                    pname = sub[0].get('fullName', {}).get('value', '')

        length = entry.get('sequence', {}).get('length', 0)

        # check for PDB cross-references
        pdb_ids = []
        for xref in entry.get('uniProtKBCrossReferences', []):
            if xref.get('database') == 'PDB':
                pdb_ids.append(xref.get('id', ''))

        # extract gene name
        gene = ''
        genes = entry.get('genes', [])
        if genes:
            gene = genes[0].get('geneName', {}).get('value', '')
            if not gene:
                ordered = genes[0].get('orderedLocusNames', [])
                if ordered:
                    gene = ordered[0].get('value', '')

        organism = entry.get('organism', {}).get('scientificName', '')

        results.append({
            'accession': acc,
            'entry_name': entry_name,
            'protein_name': pname,
            'gene_name': gene,
            'length': length,
            'has_structure': len(pdb_ids) > 0,
            'pdb_count': len(pdb_ids),
            'pdb_ids': ','.join(pdb_ids[:5]),  # first 5
            'organism': organism,
        })

    return results


def main():
    # load selected families
    sel_path = os.path.join(PROJECT_DIR, 'results', 'phase2_selected_families.tsv')
    families = []
    with open(sel_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            families.append(row)

    print(f"querying orthologs for {len(families)} families × {len(SPECIES_PANEL)} species")

    all_orthologs = []
    missing = []

    for fi, fam in enumerate(families):
        pfam_id = fam['pfam_id']
        pfam_name = fam['pfam_name']
        print(f"\n[{fi+1}/{len(families)}] {pfam_id} ({pfam_name})")

        for taxon_id, short_name, domain, species_name in SPECIES_PANEL:
            hits = queryUniProt(pfam_id, taxon_id)
            time.sleep(0.3)  # rate limit

            if hits:
                # prefer entry with PDB structure
                best = None
                for h in hits:
                    if h['has_structure']:
                        best = h
                        break
                if best is None:
                    best = hits[0]  # take first reviewed entry

                all_orthologs.append({
                    'pfam_id': pfam_id,
                    'pfam_name': pfam_name,
                    'binding_type': fam['binding_type'],
                    'functional_category': fam['functional_category'],
                    'tier': fam['tier'],
                    'species_short': short_name,
                    'species_name': species_name,
                    'domain_of_life': domain,
                    'taxon_id': taxon_id,
                    'uniprot_acc': best['accession'],
                    'entry_name': best['entry_name'],
                    'protein_name': best['protein_name'],
                    'gene_name': best['gene_name'],
                    'length': best['length'],
                    'has_structure': best['has_structure'],
                    'pdb_count': best['pdb_count'],
                    'pdb_ids': best['pdb_ids'],
                })
                status = f"  {short_name}: {best['accession']} ({best['gene_name']}) " \
                         f"len={best['length']}"
                if best['has_structure']:
                    status += f" PDB={best['pdb_ids'][:30]}"
                print(status)
            else:
                missing.append((pfam_id, short_name, species_name))
                print(f"  {short_name}: NO REVIEWED ENTRY")

    # write output
    out_path = os.path.join(PROJECT_DIR, 'results', 'phase2_orthologs.tsv')
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=[
            'pfam_id', 'pfam_name', 'binding_type', 'functional_category',
            'tier', 'species_short', 'species_name', 'domain_of_life',
            'taxon_id', 'uniprot_acc', 'entry_name', 'protein_name',
            'gene_name', 'length', 'has_structure', 'pdb_count', 'pdb_ids',
        ])
        writer.writeheader()
        for row in all_orthologs:
            writer.writerow(row)

    print(f"\n{'=' * 70}")
    print(f"ORTHOLOG SELECTION SUMMARY")
    print(f"{'=' * 70}")
    print(f"  total orthologs found: {len(all_orthologs)}")
    print(f"  missing (no reviewed entry): {len(missing)}")
    print(f"  with PDB structure: {sum(1 for o in all_orthologs if o['has_structure'])}")
    print(f"  wrote: {out_path}")

    # coverage by family
    from collections import Counter, defaultdict
    fam_coverage = defaultdict(lambda: {'bac': 0, 'arc': 0, 'euk': 0, 'total': 0})
    for o in all_orthologs:
        pid = o['pfam_id']
        fam_coverage[pid]['total'] += 1
        if o['domain_of_life'] == 'Bacteria':
            fam_coverage[pid]['bac'] += 1
        elif o['domain_of_life'] == 'Archaea':
            fam_coverage[pid]['arc'] += 1
        elif o['domain_of_life'] == 'Eukarya':
            fam_coverage[pid]['euk'] += 1

    print(f"\n  {'pfam_id':10s} {'name':20s} {'Bac':>4s} {'Arc':>4s} {'Euk':>4s} {'Tot':>4s} {'pass':>5s}")
    print(f"  {'-'*10} {'-'*20} {'-'*4} {'-'*4} {'-'*4} {'-'*4} {'-'*5}")
    pass_count = 0
    for fam in families:
        pid = fam['pfam_id']
        c = fam_coverage[pid]
        domains_hit = sum(1 for x in [c['bac'], c['arc'], c['euk']] if x > 0)
        passes = domains_hit >= 2 and c['total'] >= 5
        if passes:
            pass_count += 1
        print(f"  {pid:10s} {fam['pfam_name']:20s} {c['bac']:4d} {c['arc']:4d} "
              f"{c['euk']:4d} {c['total']:4d} {'YES' if passes else 'no':>5s}")

    print(f"\n  families passing (≥5 orthologs, ≥2 domains): {pass_count}/{len(families)}")

    if missing:
        print(f"\n  missing entries:")
        for pid, sp, spname in missing:
            print(f"    {pid} × {sp} ({spname})")

    return 0


if __name__ == '__main__':
    sys.exit(main())
