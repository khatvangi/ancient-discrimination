#!/usr/bin/env python3
"""fetch representative sequences for 87 LUCA NA-binding domains and run ProNA2020.

for each Pfam domain, queries UniProt REST API for a reviewed (Swiss-Prot)
E. coli K12 sequence. falls back to model archaea and bacteria if not found.
writes combined FASTA, runs ProNA2020, parses output into TSV.

usage:
  # step 1: fetch sequences
  python 20_run_prona2020.py fetch

  # step 2: run ProNA2020 (must use prona2020 conda env)
  # source activate prona2020 && python /storage/kiran-stuff/ProNA2020/run_prona_protlevel.py data/luca_na_representatives.fasta

  # step 3: parse results
  python 20_run_prona2020.py parse <prona_output_file>
"""

import os
import sys
import csv
import time
import json
import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')
CENSUS_FILE = os.path.join(PROJECT_DIR, 'results', 'phase1_census.tsv')
OUTPUT_FASTA = os.path.join(PROJECT_DIR, 'data', 'luca_na_representatives.fasta')
OUTPUT_TSV = os.path.join(PROJECT_DIR, 'results', 'prona2020_luca_domains.tsv')

# model organisms to search, in priority order
# taxon IDs: E. coli K12 = 83333, M. jannaschii = 243232,
# T. thermophilus HB27 = 262724, B. subtilis 168 = 224308,
# S. solfataricus P2 = 273057, H. volcanii = 309800
MODEL_ORGANISMS = [
    (83333, "Escherichia coli K12"),
    (243232, "Methanocaldococcus jannaschii"),
    (262724, "Thermus thermophilus HB27"),
    (224308, "Bacillus subtilis 168"),
    (273057, "Sulfolobus solfataricus P2"),
    (309800, "Haloferax volcanii"),
]

# broader fallback: any reviewed sequence from these taxa (bigger nets)
BROAD_TAXA = [
    (83333, "E. coli K12"),
    (2157, "Archaea"),          # all archaea
    (2, "Bacteria"),            # all bacteria
]


def loadLucaNaDomains():
    """load the 87 LUCA/preLUCA domains with RNA or DNA binding."""
    domains = []
    with open(CENSUS_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            age = row['age_class']
            if age in ('LUCA', 'preLUCA'):
                rna = row['rna_union'] == 'True'
                dna = row['dna_binding'] == 'True'
                if rna or dna:
                    domains.append({
                        'pfam_id': row['pfam_id'],
                        'age_class': age,
                        'rna_union': rna,
                        'dna_binding': dna,
                    })
    return domains


def queryUniprot(pfam_id, taxon_id, reviewed_only=True):
    """query UniProt REST API for proteins with a given Pfam domain in a given organism.
    returns list of (accession, entry_name, organism, sequence) tuples.
    """
    # use the new UniProt REST API
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    # build query
    query_parts = [f'(xref:pfam-{pfam_id})']
    if taxon_id:
        query_parts.append(f'(taxonomy_id:{taxon_id})')
    if reviewed_only:
        query_parts.append('(reviewed:true)')

    query = ' AND '.join(query_parts)

    params = {
        'query': query,
        'format': 'json',
        'size': 3,  # just need top few
        'fields': 'accession,id,organism_name,sequence,protein_name',
    }

    try:
        resp = requests.get(base_url, params=params, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            results = []
            for entry in data.get('results', []):
                acc = entry.get('primaryAccession', '')
                entry_name = entry.get('uniProtkbId', '')
                org = entry.get('organism', {}).get('scientificName', '')
                seq = entry.get('sequence', {}).get('value', '')
                prot_name = ''
                if 'proteinDescription' in entry:
                    pd = entry['proteinDescription']
                    if 'recommendedName' in pd:
                        prot_name = pd['recommendedName'].get('fullName', {}).get('value', '')
                    elif 'submissionNames' in pd and pd['submissionNames']:
                        prot_name = pd['submissionNames'][0].get('fullName', {}).get('value', '')
                results.append((acc, entry_name, org, seq, prot_name))
            return results
        else:
            return []
    except Exception as e:
        print(f"  warning: API error for {pfam_id}/{taxon_id}: {e}")
        return []


def fetchRepresentativeSequences(domains):
    """for each domain, find a representative protein sequence."""
    representatives = {}

    for i, dom in enumerate(domains):
        pfam_id = dom['pfam_id']
        print(f"[{i+1}/{len(domains)}] {pfam_id}...", end=' ', flush=True)

        found = False

        # try model organisms first (reviewed only)
        for taxon_id, org_name in MODEL_ORGANISMS:
            results = queryUniprot(pfam_id, taxon_id, reviewed_only=True)
            if results:
                acc, entry_name, org, seq, prot_name = results[0]
                # skip very short sequences (domain fragments)
                if len(seq) < 30:
                    continue
                representatives[pfam_id] = {
                    'accession': acc,
                    'entry_name': entry_name,
                    'organism': org,
                    'sequence': seq,
                    'protein_name': prot_name,
                    'seq_length': len(seq),
                }
                print(f"found {acc} ({org}, {len(seq)} aa)")
                found = True
                break
            time.sleep(0.3)  # be polite to the API

        if found:
            continue

        # fallback: any reviewed sequence
        results = queryUniprot(pfam_id, None, reviewed_only=True)
        if results:
            # prefer the shortest full-length hit (less likely to be a huge multi-domain protein)
            # but still > 50 aa
            valid = [(acc, en, org, seq, pn) for acc, en, org, seq, pn in results if len(seq) >= 50]
            if valid:
                # sort by length
                valid.sort(key=lambda x: len(x[3]))
                acc, entry_name, org, seq, prot_name = valid[0]
                representatives[pfam_id] = {
                    'accession': acc,
                    'entry_name': entry_name,
                    'organism': org,
                    'sequence': seq,
                    'protein_name': prot_name,
                    'seq_length': len(seq),
                }
                print(f"found {acc} ({org}, {len(seq)} aa) [any reviewed]")
                found = True

        if not found:
            # last resort: unreviewed
            results = queryUniprot(pfam_id, 83333, reviewed_only=False)
            if not results:
                results = queryUniprot(pfam_id, None, reviewed_only=False)
            if results:
                valid = [(acc, en, org, seq, pn) for acc, en, org, seq, pn in results if len(seq) >= 50]
                if valid:
                    valid.sort(key=lambda x: len(x[3]))
                    acc, entry_name, org, seq, prot_name = valid[0]
                    representatives[pfam_id] = {
                        'accession': acc,
                        'entry_name': entry_name,
                        'organism': org,
                        'sequence': seq,
                        'protein_name': prot_name,
                        'seq_length': len(seq),
                    }
                    print(f"found {acc} ({org}, {len(seq)} aa) [unreviewed]")
                    found = True

        if not found:
            print("NOT FOUND")

        time.sleep(0.3)

    return representatives


def writeFasta(representatives, domains, output_path):
    """write representative sequences to a multi-FASTA file."""
    with open(output_path, 'w') as f:
        for dom in domains:
            pfam_id = dom['pfam_id']
            if pfam_id in representatives:
                rep = representatives[pfam_id]
                # header format: >sp|ACC|ENTRY_NAME PFAM_ID
                header = f">sp|{rep['accession']}|{rep['entry_name']} {pfam_id} {rep['protein_name']}"
                f.write(header + '\n')
                # write sequence in 60-char lines
                seq = rep['sequence']
                for j in range(0, len(seq), 60):
                    f.write(seq[j:j+60] + '\n')

    print(f"\nwrote {len(representatives)} sequences to {output_path}")


def writeMetadata(representatives, domains, output_path):
    """write metadata TSV for tracking which sequences were used."""
    meta_path = output_path.replace('.fasta', '_metadata.tsv')
    with open(meta_path, 'w') as f:
        f.write('pfam_id\tage_class\trna_union\tdna_binding\trep_accession\tentry_name\torganism\tseq_length\tprotein_name\n')
        for dom in domains:
            pfam_id = dom['pfam_id']
            if pfam_id in representatives:
                rep = representatives[pfam_id]
                f.write(f"{pfam_id}\t{dom['age_class']}\t{dom['rna_union']}\t{dom['dna_binding']}\t"
                        f"{rep['accession']}\t{rep['entry_name']}\t{rep['organism']}\t{rep['seq_length']}\t{rep['protein_name']}\n")
            else:
                f.write(f"{pfam_id}\t{dom['age_class']}\t{dom['rna_union']}\t{dom['dna_binding']}\t"
                        f"NOT_FOUND\t\t\t\t\n")
    print(f"wrote metadata to {meta_path}")


def parsePronaOutput(output_file):
    """parse ProNA2020 output file into structured data.

    ProNA2020 output format (from control results):
    Protein                    Len  DNA_prob DNA_bind  RNA_prob RNA_bind  Pro_prob Pro_bind  Nuc_prob
    ---...
    sp|P0A9X9|CSPA_ECOLI        70    0.0373       no    0.1804       no    0.6373      YES    0.0885
    """
    results = {}
    with open(output_file) as f:
        lines = f.readlines()

    for line in lines:
        line = line.strip()
        if not line or line.startswith('---') or line.startswith('Protein'):
            continue

        # parse the fixed-width format
        parts = line.split()
        if len(parts) < 9:
            continue

        protein_id = parts[0]
        try:
            seq_len = int(parts[1])
            dna_prob = float(parts[2])
            dna_bind = parts[3]
            rna_prob = float(parts[4])
            rna_bind = parts[5]
            pro_prob = float(parts[6])
            pro_bind = parts[7]
            nuc_prob = float(parts[8])
        except (ValueError, IndexError):
            continue

        # extract accession from protein_id (e.g. "sp|P0A9X9|CSPA_ECOLI" -> "P0A9X9")
        parts_id = protein_id.split('|')
        if len(parts_id) >= 2:
            accession = parts_id[1]
        else:
            accession = protein_id

        results[accession] = {
            'protein_id': protein_id,
            'seq_length': seq_len,
            'P_DNA': dna_prob,
            'DNA_bind': dna_bind == 'YES',
            'P_RNA': rna_prob,
            'RNA_bind': rna_bind == 'YES',
            'P_protein': pro_prob,
            'P_nucleotide': nuc_prob,
        }

    return results


def classifyBinding(pdata):
    """classify predicted binding using ProNA2020's own calls.
    ProNA2020 applies a nucleotide gate: DNA_bind requires Nuc_prob > 0.5 AND DNA_prob > 0.5.
    same for RNA_bind. we use ProNA's own YES/no calls stored in DNA_bind and RNA_bind.
    """
    dna_pos = pdata['DNA_bind']
    rna_pos = pdata['RNA_bind']
    if dna_pos and rna_pos:
        return "dual"
    elif rna_pos:
        return "RNA-only"
    elif dna_pos:
        return "DNA-only"
    else:
        return "non-binder"


def buildResultsTSV(prona_results, metadata_file, output_tsv):
    """merge ProNA2020 results with metadata into final TSV.

    iterates over ALL 87 metadata rows (not just unique accessions)
    so that multi-domain proteins get their ProNA2020 scores assigned
    to each Pfam domain they represent.
    """
    # load all metadata rows (preserving order, allowing duplicate accessions)
    meta_rows = []
    with open(metadata_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['rep_accession'] != 'NOT_FOUND':
                meta_rows.append(row)

    with open(output_tsv, 'w') as f:
        f.write('pfam_id\tname\tage_class\trna_union\tdna_binding\trep_uniprot\torganism\tseq_length\t'
                'P_DNA\tP_RNA\tP_protein\tpredicted_class\n')

        for m in meta_rows:
            acc = m['rep_accession']
            if acc in prona_results:
                pdata = prona_results[acc]
                pred_class = classifyBinding(pdata)
                f.write(f"{m['pfam_id']}\t{m['protein_name']}\t{m['age_class']}\t{m['rna_union']}\t{m['dna_binding']}\t"
                        f"{acc}\t{m['organism']}\t{pdata['seq_length']}\t"
                        f"{pdata['P_DNA']:.4f}\t{pdata['P_RNA']:.4f}\t{pdata['P_protein']:.4f}\t{pred_class}\n")

    print(f"wrote results to {output_tsv}")

    # summary statistics — count per domain (all 87), not per unique protein
    classes = {}
    total = 0
    for m in meta_rows:
        acc = m['rep_accession']
        if acc in prona_results:
            pdata = prona_results[acc]
            pred_class = classifyBinding(pdata)
            classes[pred_class] = classes.get(pred_class, 0) + 1
            total += 1

    print(f"\n=== ProNA2020 RESULTS SUMMARY ===")
    print(f"total domains scored: {total}")
    for cls in ['dual', 'RNA-only', 'DNA-only', 'non-binder']:
        n = classes.get(cls, 0)
        pct = 100 * n / total if total > 0 else 0
        print(f"  {cls}: {n} ({pct:.1f}%)")

    # also report unique-protein level counts
    classes_unique = {}
    total_unique = 0
    seen_acc = set()
    for m in meta_rows:
        acc = m['rep_accession']
        if acc in prona_results and acc not in seen_acc:
            seen_acc.add(acc)
            pdata = prona_results[acc]
            pred_class = classifyBinding(pdata)
            classes_unique[pred_class] = classes_unique.get(pred_class, 0) + 1
            total_unique += 1

    print(f"\n=== UNIQUE PROTEINS (deduplicated) ===")
    print(f"total unique proteins scored: {total_unique}")
    for cls in ['dual', 'RNA-only', 'DNA-only', 'non-binder']:
        n = classes_unique.get(cls, 0)
        pct = 100 * n / total_unique if total_unique > 0 else 0
        print(f"  {cls}: {n} ({pct:.1f}%)")

    # the dual fraction assessment uses the per-domain count (87 domains)
    dual_n = classes.get('dual', 0)
    dual_pct = 100 * dual_n / total if total > 0 else 0
    na_binding_n = classes.get('dual', 0) + classes.get('RNA-only', 0) + classes.get('DNA-only', 0)
    na_binding_pct = 100 * na_binding_n / total if total > 0 else 0
    print(f"\n=== ASSESSMENT ===")
    print(f"dual-binders: {dual_n}/{total} = {dual_pct:.1f}%")
    print(f"any NA-binders (dual+RNA+DNA): {na_binding_n}/{total} = {na_binding_pct:.1f}%")
    if dual_pct > 30:
        print("→ STRONG support for ancestral generalism (>30% dual)")
    else:
        print(f"→ dual-binders below 30% threshold, but {na_binding_pct:.1f}% have some NA-binding signal")


def main():
    if len(sys.argv) < 2:
        print("usage: python 20_run_prona2020.py fetch|parse [prona_output_file]")
        sys.exit(1)

    mode = sys.argv[1]

    if mode == 'fetch':
        domains = loadLucaNaDomains()
        print(f"loaded {len(domains)} LUCA/preLUCA NA-binding domains")

        representatives = fetchRepresentativeSequences(domains)
        print(f"\nfound representatives for {len(representatives)}/{len(domains)} domains")

        writeFasta(representatives, domains, OUTPUT_FASTA)
        writeMetadata(representatives, domains, OUTPUT_FASTA)

    elif mode == 'parse':
        if len(sys.argv) < 3:
            print("usage: python 20_run_prona2020.py parse <prona_output_file>")
            sys.exit(1)

        prona_output = sys.argv[2]
        metadata_file = OUTPUT_FASTA.replace('.fasta', '_metadata.tsv')

        if not os.path.exists(metadata_file):
            print(f"error: metadata file not found: {metadata_file}")
            print("run 'fetch' step first")
            sys.exit(1)

        prona_results = parsePronaOutput(prona_output)
        print(f"parsed {len(prona_results)} ProNA2020 results")

        buildResultsTSV(prona_results, metadata_file, OUTPUT_TSV)

    else:
        print(f"unknown mode: {mode}")
        sys.exit(1)


if __name__ == '__main__':
    main()
