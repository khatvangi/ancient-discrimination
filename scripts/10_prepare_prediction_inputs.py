#!/usr/bin/env python3
"""prepare AF3/Boltz-2 input files for all predictions.

for each ortholog:
  - download protein sequence from UniProt
  - create job: protein + poly-U RNA (10-mer)
  - create job: protein + poly-dT DNA (10-mer)

output:
  data/sequences/          — individual FASTA files per ortholog
  data/substrates.fasta    — RNA and DNA substrate sequences
  results/phase2_prediction_manifest.tsv — job manifest
"""

import os
import csv
import sys
import time
import requests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')


def fetchUniProtSequence(accession):
    """fetch protein sequence from UniProt in FASTA format."""
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        lines = resp.text.strip().split('\n')
        header = lines[0]
        seq = ''.join(lines[1:])
        return header, seq
    except Exception as e:
        print(f"  ERROR fetching {accession}: {e}")
        return None, None


def main():
    # load orthologs
    orth_path = os.path.join(PROJECT_DIR, 'results', 'phase2_orthologs.tsv')
    orthologs = []
    with open(orth_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            orthologs.append(row)

    print(f"preparing prediction inputs for {len(orthologs)} orthologs")

    # create output directories
    seq_dir = os.path.join(PROJECT_DIR, 'data', 'sequences')
    os.makedirs(seq_dir, exist_ok=True)

    # substrate sequences
    rna_10mer = 'UUUUUUUUUU'  # poly-U RNA
    dna_10mer = 'TTTTTTTTTT'  # poly-dT DNA (Boltz-2 uses DNA notation)

    # download protein sequences and build manifest
    manifest = []
    downloaded = 0
    failed = 0

    for i, orth in enumerate(orthologs):
        acc = orth['uniprot_acc']
        pfam_id = orth['pfam_id']
        species = orth['species_short']

        # check if already downloaded
        fasta_path = os.path.join(seq_dir, f'{acc}.fasta')
        if os.path.exists(fasta_path):
            # read existing
            with open(fasta_path) as f:
                header = f.readline().strip()
                seq = f.read().replace('\n', '')
        else:
            header, seq = fetchUniProtSequence(acc)
            time.sleep(0.2)  # rate limit

            if seq is None:
                failed += 1
                continue

            # write individual FASTA
            with open(fasta_path, 'w') as f:
                f.write(f'{header}\n{seq}\n')
            downloaded += 1

        if (i + 1) % 20 == 0:
            print(f"  processed {i+1}/{len(orthologs)} "
                  f"(downloaded {downloaded}, failed {failed})")

        # create job entries — one for RNA, one for DNA
        job_base = f'{pfam_id}_{species}_{acc}'

        manifest.append({
            'job_id': f'{job_base}_RNA',
            'pfam_id': pfam_id,
            'pfam_name': orth['pfam_name'],
            'binding_type': orth['binding_type'],
            'functional_category': orth['functional_category'],
            'tier': orth['tier'],
            'species_short': species,
            'species_name': orth['species_name'],
            'domain_of_life': orth['domain_of_life'],
            'uniprot_acc': acc,
            'protein_name': orth['protein_name'],
            'gene_name': orth['gene_name'],
            'protein_length': orth['length'],
            'protein_fasta': fasta_path,
            'substrate_type': 'RNA',
            'substrate_seq': rna_10mer,
            'substrate_notation': 'poly-U 10mer',
        })

        manifest.append({
            'job_id': f'{job_base}_DNA',
            'pfam_id': pfam_id,
            'pfam_name': orth['pfam_name'],
            'binding_type': orth['binding_type'],
            'functional_category': orth['functional_category'],
            'tier': orth['tier'],
            'species_short': species,
            'species_name': orth['species_name'],
            'domain_of_life': orth['domain_of_life'],
            'uniprot_acc': acc,
            'protein_name': orth['protein_name'],
            'gene_name': orth['gene_name'],
            'protein_length': orth['length'],
            'protein_fasta': fasta_path,
            'substrate_type': 'DNA',
            'substrate_seq': dna_10mer,
            'substrate_notation': 'poly-dT 10mer',
        })

    print(f"\n  total sequences downloaded: {downloaded}")
    print(f"  sequences already cached: {len(orthologs) - downloaded - failed}")
    print(f"  failed downloads: {failed}")

    # write substrate FASTA
    sub_path = os.path.join(PROJECT_DIR, 'data', 'substrates.fasta')
    with open(sub_path, 'w') as f:
        f.write(f'>poly_U_10mer RNA substrate\n{rna_10mer}\n')
        f.write(f'>poly_dT_10mer DNA substrate\n{dna_10mer}\n')
    print(f"  wrote substrates: {sub_path}")

    # write manifest
    manifest_path = os.path.join(PROJECT_DIR, 'results',
                                 'phase2_prediction_manifest.tsv')
    with open(manifest_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=[
            'job_id', 'pfam_id', 'pfam_name', 'binding_type',
            'functional_category', 'tier', 'species_short', 'species_name',
            'domain_of_life', 'uniprot_acc', 'protein_name', 'gene_name',
            'protein_length', 'protein_fasta', 'substrate_type',
            'substrate_seq', 'substrate_notation',
        ])
        writer.writeheader()
        for row in manifest:
            writer.writerow(row)

    print(f"  wrote manifest: {manifest_path} ({len(manifest)} jobs)")

    # summary
    print(f"\n{'=' * 70}")
    print(f"PREDICTION MANIFEST SUMMARY")
    print(f"{'=' * 70}")
    print(f"  total jobs: {len(manifest)}")
    print(f"    RNA predictions: {sum(1 for m in manifest if m['substrate_type'] == 'RNA')}")
    print(f"    DNA predictions: {sum(1 for m in manifest if m['substrate_type'] == 'DNA')}")
    print(f"  unique proteins: {len(set(m['uniprot_acc'] for m in manifest))}")
    print(f"  families: {len(set(m['pfam_id'] for m in manifest))}")

    # size distribution
    lengths = [int(m['protein_length']) for m in manifest if m['substrate_type'] == 'RNA']
    if lengths:
        print(f"\n  protein lengths:")
        print(f"    min: {min(lengths)} aa")
        print(f"    max: {max(lengths)} aa")
        print(f"    median: {sorted(lengths)[len(lengths)//2]} aa")
        large = sum(1 for l in lengths if l > 1000)
        print(f"    >1000 aa: {large} (may need domain extraction for Boltz-2)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
