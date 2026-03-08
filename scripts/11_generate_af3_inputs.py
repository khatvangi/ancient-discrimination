#!/usr/bin/env python3
# NOTE: Phase 3 script (planned, not yet completed).
# generates input files for AF3 structural predictions.
# the full prediction campaign has not been run.
"""generate AF3 JSON input files for all 296 prediction jobs.

AF3 input format uses:
  - proteinChain: {"sequence": "...", "count": 1}
  - rnaSequence: {"sequence": "UUUUUUUUUU", "count": 1}
  - dnaSequence: {"sequence": "TTTTTTTTTT", "count": 1}

for each ortholog, creates two JSON files:
  - {pfam}_{species}_{acc}_RNA.json → protein + poly-U 10mer
  - {pfam}_{species}_{acc}_DNA.json → protein + poly-dT 10mer

uses 1 seed per job for speed (296 jobs × ~30 min = ~148 GPU-hours).
with 2 GPUs running in parallel: ~3 days wall time.

outputs:
  structures/af3_inputs/  — 296 JSON files
  structures/af3_jobs.txt — list of all job JSON paths (for SLURM array)
"""

import os
import csv
import json
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')


def readFasta(fasta_path):
    """read a single-sequence FASTA file, return sequence string."""
    with open(fasta_path) as f:
        lines = f.readlines()
    seq = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return seq


def makeAf3Json(job_name, protein_seq, substrate_type, substrate_seq, seed=1):
    """create AF3 input JSON dict.

    uses the native AlphaFold 3 JSON format (dialect='alphafold3', version=1)
    which expects: 'protein', 'rna', 'dna' keys with 'id' fields.
    """
    job = {
        'dialect': 'alphafold3',
        'version': 1,
        'name': job_name,
        'modelSeeds': [seed],
        'sequences': [
            {
                'protein': {
                    'id': 'A',
                    'sequence': protein_seq,
                }
            },
        ],
    }

    if substrate_type == 'RNA':
        job['sequences'].append({
            'rna': {
                'id': 'B',
                'sequence': substrate_seq,
            }
        })
    elif substrate_type == 'DNA':
        job['sequences'].append({
            'dna': {
                'id': 'B',
                'sequence': substrate_seq,
            }
        })

    return job


def main():
    # load manifest
    manifest_path = os.path.join(PROJECT_DIR, 'results',
                                 'phase2_prediction_manifest.tsv')
    jobs = []
    with open(manifest_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            jobs.append(row)

    print(f"generating AF3 inputs for {len(jobs)} jobs")

    # create output directory
    input_dir = os.path.join(PROJECT_DIR, 'structures', 'af3_inputs')
    os.makedirs(input_dir, exist_ok=True)

    job_list = []  # paths for SLURM array
    skipped = 0

    for i, job in enumerate(jobs):
        job_id = job['job_id']
        fasta_path = job['protein_fasta']

        # read protein sequence
        if not os.path.exists(fasta_path):
            print(f"  SKIP {job_id}: FASTA not found at {fasta_path}")
            skipped += 1
            continue

        protein_seq = readFasta(fasta_path)
        if not protein_seq:
            print(f"  SKIP {job_id}: empty sequence")
            skipped += 1
            continue

        # create AF3 JSON
        af3_input = makeAf3Json(
            job_name=job_id,
            protein_seq=protein_seq,
            substrate_type=job['substrate_type'],
            substrate_seq=job['substrate_seq'],
        )

        # write JSON
        json_path = os.path.join(input_dir, f'{job_id}.json')
        with open(json_path, 'w') as f:
            json.dump(af3_input, f, indent=2)

        job_list.append(json_path)

    # write job list for SLURM array
    joblist_path = os.path.join(PROJECT_DIR, 'structures', 'af3_jobs.txt')
    with open(joblist_path, 'w') as f:
        for p in job_list:
            f.write(p + '\n')

    print(f"\n  generated: {len(job_list)} JSON files in {input_dir}")
    print(f"  skipped:   {skipped}")
    print(f"  job list:  {joblist_path}")

    # estimate compute
    sizes = []
    for job in jobs:
        try:
            sizes.append(int(job['protein_length']))
        except (ValueError, KeyError):
            pass

    if sizes:
        median = sorted(sizes)[len(sizes) // 2]
        total_tokens = sum(s + 10 for s in sizes)  # protein + 10mer
        print(f"\n  token stats:")
        print(f"    median protein: {median} aa")
        print(f"    total tokens across all jobs: {total_tokens:,}")
        print(f"    estimated time per job: ~15-45 min (GPU-dependent)")
        print(f"    estimated total GPU-hours: ~{len(job_list) * 0.5:.0f}")
        print(f"    with 2 GPUs: ~{len(job_list) * 0.5 / 2:.0f} wall-hours "
              f"(~{len(job_list) * 0.5 / 2 / 24:.1f} days)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
