#!/usr/bin/env python3
"""generate AF3 Server JSON inputs for large proteins (>1000 aa).

these proteins may OOM on the local 24GB Titan RTX GPUs,
so we submit them to alphafoldserver.com instead.

AF3 Server format uses:
  - dialect: "alphafoldserver"
  - version: 1
  - proteinChain, rnaSequence, dnaSequence (not protein/rna/dna)
  - no 'id' field — uses 'count' instead
  - top level is a list of jobs

substrates:
  - RNA: poly-U 10mer (UUUUUUUUUU)
  - DNA: poly-dT 10mer (TTTTTTTTTT)

usage:
  python scripts/11b_generate_af3server_inputs.py [--min-length 1000]
"""

import argparse
import json
import os

import pandas as pd

PROJECT = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
MANIFEST = os.path.join(PROJECT, "results", "phase2_prediction_manifest.tsv")
SEQ_DIR = os.path.join(PROJECT, "data", "sequences")
OUTPUT_DIR = os.path.join(PROJECT, "structures", "af3server_inputs")

RNA_SUBSTRATE = "UUUUUUUUUU"
DNA_SUBSTRATE = "TTTTTTTTTT"


def readProteinSeq(fasta_path):
    """read protein sequence from single-entry FASTA."""
    seq_lines = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq_lines.append(line.strip())
    return "".join(seq_lines)


def makeServerJob(job_name, protein_seq, substrate_type, seed=1):
    """create a single AF3 Server job dict."""
    job = {
        "name": job_name,
        "modelSeeds": [seed],
        "sequences": [
            {"proteinChain": {"sequence": protein_seq, "count": 1}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }

    if substrate_type == "RNA":
        job["sequences"].append(
            {"rnaSequence": {"sequence": RNA_SUBSTRATE, "count": 1}}
        )
    elif substrate_type == "DNA":
        job["sequences"].append(
            {"dnaSequence": {"sequence": DNA_SUBSTRATE, "count": 1}}
        )

    return job


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-length", type=int, default=1000,
                        help="minimum protein length to include (default: 1000)")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(MANIFEST, sep="\t")

    # filter to large proteins
    large = df[df["protein_length"] > args.min_length].copy()
    n_proteins = large.drop_duplicates(subset=["uniprot_acc"]).shape[0]
    print(f"generating AF3 Server inputs for {n_proteins} proteins > {args.min_length} aa")
    print(f"total jobs: {len(large)} ({len(large)//2} RNA + {len(large)//2} DNA)")

    # generate individual JSON files (one job per file for easy upload)
    jobs_written = 0
    job_list = []  # for the combined batch file

    for _, row in large.iterrows():
        pfam = row["pfam_id"]
        sp = row["species_short"]
        acc = row["uniprot_acc"]
        substrate = row["substrate_type"]
        job_name = f"{pfam}_{sp}_{acc}_{substrate}"

        # read protein sequence
        fasta_path = os.path.join(SEQ_DIR, f"{acc}.fasta")
        if not os.path.exists(fasta_path):
            print(f"  WARNING: missing {fasta_path}, skipping {job_name}")
            continue

        protein_seq = readProteinSeq(fasta_path)
        tokens = len(protein_seq) + 10
        if tokens > 5000:
            print(f"  WARNING: {job_name} has {tokens} tokens (exceeds 5000 limit), skipping")
            continue

        job = makeServerJob(job_name, protein_seq, substrate)

        # write individual file (single job wrapped in list)
        out_path = os.path.join(OUTPUT_DIR, f"{job_name}.json")
        with open(out_path, "w") as f:
            json.dump([job], f, indent=2)

        job_list.append(job)
        jobs_written += 1

    # write combined batch file (all jobs in one list)
    batch_path = os.path.join(OUTPUT_DIR, "all_large_proteins.json")
    with open(batch_path, "w") as f:
        json.dump(job_list, f, indent=2)

    # write job manifest for tracking
    manifest_path = os.path.join(OUTPUT_DIR, "af3server_jobs.txt")
    with open(manifest_path, "w") as f:
        for _, row in large.iterrows():
            job_name = f"{row['pfam_id']}_{row['species_short']}_{row['uniprot_acc']}_{row['substrate_type']}"
            f.write(f"{job_name}\t{row['protein_length']}\t{row['protein_length']+10} tokens\n")

    print(f"\nwrote {jobs_written} individual JSON files to {OUTPUT_DIR}/")
    print(f"wrote combined batch: {batch_path}")
    print(f"wrote job manifest: {manifest_path}")
    print(f"\nto submit: upload individual JSONs at alphafoldserver.com")
    print(f"  or use the combined file for programmatic submission")
    print(f"  daily limit: 30 jobs/day → {jobs_written} jobs = {(jobs_written+29)//30} days")


if __name__ == "__main__":
    main()
