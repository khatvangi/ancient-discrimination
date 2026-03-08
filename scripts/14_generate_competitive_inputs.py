#!/usr/bin/env python3
"""generate AF3 Server competitive binding inputs.

competitive design: one protein + RNA + dsDNA in a single prediction.
AF3 must allocate the binding interface to one ligand or the other,
making preference directly visible within a single run.

pilot: 9 proteins across 3 families (CSD, Alba, KH) with known preferences.
full batch: all 148 manifest proteins in competitive mode.

AF3 Server format (verified working):
  - dialect: "alphafoldserver", version: 1
  - proteinChain, rnaSequence, dnaSequence with count field
  - top level is a list of jobs

usage:
  python scripts/14_generate_competitive_inputs.py --mode pilot
  python scripts/14_generate_competitive_inputs.py --mode full
  python scripts/14_generate_competitive_inputs.py --mode both
"""

import argparse
import json
import os
import sys
import time

import requests

PROJECT = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"

# ============================================================================
# nucleic acid probes — controlled variable, same for all jobs
# 15-nt stem-loop: generic, no strong sequence bias
# RNA and DNA have identical base sequence for fair comparison
# ============================================================================
RNA_PROBE = "GACUGAUUCGAUCAG"       # 15-nt ssRNA
DNA_SENSE = "GACTGATTCGATCAG"       # 15-nt sense strand
DNA_ANTISENSE = "CTGATCGAATCAGTC"   # 15-nt antisense (reverse complement)

# ============================================================================
# pilot protein definitions — 9 proteins across 3 ancient families
# each family has: bacteria/archaea/eukarya representative(s)
# crop boundaries verified against InterPro Feb 2026
# ============================================================================
PILOT_PROTEINS = [
    # --- family 1: CSD (Cold Shock Domain, PF00313) ---
    # ancestral RNA chaperone → derived DNA transcription factor
    {
        "id": "CspA_Ecoli",
        "family": "CSD",
        "uniprot": "P0A9X9",
        "organism": "E. coli",
        "domain_of_life": "Bacteria",
        "description": "cold shock protein A — RNA chaperone",
        "expected": "RNA",
        "crop": None,  # 70 aa, single domain — use full sequence
    },
    {
        "id": "CspB_Hfx",
        "family": "CSD",
        "uniprot": "A0A0D6JQU8",
        "organism": "H. massiliensis",
        "domain_of_life": "Archaea",
        "description": "cold shock protein CspB (archaeal CSD, halophile)",
        "expected": "unknown",
        "crop": None,  # 74 aa, single CSD domain
    },
    {
        "id": "YBX1_Human",
        "family": "CSD",
        "uniprot": "P67809",
        "organism": "H. sapiens",
        "domain_of_life": "Eukarya",
        "description": "Y-box binding protein 1 — DNA transcription factor",
        "expected": "DNA",
        # InterPro: CSD (PF00313) at 61-127. crop wider to include flanking.
        "crop": (52, 135),
    },
    # --- family 2: Alba (PF01918) ---
    # archaeal chromatin protein, RNA-binding ancestral → DNA-binding derived
    {
        "id": "Sso10b_Ssolf",
        "family": "Alba",
        "uniprot": "P60849",
        "organism": "S. solfataricus",
        "domain_of_life": "Archaea",
        "description": "DNA/RNA-binding protein Alba 1 (Sso10b)",
        "expected": "DNA",
        "crop": None,  # 97 aa, single Alba domain
    },
    {
        "id": "Alba_Mjann",
        "family": "Alba",
        "uniprot": "Q58620",
        "organism": "M. jannaschii",
        "domain_of_life": "Archaea",
        "description": "Alba homolog — methanogen",
        "expected": "unknown",
        "crop": None,  # 92 aa, single domain
    },
    {
        "id": "Rpp25_Human",
        "family": "Alba",
        "uniprot": "Q9BUL9",
        "organism": "H. sapiens",
        "domain_of_life": "Eukarya",
        "description": "RNase P/MRP subunit — RNA binding",
        "expected": "RNA",
        "crop": None,  # 199 aa, Alba-like fold is the whole protein
    },
    # --- family 3: KH domain (PF00013 / IPR004087) ---
    # ancient RNA-binding domain, some derived DNA-binding members
    {
        "id": "NusA_KH_Ecoli",
        "family": "KH",
        "uniprot": "P0AFF6",
        "organism": "E. coli",
        "domain_of_life": "Bacteria",
        "description": "NusA KH1+KH2 domains — RNA binding",
        "expected": "RNA",
        # InterPro: KH1 (PF13184) 199-276, KH2 (PF26594) 280-343
        # crop includes both KH domains with flanking
        "crop": (195, 350),
    },
    {
        "id": "NusA_Afulg",
        "family": "KH",
        "uniprot": "O28388",
        "organism": "A. fulgidus",
        "domain_of_life": "Archaea",
        "description": "NusA — compact archaeal form",
        "expected": "RNA",
        "crop": None,  # 139 aa — compact archaeal NusA, use full
    },
    {
        "id": "PCBP1_KH1_Human",
        "family": "KH",
        "uniprot": "Q15365",
        "organism": "H. sapiens",
        "domain_of_life": "Eukarya",
        "description": "poly(rC)-binding protein 1 — dual RNA/ssDNA binder",
        "expected": "dual",
        # InterPro: KH1 at 15-76. crop wider for flanking.
        "crop": (10, 82),
    },
]


def fetchUniprotSeq(accession, max_retries=3):
    """fetch protein sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            lines = resp.text.strip().split("\n")
            seq = "".join(l.strip() for l in lines if not l.startswith(">"))
            return seq
        except requests.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(2)
            else:
                raise RuntimeError(f"failed to fetch {accession}: {e}")


def cropSequence(seq, crop):
    """crop sequence to domain boundaries (1-indexed, inclusive)."""
    start, end = crop
    return seq[start - 1 : end]


def makeCompetitiveJob(name, protein_seq, seeds=5):
    """build one AF3 Server competitive binding JSON.

    chains:
      protein (proteinChain)
      RNA single strand (rnaSequence)
      DNA sense strand (dnaSequence)
      DNA antisense strand (dnaSequence)

    AF3 must place both nucleic acids — the one at the canonical
    binding surface "wins" the competition.
    """
    job = {
        "name": name,
        "modelSeeds": list(range(1, seeds + 1)),
        "sequences": [
            {"proteinChain": {"sequence": protein_seq, "count": 1}},
            {"rnaSequence": {"sequence": RNA_PROBE, "count": 1}},
            {"dnaSequence": {"sequence": DNA_SENSE, "count": 1}},
            {"dnaSequence": {"sequence": DNA_ANTISENSE, "count": 1}},
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }
    return job


def generatePilotInputs(outdir):
    """generate competitive inputs for 9-protein pilot."""
    os.makedirs(outdir, exist_ok=True)

    print("=" * 60)
    print("COMPETITIVE BINDING PILOT — 9 proteins × 3 families")
    print("=" * 60)
    print(f"RNA probe:     {RNA_PROBE} ({len(RNA_PROBE)} nt)")
    print(f"DNA sense:     {DNA_SENSE} ({len(DNA_SENSE)} nt)")
    print(f"DNA antisense: {DNA_ANTISENSE} ({len(DNA_ANTISENSE)} nt)")
    print(f"output:        {outdir}/")
    print("=" * 60)

    jobs = []
    manifest_rows = []

    for p in PILOT_PROTEINS:
        name = p["id"]
        acc = p["uniprot"]
        print(f"\n  {name} ({acc})...", end=" ", flush=True)

        # fetch sequence
        full_seq = fetchUniprotSeq(acc)
        print(f"{len(full_seq)} aa full", end="", flush=True)

        # crop if needed
        if p["crop"]:
            domain_seq = cropSequence(full_seq, p["crop"])
            crop_str = f"{p['crop'][0]}-{p['crop'][1]}"
            print(f" → {len(domain_seq)} aa (crop {crop_str})", end="", flush=True)
        else:
            domain_seq = full_seq
            crop_str = "full"

        # token count: protein + RNA (15) + DNA sense (15) + DNA antisense (15) = protein + 45
        tokens = len(domain_seq) + 45
        if tokens > 5000:
            print(f" SKIP (tokens={tokens} > 5000)")
            continue

        # sanity check sequence length
        if len(domain_seq) < 30:
            print(f" WARNING: only {len(domain_seq)} aa")
        if len(domain_seq) > 500:
            print(f" WARNING: {len(domain_seq)} aa — large for competitive design")

        # generate competitive job
        job_name = f"competitive_{name}"
        job = makeCompetitiveJob(job_name, domain_seq)

        # write individual JSON (wrapped in list for AF3 Server upload)
        path = os.path.join(outdir, f"{job_name}.json")
        with open(path, "w") as f:
            json.dump([job], f, indent=2)

        jobs.append(job)
        manifest_rows.append({
            "job_name": job_name,
            "id": name,
            "family": p["family"],
            "uniprot": acc,
            "organism": p["organism"],
            "domain_of_life": p["domain_of_life"],
            "expected": p["expected"],
            "full_length": len(full_seq),
            "domain_length": len(domain_seq),
            "crop": crop_str,
            "tokens": tokens,
            "description": p["description"],
        })
        print(f" → {tokens} tokens ✓")

    # write batch file (all jobs in one list)
    batch_path = os.path.join(outdir, "all_pilot_competitive.json")
    with open(batch_path, "w") as f:
        json.dump(jobs, f, indent=2)

    # write manifest TSV
    manifest_path = os.path.join(outdir, "pilot_manifest.tsv")
    if manifest_rows:
        cols = list(manifest_rows[0].keys())
        with open(manifest_path, "w") as f:
            f.write("\t".join(cols) + "\n")
            for row in manifest_rows:
                f.write("\t".join(str(row[c]) for c in cols) + "\n")

    print("\n" + "=" * 60)
    print(f"PILOT SUMMARY")
    print(f"  jobs generated:    {len(jobs)}")
    print(f"  individual JSONs:  {outdir}/competitive_*.json")
    print(f"  batch file:        {batch_path}")
    print(f"  manifest:          {manifest_path}")
    print(f"  AF3 Server quota:  {len(jobs)} jobs (fits in 1 day, limit=30)")
    print("=" * 60)

    return jobs


def generateFullInputs(outdir):
    """generate competitive inputs for all 148 manifest proteins."""
    import pandas as pd

    manifest_path = os.path.join(PROJECT, "results", "phase2_prediction_manifest.tsv")
    seq_dir = os.path.join(PROJECT, "data", "sequences")
    os.makedirs(outdir, exist_ok=True)

    df = pd.read_csv(manifest_path, sep="\t")

    # deduplicate: one job per protein (not per substrate)
    proteins = df.drop_duplicates(subset=["uniprot_acc"]).copy()
    print(f"\n{'='*60}")
    print(f"FULL COMPETITIVE BATCH — {len(proteins)} proteins × {df['pfam_id'].nunique()} families")
    print(f"{'='*60}")

    jobs = []
    manifest_rows = []

    for _, row in proteins.iterrows():
        pfam = row["pfam_id"]
        pfam_name = row["pfam_name"]
        sp = row["species_short"]
        acc = row["uniprot_acc"]
        plen = row["protein_length"]
        domain = row["domain_of_life"]

        # read sequence from local FASTA
        fasta_path = os.path.join(seq_dir, f"{acc}.fasta")
        if not os.path.exists(fasta_path):
            print(f"  SKIP {pfam}_{sp}_{acc}: missing FASTA")
            continue

        with open(fasta_path) as f:
            seq_lines = [l.strip() for l in f if not l.startswith(">")]
        protein_seq = "".join(seq_lines)

        # token count check
        tokens = len(protein_seq) + 45
        if tokens > 5000:
            print(f"  SKIP {pfam}_{sp}_{acc}: {tokens} tokens > 5000")
            continue

        # generate competitive job
        job_name = f"competitive_{pfam}_{sp}_{acc}"
        job = makeCompetitiveJob(job_name, protein_seq, seeds=5)

        path = os.path.join(outdir, f"{job_name}.json")
        with open(path, "w") as f:
            json.dump([job], f, indent=2)

        jobs.append(job)
        manifest_rows.append({
            "job_name": job_name,
            "pfam_id": pfam,
            "pfam_name": pfam_name,
            "species": sp,
            "uniprot_acc": acc,
            "domain_of_life": domain,
            "protein_length": plen,
            "tokens": tokens,
            "binding_type": row["binding_type"],
        })

    # write batch file
    batch_path = os.path.join(outdir, "all_competitive.json")
    with open(batch_path, "w") as f:
        json.dump(jobs, f, indent=2)

    # write manifest
    mpath = os.path.join(outdir, "competitive_manifest.tsv")
    if manifest_rows:
        cols = list(manifest_rows[0].keys())
        with open(mpath, "w") as f:
            f.write("\t".join(cols) + "\n")
            for row in manifest_rows:
                f.write("\t".join(str(row[c]) for c in cols) + "\n")

    n_days = (len(jobs) + 29) // 30
    print(f"\nFULL BATCH SUMMARY")
    print(f"  jobs generated:  {len(jobs)}")
    print(f"  batch file:      {batch_path}")
    print(f"  manifest:        {mpath}")
    print(f"  AF3 Server days: {n_days} (at 30 jobs/day)")
    print("=" * 60)

    return jobs


def main():
    parser = argparse.ArgumentParser(description="generate AF3 competitive binding inputs")
    parser.add_argument("--mode", choices=["pilot", "full", "both"], default="pilot",
                        help="pilot (9 proteins), full (148 proteins), or both")
    args = parser.parse_args()

    pilot_dir = os.path.join(PROJECT, "structures", "af3_competitive_pilot")
    full_dir = os.path.join(PROJECT, "structures", "af3_competitive_full")

    if args.mode in ("pilot", "both"):
        generatePilotInputs(pilot_dir)

    if args.mode in ("full", "both"):
        generateFullInputs(full_dir)


if __name__ == "__main__":
    main()
