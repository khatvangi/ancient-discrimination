#!/usr/bin/env python3
"""generate LOCAL AF3 format competitive binding inputs.

local AF3 (Docker-based, HMMER for MSA) uses dialect "alphafold3"
with protein/rna/dna keys and explicit chain id fields.

generates inputs for the 9-protein pilot and optionally the full batch.

usage:
  python scripts/14b_generate_local_competitive_inputs.py
"""

import json
import os
import sys
import time

import requests

PROJECT = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
OUTDIR = os.path.join(PROJECT, "structures", "af3_competitive_local")

# nucleic acid probes — same as server version
RNA_PROBE = "GACUGAUUCGAUCAG"       # 15-nt ssRNA
DNA_SENSE = "GACTGATTCGATCAG"       # 15-nt sense
DNA_ANTISENSE = "CTGATCGAATCAGTC"   # 15-nt antisense (reverse complement)

# pilot proteins — same definitions as 14_generate_competitive_inputs.py
PILOT_PROTEINS = [
    {"id": "CspA_Ecoli",        "uniprot": "P0A9X9",     "family": "CSD",  "crop": None},
    {"id": "CspB_Hfx",          "uniprot": "A0A0D6JQU8", "family": "CSD",  "crop": None},
    {"id": "YBX1_Human",        "uniprot": "P67809",     "family": "CSD",  "crop": (52, 135)},
    {"id": "Sso10b_Ssolf",      "uniprot": "P60849",     "family": "Alba", "crop": None},
    {"id": "Alba_Mjann",        "uniprot": "Q58620",     "family": "Alba", "crop": None},
    {"id": "Rpp25_Human",       "uniprot": "Q9BUL9",     "family": "Alba", "crop": None},
    {"id": "NusA_KH_Ecoli",     "uniprot": "P0AFF6",     "family": "KH",   "crop": (195, 350)},
    {"id": "NusA_Afulg",        "uniprot": "O28388",     "family": "KH",   "crop": None},
    {"id": "PCBP1_KH1_Human",   "uniprot": "Q15365",     "family": "KH",   "crop": (10, 82)},
]


def fetchUniprotSeq(accession, max_retries=3):
    """fetch protein sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            lines = r.text.strip().split("\n")
            return "".join(l.strip() for l in lines if not l.startswith(">"))
        except requests.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(2)
            else:
                raise RuntimeError(f"failed to fetch {accession}: {e}")


def makeLocalCompetitiveJob(name, protein_seq):
    """build one LOCAL AF3 competitive binding JSON.

    local format: dialect "alphafold3", chains have "id" fields.
    chain A = protein, B = RNA, C = DNA sense, D = DNA antisense.
    one seed, 5 diffusion samples (default).
    """
    return {
        "name": name,
        "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": ["A"], "sequence": protein_seq}},
            {"rna":     {"id": ["B"], "sequence": RNA_PROBE}},
            {"dna":     {"id": ["C"], "sequence": DNA_SENSE}},
            {"dna":     {"id": ["D"], "sequence": DNA_ANTISENSE}},
        ],
        "dialect": "alphafold3",
        "version": 1,
    }


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print(f"generating local AF3 competitive inputs → {OUTDIR}/")
    print(f"RNA:  {RNA_PROBE}")
    print(f"DNA+: {DNA_SENSE}")
    print(f"DNA-: {DNA_ANTISENSE}")
    print()

    jobs = []
    for p in PILOT_PROTEINS:
        name = p["id"]
        seq = fetchUniprotSeq(p["uniprot"])

        if p["crop"]:
            s, e = p["crop"]
            seq = seq[s - 1 : e]
            crop_info = f"crop {s}-{e}"
        else:
            crop_info = "full"

        job_name = f"competitive_{name}"
        job = makeLocalCompetitiveJob(job_name, seq)

        # write individual JSON (single job, NOT wrapped in list)
        path = os.path.join(OUTDIR, f"{job_name}.json")
        with open(path, "w") as f:
            json.dump(job, f, indent=2)

        jobs.append(job_name)
        tokens = len(seq) + 45
        print(f"  {job_name:45s}  {len(seq):4d} aa  ({crop_info})  {tokens} tokens")

    print(f"\n{len(jobs)} JSON files written to {OUTDIR}/")


if __name__ == "__main__":
    main()
