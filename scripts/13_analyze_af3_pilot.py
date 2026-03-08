#!/usr/bin/env python3
"""analyze AF3 pilot results: 6 jobs (3 proteins × RNA + DNA).

extracts ipTM, pTM, and ranking_score from summary_confidences.json,
computes the Discrimination Index (DI) for each protein, and evaluates
whether pilot results justify launching the full 296-job batch.

pilot families:
  - S4  (PF01479, bsu P37557): RNA-annotated ribosomal protein
  - HHH (PF00633, eco P0AB83): DNA-annotated but predicted generalist
  - S1  (PF00575, mja Q57840): RNA-annotated, known dual-binder

Discrimination Index = ipTM(annotated substrate) - ipTM(other substrate)
  positive → prefers annotated substrate
  near zero → generalist (thesis prediction for ancient domains)

usage:
  python scripts/13_analyze_af3_pilot.py
"""

import json
import os
import sys

PROJECT = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
OUTPUT_DIR = os.path.join(PROJECT, "structures", "af3_outputs")
RESULTS_DIR = os.path.join(PROJECT, "results")

# pilot job definitions: (pfam, species, uniprot, annotated_substrate)
PILOT_JOBS = [
    ("PF01479", "bsu", "P37557", "RNA", "S4"),
    ("PF00633", "eco", "P0AB83", "DNA", "HHH"),
    ("PF00575", "mja", "Q57840", "RNA", "S1"),
]


def loadConfidences(job_name):
    """load summary_confidences.json for a given job name."""
    conf_dir = os.path.join(OUTPUT_DIR, job_name)
    conf_file = os.path.join(conf_dir, f"{job_name}_summary_confidences.json")
    if not os.path.exists(conf_file):
        return None
    with open(conf_file) as f:
        return json.load(f)


def main():
    print("=" * 60)
    print("AF3 PILOT RESULTS ANALYSIS")
    print("=" * 60)

    results = {}  # key: (pfam, species, uniprot), value: {RNA: metrics, DNA: metrics}
    missing = []

    for pfam, sp, uni, annotated, name in PILOT_JOBS:
        results[(pfam, sp, uni)] = {"name": name, "annotated": annotated}

        for substrate in ["RNA", "DNA"]:
            job_name = f"{pfam}_{sp}_{uni}_{substrate}"
            conf = loadConfidences(job_name)

            if conf is None:
                missing.append(job_name)
                results[(pfam, sp, uni)][substrate] = None
                continue

            metrics = {
                "iptm": conf.get("iptm", 0),
                "ptm": conf.get("ptm", 0),
                "ranking_score": conf.get("ranking_score", 0),
            }

            # also extract per-chain pLDDT if available
            chain_pair_iptm = conf.get("chain_pair_iptm", {})
            if chain_pair_iptm:
                metrics["chain_pair_iptm"] = chain_pair_iptm

            results[(pfam, sp, uni)][substrate] = metrics

    # report missing jobs
    if missing:
        print(f"\nWARNING: {len(missing)} jobs not yet complete:")
        for m in missing:
            print(f"  - {m}")
        print()

    # display results table
    print(f"\n{'Protein':<10} {'Pfam':<10} {'Substrate':<10} {'ipTM':>6} {'pTM':>6} {'Ranking':>8}")
    print("-" * 60)

    for (pfam, sp, uni), data in results.items():
        for substrate in ["RNA", "DNA"]:
            if data.get(substrate) is None:
                print(f"{data['name']:<10} {pfam:<10} {substrate:<10} {'N/A':>6} {'N/A':>6} {'N/A':>8}")
            else:
                m = data[substrate]
                print(f"{data['name']:<10} {pfam:<10} {substrate:<10} {m['iptm']:>6.3f} {m['ptm']:>6.3f} {m['ranking_score']:>8.3f}")

    # compute Discrimination Index
    print(f"\n{'='*60}")
    print("DISCRIMINATION INDEX (DI)")
    print(f"{'='*60}")
    print(f"\n{'Protein':<10} {'Annotated':<10} {'ipTM(ann)':>10} {'ipTM(other)':>12} {'DI':>8} {'Interpretation'}")
    print("-" * 75)

    pilot_summary = []

    for (pfam, sp, uni), data in results.items():
        ann = data["annotated"]
        other = "DNA" if ann == "RNA" else "RNA"

        if data.get(ann) is None or data.get(other) is None:
            print(f"{data['name']:<10} {ann:<10} {'N/A':>10} {'N/A':>12} {'N/A':>8}")
            continue

        iptm_ann = data[ann]["iptm"]
        iptm_other = data[other]["iptm"]
        di = iptm_ann - iptm_other

        # interpretation
        if abs(di) < 0.05:
            interp = "GENERALIST (|DI| < 0.05)"
        elif di > 0.15:
            interp = f"SPECIALIST for {ann} (DI > 0.15)"
        elif di > 0.05:
            interp = f"mild preference for {ann}"
        elif di < -0.15:
            interp = f"UNEXPECTED: prefers {other} over annotated {ann}"
        else:
            interp = f"mild preference for {other} (unexpected)"

        print(f"{data['name']:<10} {ann:<10} {iptm_ann:>10.3f} {iptm_other:>12.3f} {di:>8.3f}   {interp}")

        pilot_summary.append({
            "protein": data["name"],
            "pfam": pfam,
            "species": sp,
            "uniprot": uni,
            "annotated_substrate": ann,
            "iptm_annotated": iptm_ann,
            "iptm_other": iptm_other,
            "ptm_annotated": data[ann]["ptm"],
            "ptm_other": data[other]["ptm"],
            "discrimination_index": di,
            "interpretation": interp,
        })

    # thesis evaluation
    print(f"\n{'='*60}")
    print("THESIS EVALUATION")
    print(f"{'='*60}")

    if pilot_summary:
        # check HHH specifically — thesis predicts generalist
        hhh = [p for p in pilot_summary if p["protein"] == "HHH"]
        if hhh:
            h = hhh[0]
            print(f"\nHHH (key test): DI = {h['discrimination_index']:.3f}")
            if abs(h["discrimination_index"]) < 0.10:
                print("  → SUPPORTS thesis: HHH shows generalist NA binding")
            else:
                print(f"  → DOES NOT support thesis: HHH shows {h['interpretation']}")

        # check S1 — known dual binder, should show low DI
        s1 = [p for p in pilot_summary if p["protein"] == "S1"]
        if s1:
            s = s1[0]
            print(f"\nS1 (positive control): DI = {s['discrimination_index']:.3f}")
            if abs(s["discrimination_index"]) < 0.15:
                print("  → EXPECTED: S1 (known dual-binder) shows low discrimination")
            else:
                print(f"  → SURPRISING: S1 shows strong preference ({s['interpretation']})")

        # check overall ipTM quality
        all_iptm = [p["iptm_annotated"] for p in pilot_summary] + [p["iptm_other"] for p in pilot_summary]
        min_iptm = min(all_iptm)
        max_iptm = max(all_iptm)
        print(f"\nipTM range: {min_iptm:.3f} - {max_iptm:.3f}")
        if min_iptm < 0.3:
            print("  WARNING: some ipTM values very low — predictions may be unreliable")
            print("  consider: are the 10mer substrates too short?")
        elif min_iptm < 0.5:
            print("  NOTE: moderate ipTM values — predictions are usable but not high-confidence")
        else:
            print("  GOOD: all ipTM values > 0.5 — predictions are reliable")

    # go/no-go decision
    print(f"\n{'='*60}")
    print("GO/NO-GO DECISION")
    print(f"{'='*60}")

    if not pilot_summary:
        print("\nCANNOT DECIDE: no pilot results available yet")
        print("run this script again after pilot completes")
    else:
        # sanity check: annotated ipTM > 0.3 for at least 2/3 proteins
        good_ann = sum(1 for p in pilot_summary if p["iptm_annotated"] > 0.3)
        if good_ann >= 2:
            print(f"\nGO: {good_ann}/3 proteins have annotated ipTM > 0.3")
            print("proceed with full 296-job batch")
        else:
            print(f"\nCAUTION: only {good_ann}/3 proteins have annotated ipTM > 0.3")
            print("consider adjusting substrate length or using different substrates")

    # save summary JSON
    os.makedirs(RESULTS_DIR, exist_ok=True)
    summary_path = os.path.join(RESULTS_DIR, "pilot_af3_summary.json")
    with open(summary_path, "w") as f:
        json.dump({
            "pilot_jobs": PILOT_JOBS,
            "results": pilot_summary,
            "missing": missing,
        }, f, indent=2, default=str)
    print(f"\nsummary saved: {summary_path}")


if __name__ == "__main__":
    main()
