#!/usr/bin/env python3
"""
conservative substitution analysis of "neither" ancestral residues.

for every NA-contacting alignment column across all families, classifies
the ancestral-to-modern residue relationship as:
  - identical:    same residue
  - conservative: same physicochemical group (preserves contact capacity)
  - radical:      different physicochemical group (contact type changed)
  - absent:       no modern residue to compare (gap/missing)

physicochemical groups (from user specification):
  acidic H-bond donors:   {D, E}
  basic H-bond donors:    {K, R}
  polar/H-bond:           {S, T, N, Q}
  aromatic/stacking:       {F, Y, W, H}
  small hydrophobic:       {A, V, I, L, M}
  special:                 {G, P, C}

produces:
  results/asr/{family}/substitution_analysis.tsv  (per-family)
  results/asr/neither_decomposition.tsv           (aggregate)
"""

import os
import glob
import pandas as pd
import numpy as np

BASE = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/results/asr"
MASTER = os.path.join(BASE, "convergence_master.tsv")

# ── physicochemical groups ──────────────────────────────────────────────────

PHYSICO_GROUPS = {
    "D": "acidic", "E": "acidic",
    "K": "basic",  "R": "basic",
    "S": "polar",  "T": "polar",  "N": "polar",  "Q": "polar",
    "F": "aromatic", "Y": "aromatic", "W": "aromatic", "H": "aromatic",
    "A": "hydrophobic", "V": "hydrophobic", "I": "hydrophobic",
    "L": "hydrophobic", "M": "hydrophobic",
    "G": "special", "P": "special", "C": "special",
}

# highlight families for detailed output
HIGHLIGHT = ["PF00013", "PF00270", "PF03129"]


def classify_substitution(ancestral, modern):
    """classify the relationship between ancestral and modern residues."""
    if not ancestral or not modern or ancestral == "nan" or modern == "nan":
        return "absent"
    ancestral = ancestral.strip().upper()
    modern = modern.strip().upper()
    if not ancestral or not modern or ancestral == "-" or modern == "-":
        return "absent"
    if ancestral == modern:
        return "identical"
    anc_group = PHYSICO_GROUPS.get(ancestral, "unknown")
    mod_group = PHYSICO_GROUPS.get(modern, "unknown")
    if anc_group == "unknown" or mod_group == "unknown":
        return "absent"
    if anc_group == mod_group:
        return "conservative"
    return "radical"


def simplify_contact_type(classification):
    """simplify contact classification to broad categories."""
    if not classification or classification == "nan" or str(classification).strip() == "":
        return "none"
    c = str(classification).lower()
    if "2oh" in c:
        return "2oh"
    elif "base" in c and "backbone" in c:
        return "base+bb"
    elif "base" in c:
        return "base"
    elif "backbone" in c:
        return "backbone"
    elif "sugar" in c:
        return "sugar"
    return "other"


def analyze_family(pfam_id):
    """analyze substitution patterns for a single family."""
    fam_dir = os.path.join(BASE, pfam_id.lower())
    conv_file = os.path.join(fam_dir, f"{pfam_id}_convergence_summary.tsv")

    if not os.path.exists(conv_file):
        return None

    try:
        df = pd.read_csv(conv_file, sep="\t")
    except Exception:
        return None

    if len(df) == 0:
        return None

    results = []
    for _, row in df.iterrows():
        col_type = str(row.get("col_type", "")).strip()
        rna_aa = str(row.get("rna_aa", "")).strip()
        dna_aa = str(row.get("dna_aa", "")).strip()
        rna_class = str(row.get("rna_classification", "")).strip()
        dna_class = str(row.get("dna_classification", "")).strip()
        anc_state = str(row.get("Node1_state", "")).strip()
        anc_pp = float(row.get("Node1_pp", 0))

        # handle nan strings
        if rna_aa == "nan": rna_aa = ""
        if dna_aa == "nan": dna_aa = ""
        if rna_class == "nan": rna_class = ""
        if dna_class == "nan": dna_class = ""
        if anc_state == "nan": anc_state = ""

        # for each column, we need to check against the relevant modern residue
        # and determine what "match" means
        comparisons = []

        if col_type in ("RNA-only", "shared") and rna_aa:
            comparisons.append({
                "target_type": "RNA",
                "modern_aa": rna_aa,
                "contact_type": simplify_contact_type(rna_class),
            })

        if col_type in ("DNA-only", "shared") and dna_aa:
            comparisons.append({
                "target_type": "DNA",
                "modern_aa": dna_aa,
                "contact_type": simplify_contact_type(dna_class),
            })

        for comp in comparisons:
            sub_class = classify_substitution(anc_state, comp["modern_aa"])

            results.append({
                "pfam_id": pfam_id,
                "aln_col": row["aln_col"],
                "col_type": col_type,
                "target_type": comp["target_type"],
                "modern_aa": comp["modern_aa"],
                "ancestral_aa": anc_state,
                "ancestral_pp": anc_pp,
                "contact_type": comp["contact_type"],
                "substitution_class": sub_class,
                "anc_group": PHYSICO_GROUPS.get(anc_state.upper(), "?") if anc_state else "?",
                "mod_group": PHYSICO_GROUPS.get(comp["modern_aa"].upper(), "?") if comp["modern_aa"] else "?",
            })

    return pd.DataFrame(results) if results else None


def print_family_detail(fam_df, pfam_id):
    """print detailed per-position analysis for a highlighted family."""
    print(f"\n{'='*80}")
    print(f"DETAILED: {pfam_id}")
    print(f"{'='*80}")
    print(f"\n{'Col':>4} {'Type':<8} {'Tgt':<4} {'Modern':<8} {'Ancestor':<10} {'PP':>5} "
          f"{'Contact':<10} {'SubClass':<13} {'Groups'}")
    print("-" * 95)

    for _, r in fam_df.sort_values("aln_col").iterrows():
        # mark "neither" rows that are radical — these are the key ones
        marker = ""
        if r["substitution_class"] == "radical":
            marker = " ★"
        elif r["substitution_class"] == "conservative":
            marker = " ~"

        print(f"{r['aln_col']:>4} {r['col_type']:<8} {r['target_type']:<4} "
              f"{r['modern_aa']:<8} {r['ancestral_aa']:<10} {r['ancestral_pp']:>5.2f} "
              f"{r['contact_type']:<10} {r['substitution_class']:<13} "
              f"{r['anc_group']}→{r['mod_group']}{marker}")


def make_decomposition_table(all_df, label="ALL FAMILIES"):
    """create the neither decomposition table by contact type."""
    # filter to "neither" positions: where substitution is conservative or radical
    # (identical = "matches", absent = skip)
    neither_df = all_df[all_df["substitution_class"].isin(["conservative", "radical"])]

    print(f"\n{'='*80}")
    print(f"NEITHER DECOMPOSITION: {label}")
    print(f"{'='*80}")

    # overall breakdown first
    total = len(all_df[all_df["substitution_class"] != "absent"])
    n_identical = (all_df["substitution_class"] == "identical").sum()
    n_conservative = (all_df["substitution_class"] == "conservative").sum()
    n_radical = (all_df["substitution_class"] == "radical").sum()

    print(f"\noverall (n={total}):")
    print(f"  identical:    {n_identical:>4} ({100*n_identical/total:.1f}%)")
    print(f"  conservative: {n_conservative:>4} ({100*n_conservative/total:.1f}%)")
    print(f"  radical:      {n_radical:>4} ({100*n_radical/total:.1f}%)")
    print(f"  → 'neither' = conservative + radical = "
          f"{n_conservative+n_radical} ({100*(n_conservative+n_radical)/total:.1f}%)")

    # by contact type
    print(f"\n{'Contact type':<15} {'Identical':>10} {'Conservative':>14} {'Radical':>10} "
          f"{'Total':>8} {'%Radical':>10}")
    print("-" * 75)

    contact_types = ["backbone", "base", "2oh", "base+bb", "sugar"]
    for ct in contact_types:
        ct_df = all_df[(all_df["contact_type"] == ct) & (all_df["substitution_class"] != "absent")]
        if len(ct_df) == 0:
            continue
        n_id = (ct_df["substitution_class"] == "identical").sum()
        n_con = (ct_df["substitution_class"] == "conservative").sum()
        n_rad = (ct_df["substitution_class"] == "radical").sum()
        n_tot = len(ct_df)
        pct_rad = 100 * n_rad / n_tot if n_tot > 0 else 0

        print(f"{ct:<15} {n_id:>10} ({100*n_id/n_tot:>4.0f}%) "
              f"{n_con:>7} ({100*n_con/n_tot:>4.0f}%) "
              f"{n_rad:>5} ({100*n_rad/n_tot:>4.0f}%) "
              f"{n_tot:>6} {pct_rad:>9.1f}%")

    # total row
    all_valid = all_df[all_df["substitution_class"] != "absent"]
    n_id = (all_valid["substitution_class"] == "identical").sum()
    n_con = (all_valid["substitution_class"] == "conservative").sum()
    n_rad = (all_valid["substitution_class"] == "radical").sum()
    n_tot = len(all_valid)
    print("-" * 75)
    print(f"{'ALL':<15} {n_id:>10} ({100*n_id/n_tot:>4.0f}%) "
          f"{n_con:>7} ({100*n_con/n_tot:>4.0f}%) "
          f"{n_rad:>5} ({100*n_rad/n_tot:>4.0f}%) "
          f"{n_tot:>6} {100*n_rad/n_tot:>9.1f}%")

    # by target type (RNA vs DNA comparisons)
    print(f"\nby target type:")
    for tt in ["RNA", "DNA"]:
        tt_df = all_valid[all_valid["target_type"] == tt]
        if len(tt_df) == 0:
            continue
        n_id = (tt_df["substitution_class"] == "identical").sum()
        n_con = (tt_df["substitution_class"] == "conservative").sum()
        n_rad = (tt_df["substitution_class"] == "radical").sum()
        n_tot = len(tt_df)
        print(f"  {tt}: identical={100*n_id/n_tot:.1f}%, "
              f"conservative={100*n_con/n_tot:.1f}%, "
              f"radical={100*n_rad/n_tot:.1f}% (n={n_tot})")

    return all_valid


def main():
    # load master to get family list
    master = pd.read_csv(MASTER, sep="\t")
    ok_families = master[master["status"] == "OK"]["pfam_id"].tolist()

    print(f"analyzing {len(ok_families)} families...\n")

    # analyze all families
    all_results = []
    for pfam_id in sorted(ok_families):
        fam_df = analyze_family(pfam_id)
        if fam_df is not None and len(fam_df) > 0:
            all_results.append(fam_df)

            # save per-family file
            out_path = os.path.join(BASE, pfam_id.lower(), "substitution_analysis.tsv")
            fam_df.to_csv(out_path, sep="\t", index=False)

            n_valid = len(fam_df[fam_df["substitution_class"] != "absent"])
            n_id = (fam_df["substitution_class"] == "identical").sum()
            n_con = (fam_df["substitution_class"] == "conservative").sum()
            n_rad = (fam_df["substitution_class"] == "radical").sum()
            di = master[master["pfam_id"] == pfam_id]["DI"].values[0]

            print(f"  {pfam_id} (DI={di:.3f}): {n_valid} comparisons → "
                  f"identical={n_id}, conservative={n_con}, radical={n_rad}")

    if not all_results:
        print("ERROR: no results!")
        return

    all_df = pd.concat(all_results, ignore_index=True)

    # ── detailed output for highlighted families ────────────────────────────

    for pfam_id in HIGHLIGHT:
        fam_df = all_df[all_df["pfam_id"] == pfam_id]
        if len(fam_df) > 0:
            print_family_detail(fam_df, pfam_id)

    # ── aggregate decomposition tables ──────────────────────────────────────

    # all families combined
    make_decomposition_table(all_df, "ALL 17 FAMILIES")

    # just the 3 highlighted families
    highlight_df = all_df[all_df["pfam_id"].isin(HIGHLIGHT)]
    make_decomposition_table(highlight_df, "3 SELECTED FAMILIES (PF00013 + PF00270 + PF03129)")

    # by DI category
    for cat_name, di_min, di_max in [
        ("generalists (DI<0.10)", 0, 0.10),
        ("moderates (0.10-0.30)", 0.10, 0.30),
        ("specialists (DI>0.30)", 0.30, 1.0),
    ]:
        cat_fams = master[(master["DI"] >= di_min) & (master["DI"] < di_max) &
                          (master["status"] == "OK")]["pfam_id"].tolist()
        cat_df = all_df[all_df["pfam_id"].isin(cat_fams)]
        if len(cat_df) > 0:
            make_decomposition_table(cat_df, f"{cat_name} (n={len(cat_fams)} families)")

    # ── save aggregate ──────────────────────────────────────────────────────

    out_path = os.path.join(BASE, "neither_decomposition.tsv")
    # build summary table
    summary_rows = []
    for ct in ["backbone", "base", "2oh", "base+bb", "sugar", "ALL"]:
        if ct == "ALL":
            ct_df = all_df[all_df["substitution_class"] != "absent"]
        else:
            ct_df = all_df[(all_df["contact_type"] == ct) &
                          (all_df["substitution_class"] != "absent")]
        if len(ct_df) == 0:
            continue
        n_id = (ct_df["substitution_class"] == "identical").sum()
        n_con = (ct_df["substitution_class"] == "conservative").sum()
        n_rad = (ct_df["substitution_class"] == "radical").sum()
        n_tot = len(ct_df)
        summary_rows.append({
            "contact_type": ct,
            "n_total": n_tot,
            "n_identical": n_id,
            "pct_identical": round(100 * n_id / n_tot, 1),
            "n_conservative": n_con,
            "pct_conservative": round(100 * n_con / n_tot, 1),
            "n_radical": n_rad,
            "pct_radical": round(100 * n_rad / n_tot, 1),
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(out_path, sep="\t", index=False)
    print(f"\nsaved: {out_path}")

    # ── key interpretation ──────────────────────────────────────────────────

    all_valid = all_df[all_df["substitution_class"] != "absent"]
    n_tot = len(all_valid)
    n_rad = (all_valid["substitution_class"] == "radical").sum()
    n_con = (all_valid["substitution_class"] == "conservative").sum()
    n_id = (all_valid["substitution_class"] == "identical").sum()

    # 2'OH specific
    oh2 = all_valid[all_valid["contact_type"] == "2oh"]
    oh2_rad = (oh2["substitution_class"] == "radical").sum()
    bb = all_valid[all_valid["contact_type"] == "backbone"]
    bb_rad = (bb["substitution_class"] == "radical").sum()

    print("\n" + "=" * 80)
    print("KEY INTERPRETATION")
    print("=" * 80)

    pct_rad = 100 * n_rad / n_tot
    pct_con = 100 * n_con / n_tot
    pct_id = 100 * n_id / n_tot

    if pct_rad > pct_con:
        print(f"\n→ RADICAL dominates ({pct_rad:.0f}% radical vs {pct_con:.0f}% conservative)")
        print("  the 'neither' positions represent GENUINELY DIFFERENT residues.")
        print("  the ancestral binding mode is extinct — modern specialists invented")
        print("  new contact residues de novo.")
    elif pct_con > pct_rad:
        print(f"\n→ CONSERVATIVE dominates ({pct_con:.0f}% conservative vs {pct_rad:.0f}% radical)")
        print("  the 'neither' count is inflated by conservative substitutions.")
        print("  the ancestor used a recognizable variant of the modern binding mode.")
    else:
        print(f"\n→ MIXED ({pct_rad:.0f}% radical, {pct_con:.0f}% conservative)")

    if len(oh2) > 0 and len(bb) > 0:
        oh2_pct_rad = 100 * oh2_rad / len(oh2) if len(oh2) > 0 else 0
        bb_pct_rad = 100 * bb_rad / len(bb) if len(bb) > 0 else 0
        print(f"\n  2'OH contacts:     {oh2_pct_rad:.0f}% radical (n={len(oh2)})")
        print(f"  backbone contacts: {bb_pct_rad:.0f}% radical (n={len(bb)})")

        if oh2_pct_rad > bb_pct_rad + 5:
            print(f"\n  ★ 2'OH contacts are MORE radical than backbone contacts.")
            print(f"    → ancestor preserved nonspecific binding (backbone)")
            print(f"      but LACKED discriminating contacts (2'OH, groove reading)")
            print(f"    → strongest evidence for 'generalism = absence of specificity'")
        elif bb_pct_rad > oh2_pct_rad + 5:
            print(f"\n  ★ backbone contacts are MORE radical than 2'OH contacts.")
            print(f"    → unexpected: ancestor diverged in structural contacts")
        else:
            print(f"\n  ~ 2'OH and backbone contacts show similar radical rates")

    print(f"\n  overall: {pct_id:.0f}% identical, {pct_con:.0f}% conservative, {pct_rad:.0f}% radical")
    print(f"  'neither' = {pct_con + pct_rad:.0f}% of all contact positions")
    print(f"  → of the 'neither' positions, {100*n_rad/(n_rad+n_con):.0f}% are radical changes")


if __name__ == "__main__":
    main()
