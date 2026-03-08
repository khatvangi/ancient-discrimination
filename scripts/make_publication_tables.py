#!/usr/bin/env python3
"""
publication tables for the ancient discrimination paper.

table 1: LUCA NA-binding domain census (summary)
table 2: NAS-Bench + ASR results for 18 LUCA families (core table)
table 3: ancestral substitution analysis by contact type
table 4: cross-tool comparison matrix (supplementary)
table 5: modern control validation (supplementary)
"""

import os
import json
import glob
import numpy as np
import pandas as pd

BASE = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
RESULTS = os.path.join(BASE, "results")
ASR = os.path.join(RESULTS, "asr")
TABDIR = os.path.join(RESULTS, "tables")
os.makedirs(TABDIR, exist_ok=True)


# ============================================================================
# TABLE 1 — LUCA NA-BINDING DOMAIN CENSUS
# ============================================================================

def make_table_1():
    """summary table of the Phase 1 census."""
    print("\n=== TABLE 1: LUCA NA-Binding Domain Census ===")

    # load intersection counts
    with open(os.path.join(RESULTS, "phase1_intersections.json")) as f:
        data = json.load(f)

    # load PDB census for co-crystal counts
    pdb = pd.read_csv(os.path.join(RESULTS, "luca_pdb_census.tsv"), sep="\t")
    n_both_pdb = (pdb["has_both"] == True).sum()

    # load nasbench for paired SI count
    nb = pd.read_csv(os.path.join(RESULTS, "nasbench_full_luca.tsv"), sep="\t")
    n_paired = len(nb[nb["generalism_mode"].notna() & (nb["generalism_mode"] != "incomplete")])

    # load convergence master for ASR count
    conv = pd.read_csv(os.path.join(ASR, "convergence_master.tsv"), sep="\t")
    n_asr = len(conv[conv["status"] == "OK"])

    # total LUCA domains (robust set from Wehbi)
    census = pd.read_csv(os.path.join(RESULTS, "phase1_census.tsv"), sep="\t")
    n_luca = len(census[census["age_class"] == "LUCA"])
    n_preluca = len(census[census["age_class"] == "preLUCA"])
    n_total = n_luca + n_preluca  # LUCA-age or older

    # RNA/DNA counts from the union upper-bound estimate
    ub = data["union_upper_bound"]
    n_rna = ub["luca_rna_binding_count"]
    n_dna = ub["luca_dna_binding_count"]

    rows = [
        ("Total LUCA-age Pfam domains", n_luca, f"{100*n_luca/(n_luca+n_preluca):.1f}"),
        ("Total pre-LUCA Pfam domains", n_preluca, f"{100*n_preluca/(n_luca+n_preluca):.1f}"),
        ("RNA-binding (LUCA)", n_rna, f"{100*n_rna/n_luca:.1f}"),
        ("DNA-binding (LUCA)", n_dna, f"{100*n_dna/n_luca:.1f}"),
        ("RNA:DNA ratio", f"{n_rna/n_dna:.1f}:1", "—"),
        ("Both RNA + DNA co-crystals available", n_both_pdb, f"{100*n_both_pdb/n_luca:.1f}"),
        ("Paired SI computed (NAS-Bench)", n_paired, f"{100*n_paired/n_luca:.1f}"),
        ("Ancestral state reconstruction completed", n_asr, f"{100*n_asr/n_luca:.1f}"),
    ]

    df = pd.DataFrame(rows, columns=["Category", "Count", "% of LUCA"])
    outpath = os.path.join(TABDIR, "table1_census.tsv")
    df.to_csv(outpath, sep="\t", index=False)
    print(f"  saved: {outpath}")
    print(df.to_string(index=False))
    return df


# ============================================================================
# TABLE 2 — NAS-BENCH + ASR RESULTS FOR 18 LUCA FAMILIES
# ============================================================================

def make_table_2():
    """core data table: SI, DI, ASR root states, ProNA2020, DRBP-EDP."""
    print("\n=== TABLE 2: NAS-Bench + ASR Results ===")

    nb = pd.read_csv(os.path.join(RESULTS, "nasbench_full_luca.tsv"), sep="\t")
    prona = pd.read_csv(os.path.join(RESULTS, "prona2020_luca_domains.tsv"), sep="\t")
    drbp = pd.read_csv(os.path.join(RESULTS, "drbp_edp_luca_domains.tsv"), sep="\t")
    conv = pd.read_csv(os.path.join(ASR, "convergence_master.tsv"), sep="\t")

    # merge prona and drbp predictions
    prona_map = prona.set_index("pfam_id")[["P_DNA", "P_RNA", "predicted_class"]].to_dict("index")
    drbp_map = drbp.set_index("pfam_id")[["stage2_pred"]].to_dict("index")
    conv_map = conv.set_index("pfam_id")[
        ["root_rna_pct", "root_dna_pct", "root_neither_pct", "ancestral_bias", "status"]
    ].to_dict("index")

    # pfam short names (from prona which has 'name' column)
    name_map = prona.set_index("pfam_id")["name"].to_dict()

    rows = []
    for _, row in nb.iterrows():
        pid = row["pfam_id"]
        if row["generalism_mode"] == "incomplete" or pd.isna(row["generalism_mode"]):
            continue

        # prona2020 predictions
        p = prona_map.get(pid, {})
        p_rna = p.get("P_RNA", np.nan)
        p_dna = p.get("P_DNA", np.nan)
        p_cls = p.get("predicted_class", "—")

        # drbp predictions
        d = drbp_map.get(pid, {})
        d_pred = d.get("stage2_pred", "—")

        # asr root states
        c = conv_map.get(pid, {})
        root_rna = c.get("root_rna_pct", np.nan)
        root_dna = c.get("root_dna_pct", np.nan)
        root_neither = c.get("root_neither_pct", np.nan)
        anc_bias = c.get("ancestral_bias", "—")
        asr_status = c.get("status", "—")

        rows.append({
            "pfam_id": pid,
            "name": name_map.get(pid, "—"),
            "age": row["age"],
            "SI_RNA": round(row["rna_SI"], 3) if pd.notna(row.get("rna_SI")) else "—",
            "SI_DNA": round(row["dna_SI"], 3) if pd.notna(row.get("dna_SI")) else "—",
            "DI": round(row["DI"], 3),
            "category": row["generalism_mode"],
            "ProNA_P_RNA": round(p_rna, 3) if pd.notna(p_rna) else "—",
            "ProNA_P_DNA": round(p_dna, 3) if pd.notna(p_dna) else "—",
            "ProNA_class": p_cls,
            "DRBP_pred": d_pred,
            "root_RNA%": round(root_rna, 1) if pd.notna(root_rna) else "—",
            "root_DNA%": round(root_dna, 1) if pd.notna(root_dna) else "—",
            "root_neither%": round(root_neither, 1) if pd.notna(root_neither) else "—",
            "ancestral_bias": anc_bias if anc_bias else "—",
        })

    df = pd.DataFrame(rows)
    df = df.sort_values("DI")

    outpath = os.path.join(TABDIR, "table2_nasbench_asr.tsv")
    df.to_csv(outpath, sep="\t", index=False)
    print(f"  saved: {outpath}")
    print(f"  {len(df)} families")
    return df


# ============================================================================
# TABLE 3 — ANCESTRAL SUBSTITUTION ANALYSIS
# ============================================================================

def make_table_3():
    """substitution classification by contact type — figure 2A in table form."""
    print("\n=== TABLE 3: Ancestral Substitution Analysis ===")

    # load all substitution data
    all_rows = []
    for f in sorted(glob.glob(os.path.join(ASR, "pf*/substitution_analysis.tsv"))):
        df = pd.read_csv(f, sep="\t")
        all_rows.append(df)
    all_df = pd.concat(all_rows, ignore_index=True)
    all_df = all_df[all_df["substitution_class"] != "absent"].copy()

    # load master for DI categories
    master = pd.read_csv(os.path.join(ASR, "convergence_master.tsv"), sep="\t")
    master = master[master["status"] == "OK"]
    di_map = master.set_index("pfam_id")["DI"].to_dict()
    all_df["DI"] = all_df["pfam_id"].map(di_map)

    def di_cat(di):
        if pd.isna(di): return "unknown"
        if di < 0.10: return "generalist"
        if di < 0.30: return "moderate"
        return "specialist"

    all_df["di_cat"] = all_df["DI"].apply(di_cat)

    contact_order = ["backbone", "base", "2oh", "base+bb", "sugar", "ALL"]
    contact_labels = ["Backbone", "Base", "2'OH", "Base+BB", "Sugar", "ALL"]

    rows = []

    for ct, ct_label in zip(contact_order, contact_labels):
        if ct == "ALL":
            subset = all_df
        else:
            subset = all_df[all_df["contact_type"] == ct]

        n = len(subset)
        if n == 0:
            continue
        n_id = (subset["substitution_class"] == "identical").sum()
        n_con = (subset["substitution_class"] == "conservative").sum()
        n_rad = (subset["substitution_class"] == "radical").sum()

        rows.append({
            "contact_type": ct_label,
            "DI_category": "all",
            "n_positions": n,
            "n_identical": n_id,
            "pct_identical": round(100 * n_id / n, 1),
            "n_conservative": n_con,
            "pct_conservative": round(100 * n_con / n, 1),
            "n_radical": n_rad,
            "pct_radical": round(100 * n_rad / n, 1),
        })

        # subrows by DI category
        for cat in ["generalist", "moderate", "specialist"]:
            cat_sub = subset[subset["di_cat"] == cat]
            nc = len(cat_sub)
            if nc == 0:
                continue
            nc_id = (cat_sub["substitution_class"] == "identical").sum()
            nc_con = (cat_sub["substitution_class"] == "conservative").sum()
            nc_rad = (cat_sub["substitution_class"] == "radical").sum()
            rows.append({
                "contact_type": f"  {ct_label}",
                "DI_category": cat,
                "n_positions": nc,
                "n_identical": nc_id,
                "pct_identical": round(100 * nc_id / nc, 1),
                "n_conservative": nc_con,
                "pct_conservative": round(100 * nc_con / nc, 1),
                "n_radical": nc_rad,
                "pct_radical": round(100 * nc_rad / nc, 1),
            })

    df = pd.DataFrame(rows)
    outpath = os.path.join(TABDIR, "table3_substitutions.tsv")
    df.to_csv(outpath, sep="\t", index=False)
    print(f"  saved: {outpath}")
    print(df.to_string(index=False))
    return df


# ============================================================================
# TABLE 4 — CROSS-TOOL COMPARISON MATRIX (SUPPLEMENTARY)
# ============================================================================

def make_table_4():
    """full prediction matrix: NAS-Bench + ProNA2020 + DRBP-EDP per family."""
    print("\n=== TABLE 4: Cross-Tool Comparison ===")

    nb = pd.read_csv(os.path.join(RESULTS, "nasbench_full_luca.tsv"), sep="\t")
    prona = pd.read_csv(os.path.join(RESULTS, "prona2020_luca_domains.tsv"), sep="\t")
    drbp = pd.read_csv(os.path.join(RESULTS, "drbp_edp_luca_domains.tsv"), sep="\t")

    # nasbench classification
    nb_class = {}
    for _, row in nb.iterrows():
        di = row["DI"]
        delta = row.get("delta_SI", 0)
        if pd.isna(di):
            nb_class[row["pfam_id"]] = "incomplete"
        elif di < 0.10:
            nb_class[row["pfam_id"]] = "dual"
        elif pd.notna(delta) and delta > 0:
            nb_class[row["pfam_id"]] = "RNA-biased"
        else:
            nb_class[row["pfam_id"]] = "DNA-biased"

    prona_map = prona.set_index("pfam_id").to_dict("index")
    drbp_map = drbp.set_index("pfam_id").to_dict("index")
    name_map = prona.set_index("pfam_id")["name"].to_dict()

    # all families from nasbench
    families = nb.sort_values("DI")["pfam_id"].tolist()

    rows = []
    for pid in families:
        p = prona_map.get(pid, {})
        d = drbp_map.get(pid, {})
        nbr = nb[nb["pfam_id"] == pid]

        row_data = {
            "pfam_id": pid,
            "name": name_map.get(pid, "—"),
            "DI": round(nbr.iloc[0]["DI"], 3) if len(nbr) > 0 and pd.notna(nbr.iloc[0]["DI"]) else "—",
            "NASBench_class": nb_class.get(pid, "—"),
            "NASBench_SI_RNA": round(nbr.iloc[0]["rna_SI"], 3) if len(nbr) > 0 and pd.notna(nbr.iloc[0].get("rna_SI")) else "—",
            "NASBench_SI_DNA": round(nbr.iloc[0]["dna_SI"], 3) if len(nbr) > 0 and pd.notna(nbr.iloc[0].get("dna_SI")) else "—",
            "ProNA_P_RNA": round(p.get("P_RNA", np.nan), 3) if pd.notna(p.get("P_RNA", np.nan)) else "—",
            "ProNA_P_DNA": round(p.get("P_DNA", np.nan), 3) if pd.notna(p.get("P_DNA", np.nan)) else "—",
            "ProNA_class": p.get("predicted_class", "—"),
            "DRBP_stage1": d.get("stage1_pred", "—"),
            "DRBP_stage2": d.get("stage2_pred", "—"),
        }

        # agreement check
        classes = []
        nb_c = nb_class.get(pid, "")
        if "RNA" in nb_c: classes.append("RNA")
        elif "DNA" in nb_c: classes.append("DNA")
        elif nb_c == "dual": classes.append("dual")

        p_c = str(p.get("predicted_class", ""))
        if "RNA" in p_c: classes.append("RNA")
        elif "DNA" in p_c: classes.append("DNA")
        elif "dual" in p_c: classes.append("dual")

        d_c = str(d.get("stage2_pred", ""))
        if "RNA" in d_c: classes.append("RNA")
        elif "DNA" in d_c: classes.append("DNA")

        if len(classes) >= 2:
            row_data["agreement"] = "agree" if len(set(classes)) == 1 else "disagree"
        else:
            row_data["agreement"] = "—"

        rows.append(row_data)

    df = pd.DataFrame(rows)
    outpath = os.path.join(TABDIR, "table4_cross_tool.tsv")
    df.to_csv(outpath, sep="\t", index=False)
    print(f"  saved: {outpath}")
    print(f"  {len(df)} families, agreement: {(df['agreement']=='agree').sum()}/{(df['agreement']!='—').sum()}")
    return df


# ============================================================================
# TABLE 5 — MODERN CONTROL VALIDATION (SUPPLEMENTARY)
# ============================================================================

def make_table_5():
    """modern RNA/DNA specialists with SI values and known biology."""
    print("\n=== TABLE 5: Modern Control Validation ===")

    modern = pd.read_csv(os.path.join(RESULTS, "nasbench_modern_controls.tsv"), sep="\t")

    rows = []
    for _, row in modern.iterrows():
        si = row.get("SI", np.nan)
        rows.append({
            "name": row["name"],
            "PDB": row["pdb_id"],
            "category": row["category"],
            "target_NA": row["target_na"],
            "description": row["description"],
            "SI": round(si, 3) if pd.notna(si) else "FAIL",
            "bb_frac": round(row["bb_frac"], 3) if pd.notna(row.get("bb_frac")) else "—",
            "base_frac": round(row["base_frac"], 3) if pd.notna(row.get("base_frac")) else "—",
            "2oh_frac": round(row["oh_frac"], 3) if pd.notna(row.get("oh_frac")) else "—",
            "n_contacts": int(row["n_contacts"]) if pd.notna(row.get("n_contacts")) and row.get("n_contacts", 0) > 0 else "—",
        })

    df = pd.DataFrame(rows)
    outpath = os.path.join(TABDIR, "table5_modern_controls.tsv")
    df.to_csv(outpath, sep="\t", index=False)
    print(f"  saved: {outpath}")
    print(df.to_string(index=False))
    return df


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("generating publication tables...")
    print(f"output: {TABDIR}/\n")

    for name, func in [
        ("Table 1", make_table_1),
        ("Table 2", make_table_2),
        ("Table 3", make_table_3),
        ("Table 4", make_table_4),
        ("Table 5", make_table_5),
    ]:
        try:
            func()
        except Exception as e:
            print(f"  ERROR generating {name}: {e}")
            import traceback
            traceback.print_exc()

    print("\nDone.")


if __name__ == "__main__":
    main()
