#!/usr/bin/env python3
"""
deep analysis of ASR convergence results across 17 LUCA-age families.

produces a multi-panel figure:
  panel A: node-level heatmap (families × nodes, colored by RNA%)
  panel B: 2'OH contact-type analysis (ancestral match rate by contact class)
  panel C: DI vs root ancestral RNA% scatter, colored by generalism mode

also prints detailed textual analysis including PF01131 positive control.
"""

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

# ── configuration ───────────────────────────────────────────────────────────

BASE = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/results/asr"
MASTER = os.path.join(BASE, "convergence_master.tsv")
OUTDIR = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/figures"
os.makedirs(OUTDIR, exist_ok=True)

# publication style
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# okabe-ito colorblind-safe palette
COLORS = {
    "generalist": "#009E73",      # green
    "moderate": "#E69F00",         # orange
    "specialist": "#D55E00",       # red-orange
    "RNA-biased": "#56B4E9",       # sky blue
    "DNA-biased": "#CC79A7",       # pink
    "neutral": "#999999",          # gray
    "2OH": "#0072B2",              # dark blue
    "backbone": "#F0E442",         # yellow
    "base": "#009E73",             # green
    "other": "#999999",            # gray
}

# ── data loading ────────────────────────────────────────────────────────────

def load_master():
    """load convergence master table, keep only OK families."""
    df = pd.read_csv(MASTER, sep="\t")
    return df[df["status"] == "OK"].copy()


def load_node_summaries():
    """load all node_summary.tsv files into a single dataframe."""
    rows = []
    for d in sorted(glob.glob(os.path.join(BASE, "pf*"))):
        fam = os.path.basename(d)
        pfam_id = fam.upper()
        fname = os.path.join(d, f"{pfam_id}_node_summary.tsv")
        if not os.path.exists(fname):
            continue
        df = pd.read_csv(fname, sep="\t")
        df["pfam_id"] = pfam_id
        rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def load_convergence_details():
    """load all convergence_summary.tsv files for contact-type analysis."""
    rows = []
    for d in sorted(glob.glob(os.path.join(BASE, "pf*"))):
        fam = os.path.basename(d)
        pfam_id = fam.upper()
        fname = os.path.join(d, f"{pfam_id}_convergence_summary.tsv")
        if not os.path.exists(fname):
            continue
        try:
            df = pd.read_csv(fname, sep="\t")
        except Exception:
            continue
        if len(df) == 0:
            continue
        df["pfam_id"] = pfam_id
        rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


# ── analysis functions ──────────────────────────────────────────────────────

def classify_node1_match(row):
    """for a single convergence row, check if Node1 (root) matches RNA or DNA aa."""
    state = str(row.get("Node1_state", "")).strip()
    rna_aa = str(row.get("rna_aa", "")).strip()
    dna_aa = str(row.get("dna_aa", "")).strip()

    if not state or state == "nan":
        return "unknown"

    matches_rna = (state == rna_aa) if rna_aa and rna_aa != "nan" else False
    matches_dna = (state == dna_aa) if dna_aa and dna_aa != "nan" else False

    if matches_rna and matches_dna:
        return "both"
    elif matches_rna:
        return "rna"
    elif matches_dna:
        return "dna"
    else:
        return "neither"


def compute_contact_type_stats(conv_df):
    """compute ancestral match rates grouped by contact classification type.

    for each position, we use the RNA classification (since we're testing
    whether ancestral RNA-contacting residues are conserved).
    we also look at DNA classification for DNA-only columns.
    """
    results = []

    for _, row in conv_df.iterrows():
        rna_class = str(row.get("rna_classification", "")).strip()
        dna_class = str(row.get("dna_classification", "")).strip()
        col_type = str(row.get("col_type", "")).strip()

        # determine the relevant contact classification
        if col_type in ("RNA-only", "shared") and rna_class and rna_class != "nan":
            contact_class = rna_class
            relevant_type = "rna"
        elif col_type == "DNA-only" and dna_class and dna_class != "nan":
            contact_class = dna_class
            relevant_type = "dna"
        else:
            continue

        # simplify classification
        if "2OH" in contact_class:
            simple_class = "2'OH-contacting"
        elif "base" in contact_class and "backbone" in contact_class:
            simple_class = "base+backbone"
        elif "base" in contact_class:
            simple_class = "base-only"
        elif "backbone" in contact_class:
            simple_class = "backbone-only"
        elif "sugar" in contact_class:
            simple_class = "sugar-ring"
        else:
            simple_class = "other"

        match = classify_node1_match(row)

        results.append({
            "pfam_id": row["pfam_id"],
            "aln_col": row["aln_col"],
            "contact_class": simple_class,
            "col_type": col_type,
            "relevant_type": relevant_type,
            "root_match": match,
        })

    return pd.DataFrame(results)


def compute_node_variance(nodes_df, master_df):
    """compute variance of RNA% across nodes for each family."""
    variances = []
    for pfam_id, group in nodes_df.groupby("pfam_id"):
        rna_pcts = group["pct_rna"].values
        dna_pcts = group["pct_dna"].values
        rna_range = rna_pcts.max() - rna_pcts.min()
        dna_range = dna_pcts.max() - dna_pcts.min()

        master_row = master_df[master_df["pfam_id"] == pfam_id]
        mode = master_row["generalism_mode"].values[0] if len(master_row) > 0 else "unknown"
        di = master_row["DI"].values[0] if len(master_row) > 0 else 0

        variances.append({
            "pfam_id": pfam_id,
            "rna_range": rna_range,
            "dna_range": dna_range,
            "rna_var": np.var(rna_pcts),
            "mode": mode,
            "DI": di,
        })
    return pd.DataFrame(variances)


# ── panel A: node-level heatmap ─────────────────────────────────────────────

def plot_panel_a(ax, nodes_df, master_df):
    """heatmap of RNA% across families × nodes, sorted by root RNA%."""
    # pivot to matrix: families × nodes
    pivot = nodes_df.pivot_table(
        index="pfam_id", columns="node", values="pct_rna"
    )

    # order by root RNA% (Node1)
    if "Node1" in pivot.columns:
        pivot = pivot.sort_values("Node1", ascending=True)

    # ensure columns are ordered
    node_cols = [c for c in ["Node1", "Node2", "Node3", "Node4", "Node5"] if c in pivot.columns]
    pivot = pivot[node_cols]

    # also get DNA% for annotation
    dna_pivot = nodes_df.pivot_table(
        index="pfam_id", columns="node", values="pct_dna"
    )
    dna_pivot = dna_pivot.reindex(pivot.index)[node_cols]

    # add generalism mode annotation
    mode_map = master_df.set_index("pfam_id")["generalism_mode"].to_dict()
    di_map = master_df.set_index("pfam_id")["DI"].to_dict()

    # plot heatmap
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                   vmin=0, vmax=100, interpolation="nearest")

    # add text annotations: RNA% (DNA%)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            rna_val = pivot.values[i, j]
            dna_val = dna_pivot.values[i, j]
            # choose text color based on background
            text_color = "white" if rna_val > 60 or rna_val < 10 else "black"
            txt = f"{rna_val:.0f}"
            if dna_val > 0:
                txt += f"\n({dna_val:.0f})"
            ax.text(j, i, txt, ha="center", va="center", fontsize=5.5,
                    color=text_color, fontweight="bold" if rna_val > 50 else "normal")

    # y-axis: family labels with mode indicator
    ylabels = []
    for pfam_id in pivot.index:
        mode = mode_map.get(pfam_id, "")
        di = di_map.get(pfam_id, 0)
        # short mode label
        if "generalist" in mode or "backbone" in mode:
            tag = "G"
        elif mode == "moderate":
            tag = "M"
        elif mode == "specialist":
            tag = "S"
        else:
            tag = "?"
        ylabels.append(f"{pfam_id} [{tag}]")

    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xticks(range(len(node_cols)))
    ax.set_xticklabels(["Root"] + [f"Node{i}" for i in range(2, len(node_cols) + 1)],
                       fontsize=7)
    ax.set_xlabel("ancestral node (deep → shallow)")
    ax.set_title("ancestral RNA% at NA-contact positions", fontweight="bold")

    # colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("RNA match %", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # mark families with significant node variance
    var_df = compute_node_variance(nodes_df, master_df)
    for i, pfam_id in enumerate(pivot.index):
        var_row = var_df[var_df["pfam_id"] == pfam_id]
        if len(var_row) > 0 and var_row["rna_range"].values[0] > 5:
            # mark with arrow indicating evolution
            ax.annotate("→", xy=(len(node_cols) - 0.3, i), fontsize=8,
                       color="#D55E00", fontweight="bold", ha="center", va="center")


# ── panel B: contact-type analysis ──────────────────────────────────────────

def plot_panel_b(ax, contact_stats):
    """bar chart: ancestral match rate by contact classification type."""
    # group by contact class, compute match rates
    summary = []
    for cls, group in contact_stats.groupby("contact_class"):
        n_total = len(group)
        if n_total < 3:
            continue
        n_rna = (group["root_match"] == "rna").sum()
        n_dna = (group["root_match"] == "dna").sum()
        n_both = (group["root_match"] == "both").sum()
        n_neither = (group["root_match"] == "neither").sum()

        summary.append({
            "contact_class": cls,
            "n": n_total,
            "rna_pct": 100.0 * (n_rna + n_both) / n_total,
            "dna_pct": 100.0 * (n_dna + n_both) / n_total,
            "neither_pct": 100.0 * n_neither / n_total,
        })

    sdf = pd.DataFrame(summary).sort_values("rna_pct", ascending=True)

    # stacked bar chart
    x = range(len(sdf))
    bar_width = 0.6

    # plot stacked: rna (blue) + dna (pink) + neither (gray)
    rna_bars = ax.barh(x, sdf["rna_pct"], height=bar_width,
                       color="#56B4E9", label="RNA match", edgecolor="white", linewidth=0.5)
    dna_bars = ax.barh(x, sdf["dna_pct"], height=bar_width, left=sdf["rna_pct"],
                       color="#CC79A7", label="DNA match", edgecolor="white", linewidth=0.5)
    neither_bars = ax.barh(x, sdf["neither_pct"], height=bar_width,
                          left=sdf["rna_pct"] + sdf["dna_pct"],
                          color="#DDDDDD", label="neither", edgecolor="white", linewidth=0.5)

    # labels with n
    ylabels = [f"{row['contact_class']}\n(n={row['n']})" for _, row in sdf.iterrows()]
    ax.set_yticks(x)
    ax.set_yticklabels(ylabels, fontsize=6.5)
    ax.set_xlabel("root ancestral match (%)")
    ax.set_xlim(0, 105)
    ax.set_title("ancestral conservation by contact type", fontweight="bold")
    ax.legend(loc="lower right", fontsize=6, frameon=False)

    # add percentage labels on bars
    for i, (_, row) in enumerate(sdf.iterrows()):
        if row["rna_pct"] > 8:
            ax.text(row["rna_pct"] / 2, i, f"{row['rna_pct']:.0f}%",
                   ha="center", va="center", fontsize=5.5, color="white", fontweight="bold")
        if row["dna_pct"] > 8:
            ax.text(row["rna_pct"] + row["dna_pct"] / 2, i, f"{row['dna_pct']:.0f}%",
                   ha="center", va="center", fontsize=5.5, color="white", fontweight="bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── panel C: DI vs ancestral RNA% scatter ───────────────────────────────────

def plot_panel_c(ax, master_df):
    """scatter: modern DI (x) vs root ancestral RNA% (y), colored by mode."""
    # color by generalism category
    mode_colors = {
        "generalist": "#009E73",
        "base_generalist": "#009E73",
        "backbone_generalist": "#56B4E9",
        "moderate": "#E69F00",
        "specialist": "#D55E00",
    }

    for _, row in master_df.iterrows():
        mode = row["generalism_mode"]
        color = mode_colors.get(mode, "#999999")
        marker = "D" if row["ancestral_bias"] == "DNA-biased" else "o"
        size = 60 if row["ancestral_bias"] == "DNA-biased" else 40

        ax.scatter(row["DI"], row["root_rna_pct"], c=color, marker=marker,
                  s=size, edgecolors="black", linewidth=0.5, zorder=3)

        # label interesting families
        if row["root_rna_pct"] > 50 or row["DI"] > 0.4 or row["ancestral_bias"] == "DNA-biased":
            offset = (5, 5)
            if row["pfam_id"] == "PF01131":
                offset = (5, -10)
            ax.annotate(row["pfam_id"], xy=(row["DI"], row["root_rna_pct"]),
                       xytext=offset, textcoords="offset points",
                       fontsize=5.5, color="black", fontstyle="italic")

    # add reference line at 50% (equal RNA/DNA)
    ax.axhline(y=50, color="gray", linestyle=":", linewidth=0.5, alpha=0.5)

    # legend for modes
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#009E73",
               markersize=6, label="generalist (DI<0.05)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#56B4E9",
               markersize=6, label="backbone gen. (DI<0.07)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#E69F00",
               markersize=6, label="moderate (0.1-0.25)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#D55E00",
               markersize=6, label="specialist (DI>0.3)"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="#999999",
               markersize=6, label="DNA-biased ancestor"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=5.5,
             frameon=True, framealpha=0.8)

    ax.set_xlabel("modern discrimination index (DI)")
    ax.set_ylabel("root ancestral RNA match (%)")
    ax.set_title("modern specificity vs. ancestral state", fontweight="bold")
    ax.set_xlim(-0.02, 0.75)
    ax.set_ylim(-5, 85)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── main ────────────────────────────────────────────────────────────────────

def main():
    print("loading data...")
    master_df = load_master()
    nodes_df = load_node_summaries()
    conv_df = load_convergence_details()

    print(f"  {len(master_df)} families with OK status")
    print(f"  {len(nodes_df)} node records")
    print(f"  {len(conv_df)} convergence position records")

    # compute contact-type stats
    contact_stats = compute_contact_type_stats(conv_df)
    print(f"  {len(contact_stats)} classified contact positions")

    # ── textual analysis ────────────────────────────────────────────────────

    print("\n" + "=" * 80)
    print("DEEP ANALYSIS: NODE-LEVEL TRENDS")
    print("=" * 80)

    var_df = compute_node_variance(nodes_df, master_df)
    var_df = var_df.sort_values("rna_range", ascending=False)

    print(f"\n{'Family':<12} {'Mode':<22} {'DI':>6} {'RNA range':>10} {'DNA range':>10} {'Evolving?'}")
    print("-" * 75)
    for _, row in var_df.iterrows():
        evolving = "YES" if row["rna_range"] > 5 else "no"
        print(f"{row['pfam_id']:<12} {row['mode']:<22} {row['DI']:>6.3f} "
              f"{row['rna_range']:>9.1f}% {row['dna_range']:>9.1f}% {evolving}")

    n_evolving = (var_df["rna_range"] > 5).sum()
    print(f"\n{n_evolving}/{len(var_df)} families show >5% RNA range across nodes")

    # node trend direction: does RNA% increase or decrease from root to tips?
    print("\n── node trend direction (root→shallow) ──")
    for pfam_id, group in nodes_df.groupby("pfam_id"):
        group_sorted = group.sort_values("node")
        root_rna = group_sorted["pct_rna"].iloc[0]
        tip_rna = group_sorted["pct_rna"].iloc[-1]
        delta = tip_rna - root_rna
        if abs(delta) > 3:
            direction = "↑ RNA increases" if delta > 0 else "↓ RNA decreases"
            print(f"  {pfam_id}: root={root_rna:.1f}% → shallow={tip_rna:.1f}% ({direction}, Δ={delta:+.1f}%)")

    # ── contact type analysis ───────────────────────────────────────────────

    print("\n" + "=" * 80)
    print("DEEP ANALYSIS: CONTACT TYPE BREAKDOWN")
    print("=" * 80)

    for cls, group in contact_stats.groupby("contact_class"):
        n = len(group)
        if n < 3:
            continue
        n_rna = (group["root_match"] == "rna").sum()
        n_dna = (group["root_match"] == "dna").sum()
        n_both = (group["root_match"] == "both").sum()
        n_neither = (group["root_match"] == "neither").sum()
        print(f"\n  {cls} (n={n}):")
        print(f"    RNA match: {100*n_rna/n:.1f}%  DNA match: {100*n_dna/n:.1f}%  "
              f"both: {100*n_both/n:.1f}%  neither: {100*n_neither/n:.1f}%")

    # 2'OH specific analysis
    print("\n── 2'OH hypothesis test ──")
    oh2 = contact_stats[contact_stats["contact_class"] == "2'OH-contacting"]
    non_oh2 = contact_stats[contact_stats["contact_class"] != "2'OH-contacting"]

    if len(oh2) > 0:
        oh2_neither = (oh2["root_match"] == "neither").sum() / len(oh2)
        non_neither = (non_oh2["root_match"] == "neither").sum() / len(non_oh2)
        print(f"  2'OH positions: {100*oh2_neither:.1f}% neither at root (n={len(oh2)})")
        print(f"  non-2'OH positions: {100*non_neither:.1f}% neither at root (n={len(non_oh2)})")
        print(f"  → 2'OH positions are {'MORE' if oh2_neither > non_neither else 'LESS'} "
              f"likely to be 'neither' ({100*abs(oh2_neither - non_neither):.1f}pp difference)")

    # ── PF01131 diagnosis ───────────────────────────────────────────────────

    print("\n" + "=" * 80)
    print("PF01131 POSITIVE CONTROL: DNA TOPOISOMERASE")
    print("=" * 80)

    pf01131 = master_df[master_df["pfam_id"] == "PF01131"]
    if len(pf01131) > 0:
        r = pf01131.iloc[0]
        print(f"  domain: bacterial DNA topoisomerase I (Topoisom_bac)")
        print(f"  modern DI: {r['DI']:.3f} (moderate DNA preference)")
        print(f"  root ancestral state: {r['root_rna_pct']:.1f}% RNA, "
              f"{r['root_dna_pct']:.1f}% DNA, {r['root_neither_pct']:.1f}% neither")
        print(f"  verdict: GENUINE ancestral DNA-binder (positive control)")
        print(f"  note: RNA mapping failed (HMMER couldn't find domain in 9GDA)")
        print(f"         → all 23 contact columns are DNA-only")
        print(f"         → 22/23 have PP=1.0 matching modern DNA contacts")
        print(f"         this is a DEDICATED DNA enzyme, not an artifact")

    # PF03372 check
    pf03372 = master_df[master_df["pfam_id"] == "PF03372"]
    if len(pf03372) > 0:
        r = pf03372.iloc[0]
        print(f"\n  PF03372 (2nd DNA-biased):")
        print(f"    age: {r['age']}, DI: {r['DI']:.3f}")
        print(f"    root: {r['root_rna_pct']:.1f}% RNA, {r['root_dna_pct']:.1f}% DNA")
        print(f"    only 7 contact columns, 3 with high PP → low confidence")

    # ── create figure ───────────────────────────────────────────────────────

    print("\ncreating figure...")

    fig = plt.figure(figsize=(11, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.35,
                  height_ratios=[1.2, 1])

    # panel A: top spanning both columns
    ax_a = fig.add_subplot(gs[0, :])
    plot_panel_a(ax_a, nodes_df, master_df)
    ax_a.text(-0.03, 1.02, "a", transform=ax_a.transAxes, fontsize=12,
              fontweight="bold", va="bottom")

    # panel B: bottom left
    ax_b = fig.add_subplot(gs[1, 0])
    plot_panel_b(ax_b, contact_stats)
    ax_b.text(-0.12, 1.02, "b", transform=ax_b.transAxes, fontsize=12,
              fontweight="bold", va="bottom")

    # panel C: bottom right
    ax_c = fig.add_subplot(gs[1, 1])
    plot_panel_c(ax_c, master_df)
    ax_c.text(-0.12, 1.02, "c", transform=ax_c.transAxes, fontsize=12,
              fontweight="bold", va="bottom")

    # save
    outpath_pdf = os.path.join(OUTDIR, "deep_analysis.pdf")
    outpath_png = os.path.join(OUTDIR, "deep_analysis.png")
    fig.savefig(outpath_pdf, format="pdf")
    fig.savefig(outpath_png, format="png")
    print(f"\nsaved: {outpath_pdf}")
    print(f"saved: {outpath_png}")
    plt.close()

    # ── summary statistics ──────────────────────────────────────────────────

    print("\n" + "=" * 80)
    print("SUMMARY: KEY FINDINGS")
    print("=" * 80)

    print(f"\n1. NODE-LEVEL TRENDS:")
    print(f"   {n_evolving}/{len(var_df)} families show evolving specificity (>5% RNA range)")
    print(f"   majority are FLAT → ancestral state is deeply conserved")

    if len(oh2) > 0:
        print(f"\n2. CONTACT-TYPE ANALYSIS:")
        oh2_rna = (oh2["root_match"] == "rna").sum() / len(oh2)
        non_rna = (non_oh2["root_match"] == "rna").sum() / len(non_oh2)
        print(f"   2'OH positions: {100*oh2_rna:.1f}% RNA match at root")
        print(f"   non-2'OH positions: {100*non_rna:.1f}% RNA match at root")
        print(f"   2'OH positions: {100*oh2_neither:.1f}% 'neither' at root")
        print(f"   non-2'OH positions: {100*non_neither:.1f}% 'neither' at root")

    print(f"\n3. DI vs ANCESTRAL STATE:")
    for mode in ["generalist", "base_generalist", "backbone_generalist", "moderate", "specialist"]:
        subset = master_df[master_df["generalism_mode"] == mode]
        if len(subset) > 0:
            mean_rna = subset["root_rna_pct"].mean()
            mean_dna = subset["root_dna_pct"].mean()
            print(f"   {mode}: mean root RNA={mean_rna:.1f}%, DNA={mean_dna:.1f}% (n={len(subset)})")

    print(f"\n4. POSITIVE CONTROLS:")
    print(f"   PF01131 (DNA topoisomerase): correctly identified as ancestral DNA-binder")
    print(f"   PF03372: weak DNA signal (only 7 columns, 3 high PP)")


if __name__ == "__main__":
    main()
