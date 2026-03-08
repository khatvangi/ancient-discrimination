#!/usr/bin/env python3
"""
publication figures for the ancient discrimination paper.

figure 1: NAS-Bench validation (modern specialists + LUCA comparison)
figure 2: ancestral contact chemistry (main result)
figure 3: lysine claw mechanism (PF03129 case study)
figure 4: cross-tool disagreement heatmap
figure 5: KH phylogeny (simplified: P_RNA vs P_DNA scatter)
"""

import os
import sys
import glob
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats

# ── global style ────────────────────────────────────────────────────────────

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

# color palette
RNA_COLOR = "#E74C3C"
DNA_COLOR = "#3498DB"
DUAL_COLOR = "#9B59B6"
GENERALIST_COLOR = "#3498DB"
MODERATE_COLOR = "#F39C12"
SPECIALIST_COLOR = "#E74C3C"
IDENTICAL_COLOR = "#27AE60"
CONSERVATIVE_COLOR = "#F1C40F"
RADICAL_COLOR = "#E74C3C"

BASE = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination"
RESULTS = os.path.join(BASE, "results")
ASR = os.path.join(RESULTS, "asr")
FIGDIR = os.path.join(BASE, "figures")
os.makedirs(FIGDIR, exist_ok=True)


def save_fig(fig, name):
    """save figure as png and pdf."""
    fig.savefig(os.path.join(FIGDIR, f"{name}.png"), format="png")
    fig.savefig(os.path.join(FIGDIR, f"{name}.pdf"), format="pdf")
    print(f"  saved: figures/{name}.png + .pdf")
    plt.close(fig)


# ============================================================================
# FIGURE 2 — ANCESTRAL CONTACT CHEMISTRY (priority 1)
# ============================================================================

def make_figure_2():
    """main result: ancestral residues are radically different from modern specialists."""
    print("\n=== FIGURE 2: Ancestral Contact Chemistry ===")

    # load all substitution analysis data
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
    all_df["di_cat"] = pd.cut(all_df["DI"], bins=[0, 0.10, 0.25, 1.0],
                               labels=["Generalist\n(DI<0.10)",
                                       "Moderate\n(0.10-0.25)",
                                       "Specialist\n(DI>0.25)"])

    fig, axes = plt.subplots(1, 3, figsize=(7, 3))

    # ── panel A: stacked bar by contact type ────────────────────────────────

    ax = axes[0]
    contact_order = ["backbone", "base", "2oh", "base+bb", "ALL"]
    contact_labels = ["Backbone", "Base", "2'OH", "Base+BB", "ALL"]

    pct_data = {"Identical": [], "Conservative": [], "Radical": []}
    n_data = []

    for ct in contact_order:
        if ct == "ALL":
            subset = all_df
        else:
            subset = all_df[all_df["contact_type"] == ct]
        n = len(subset)
        n_data.append(n)
        for cls in ["identical", "conservative", "radical"]:
            key = cls.capitalize()
            pct_data[key].append(100 * (subset["substitution_class"] == cls).sum() / n if n > 0 else 0)

    y_pos = np.arange(len(contact_order))
    left = np.zeros(len(contact_order))

    for cls, color in [("Identical", IDENTICAL_COLOR),
                        ("Conservative", CONSERVATIVE_COLOR),
                        ("Radical", RADICAL_COLOR)]:
        vals = pct_data[cls]
        bars = ax.barh(y_pos, vals, left=left, height=0.6, color=color,
                       edgecolor="white", linewidth=0.3, label=cls)
        # add percentage labels
        for i, (v, l) in enumerate(zip(vals, left)):
            if v > 10:
                ax.text(l + v / 2, i, f"{v:.0f}%", ha="center", va="center",
                       fontsize=5.5, fontweight="bold",
                       color="white" if cls != "Conservative" else "black")
        left += vals

    # add n= labels on right
    for i, n in enumerate(n_data):
        ax.text(102, i, f"n={n}", ha="left", va="center", fontsize=5.5, color="gray")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(contact_labels, fontsize=7)
    ax.set_xlabel("Root ancestral match (%)")
    ax.set_xlim(0, 115)
    ax.set_title("a", fontweight="bold", loc="left", fontsize=12)
    ax.legend(loc="lower right", fontsize=6, frameon=False)

    # ── panel B: grouped bars by DI category ────────────────────────────────

    ax = axes[1]
    cats = ["Generalist\n(DI<0.10)", "Moderate\n(0.10-0.25)", "Specialist\n(DI>0.25)"]
    cat_colors = [GENERALIST_COLOR, MODERATE_COLOR, SPECIALIST_COLOR]

    all_radical = []
    oh2_radical = []

    for cat in cats:
        cat_df = all_df[all_df["di_cat"] == cat]
        n_total = len(cat_df)
        n_rad = (cat_df["substitution_class"] == "radical").sum()
        all_radical.append(100 * n_rad / n_total if n_total > 0 else 0)

        oh2_df = cat_df[cat_df["contact_type"] == "2oh"]
        n_oh2 = len(oh2_df)
        n_oh2_rad = (oh2_df["substitution_class"] == "radical").sum()
        oh2_radical.append(100 * n_oh2_rad / n_oh2 if n_oh2 > 0 else 0)

    x = np.arange(len(cats))
    w = 0.35
    bars1 = ax.bar(x - w / 2, all_radical, w, color=[c + "80" for c in cat_colors],
                   edgecolor=cat_colors, linewidth=1, label="All contacts")
    bars2 = ax.bar(x + w / 2, oh2_radical, w, color=cat_colors,
                   edgecolor="black", linewidth=0.5, label="2'OH only")

    # value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 1,
                   f"{h:.0f}%", ha="center", va="bottom", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=6.5)
    ax.set_ylabel("% Radical substitutions")
    ax.set_ylim(0, 85)
    ax.set_title("b", fontweight="bold", loc="left", fontsize=12)
    ax.legend(loc="upper right", fontsize=6, frameon=False)

    # ── panel C: RNA vs DNA direction ───────────────────────────────────────

    ax = axes[2]
    rna_df = all_df[all_df["target_type"] == "RNA"]
    dna_df = all_df[all_df["target_type"] == "DNA"]

    rna_pcts = []
    dna_pcts = []
    for cls in ["identical", "conservative", "radical"]:
        rna_pcts.append(100 * (rna_df["substitution_class"] == cls).sum() / len(rna_df))
        dna_pcts.append(100 * (dna_df["substitution_class"] == cls).sum() / len(dna_df))

    x = np.arange(3)
    w = 0.35
    labels_x = ["Identical", "Conservative", "Radical"]
    colors = [IDENTICAL_COLOR, CONSERVATIVE_COLOR, RADICAL_COLOR]

    for i, (lbl, col) in enumerate(zip(labels_x, colors)):
        ax.bar(i - w / 2, rna_pcts[i], w, color=col, alpha=0.5,
               edgecolor=col, linewidth=1)
        ax.bar(i + w / 2, dna_pcts[i], w, color=col,
               edgecolor="black", linewidth=0.5)
        # labels
        ax.text(i - w / 2, rna_pcts[i] + 1, f"{rna_pcts[i]:.0f}",
               ha="center", fontsize=5.5, color=RNA_COLOR)
        ax.text(i + w / 2, dna_pcts[i] + 1, f"{dna_pcts[i]:.0f}",
               ha="center", fontsize=5.5, color=DNA_COLOR)

    ax.set_xticks(x)
    ax.set_xticklabels(labels_x, fontsize=7)
    ax.set_ylabel("% of contact positions")
    ax.set_ylim(0, 75)
    ax.set_title("c", fontweight="bold", loc="left", fontsize=12)

    # custom legend for RNA/DNA
    legend_elements = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="gray",
               alpha=0.5, markersize=8, label=f"vs RNA specialist (n={len(rna_df)})"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="gray",
               markersize=8, label=f"vs DNA specialist (n={len(dna_df)})"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=5.5, frameon=False)

    fig.tight_layout()
    save_fig(fig, "fig2_ancestral_chemistry")


# ============================================================================
# FIGURE 3 — LYSINE CLAW MECHANISM (priority 2)
# ============================================================================

def make_figure_3():
    """PF03129 case study: ancestral lysine cluster -> modern diverse contacts."""
    print("\n=== FIGURE 3: Lysine Claw Mechanism ===")

    # load PF03129 convergence data
    conv = pd.read_csv(os.path.join(ASR, "pf03129/PF03129_convergence_summary.tsv"), sep="\t")

    # key columns for the lysine cluster
    key_cols = [67, 68, 70, 71]
    # context columns
    all_cols = [66, 67, 68, 69, 70, 71, 72, 73, 74]

    # physicochemical color map
    aa_colors = {
        "K": "#3498DB", "R": "#3498DB",
        "D": "#E74C3C", "E": "#E74C3C",
        "F": "#9B59B6", "Y": "#9B59B6", "W": "#9B59B6", "H": "#9B59B6",
        "S": "#27AE60", "T": "#27AE60", "N": "#27AE60", "Q": "#27AE60",
        "A": "#BDC3C7", "V": "#BDC3C7", "I": "#BDC3C7", "L": "#BDC3C7", "M": "#BDC3C7",
        "G": "#F39C12", "P": "#F39C12", "C": "#F39C12",
    }

    fig = plt.figure(figsize=(7, 4.5))
    gs = GridSpec(1, 2, figure=fig, wspace=0.4, width_ratios=[1.2, 1])

    # ── panel A: alignment view ─────────────────────────────────────────────

    ax = fig.add_subplot(gs[0])

    col_data = []
    for col in all_cols:
        row = conv[conv["aln_col"] == col]
        if len(row) == 0:
            continue
        r = row.iloc[0]
        rna_aa = str(r.get("rna_aa", "")).strip()
        dna_aa = str(r.get("dna_aa", "")).strip()
        if rna_aa == "nan": rna_aa = "-"
        if dna_aa == "nan": dna_aa = "-"

        rna_cls = str(r.get("rna_classification", ""))
        dna_cls = str(r.get("dna_classification", ""))
        if "2OH" in rna_cls.upper():
            ct = "2'OH"
        elif "base" in rna_cls.lower() and "backbone" in rna_cls.lower():
            ct = "base+bb"
        elif "base" in rna_cls.lower():
            ct = "base"
        elif "backbone" in rna_cls.lower():
            ct = "bb"
        elif "backbone" in dna_cls.lower():
            ct = "bb"
        elif "base" in dna_cls.lower():
            ct = "base"
        else:
            ct = "?"

        col_data.append({
            "col": col, "rna_aa": rna_aa, "dna_aa": dna_aa,
            "root": str(r["Node1_state"]),
            "root_pp": float(r["Node1_pp"]),
            "node2": str(r["Node2_state"]),
            "node3": str(r["Node3_state"]),
            "contact": ct,
            "col_type": str(r["col_type"]),
            "is_key": col in key_cols,
        })

    n_cols = len(col_data)
    row_labels = ["Modern RNA", "Modern DNA", "", "Node 3", "Node 2", "Root (LUCA)"]
    n_rows = len(row_labels)

    for j, cd in enumerate(col_data):
        residues = [cd["rna_aa"], cd["dna_aa"], "",
                    cd["node3"], cd["node2"], cd["root"]]

        for i, res in enumerate(residues):
            if res == "" or res == "-" or res == "nan":
                continue

            color = aa_colors.get(res, "#EEEEEE")
            alpha = 1.0 if cd["is_key"] else 0.6
            edgewidth = 1.5 if cd["is_key"] else 0.5
            edgecol = "black" if cd["is_key"] else "#888888"

            rect = plt.Rectangle((j - 0.4, n_rows - 1 - i - 0.4), 0.8, 0.8,
                                facecolor=color, edgecolor=edgecol,
                                linewidth=edgewidth, alpha=alpha, zorder=2)
            ax.add_patch(rect)
            ax.text(j, n_rows - 1 - i, res, ha="center", va="center",
                   fontsize=9, fontweight="bold", color="white", zorder=3)

    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(-2.5, n_rows - 0.3)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels([str(cd["col"]) for cd in col_data], fontsize=6.5)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(list(reversed(row_labels)), fontsize=7)
    ax.set_xlabel("Alignment column")

    for j, cd in enumerate(col_data):
        ax.text(j, -1.2, cd["contact"], ha="center", va="center",
               fontsize=5.5, fontstyle="italic", color="#666666")
        ax.text(j, -2.0, f"{cd['root_pp']:.2f}", ha="center", va="center",
               fontsize=5.5, color="#666666")

    ax.text(-0.6, -1.2, "Contact:", ha="right", va="center", fontsize=5.5,
           fontstyle="italic", color="#666666")
    ax.text(-0.6, -2.0, "Root PP:", ha="right", va="center", fontsize=5.5,
           color="#666666")

    ax.axhline(y=n_rows - 3.5, color="black", linewidth=1, linestyle="--", alpha=0.5)

    # physicochemical legend
    legend_items = [
        ("K/R", "#3498DB", "Basic (electrostatic)"),
        ("D/E", "#E74C3C", "Acidic"),
        ("F/Y/W/H", "#9B59B6", "Aromatic (stacking)"),
        ("S/T/N/Q", "#27AE60", "Polar"),
        ("A-M", "#BDC3C7", "Hydrophobic"),
        ("G/P/C", "#F39C12", "Special"),
    ]
    for k, (aas, color, label) in enumerate(legend_items):
        y = n_rows - 1 - k * 0.55
        ax.plot(n_cols + 0.3, y, "s", color=color, markersize=5)
        ax.text(n_cols + 0.7, y, f"{aas}: {label}", fontsize=5, va="center")

    ax.set_xlim(-0.8, n_cols + 4.0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False)
    ax.set_title("a", fontweight="bold", loc="left", fontsize=12)

    # ── panel B: schematic diagram ──────────────────────────────────────────

    ax = fig.add_subplot(gs[1])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("b", fontweight="bold", loc="left", fontsize=12)

    # ancestral mode (left side)
    ax.text(2.5, 9.5, "Ancestral mode", ha="center", fontsize=9, fontweight="bold")
    ax.text(2.5, 9.0, "(LUCA)", ha="center", fontsize=7, fontstyle="italic", color="gray")

    # phosphate backbone
    ax.plot([1.0, 1.0], [2, 8], color="#888888", linewidth=3, zorder=1)
    for y in [3, 4.5, 6, 7.5]:
        ax.plot(1.0, y, "o", color="#E74C3C", markersize=8, zorder=2)
        ax.text(0.3, y, "P", ha="center", va="center", fontsize=6, color="#E74C3C")

    # ancestral Lys residues
    for y, col_num in zip([3, 4.5, 6, 7.5], [67, 68, 70, 71]):
        ax.annotate("", xy=(1.5, y), xytext=(3.5, y),
                    arrowprops=dict(arrowstyle="->", color="#3498DB", lw=1.5))
        ax.plot(3.8, y, "s", color="#3498DB", markersize=12, zorder=2)
        ax.text(3.8, y, "K", ha="center", va="center", fontsize=8,
               fontweight="bold", color="white", zorder=3)
        ax.text(4.3, y + 0.25, "+", fontsize=8, color="#3498DB", fontweight="bold")

    ax.text(2.5, 1.5, "Electrostatic grip\n(nonspecific)", ha="center",
           fontsize=7, fontstyle="italic", color="#3498DB")

    # modern mode (right side)
    ax.text(7.5, 9.5, "Modern mode", ha="center", fontsize=9, fontweight="bold")
    ax.text(7.5, 9.0, "(tRNA synthetase)", ha="center", fontsize=7,
           fontstyle="italic", color="gray")

    ax.plot([5.5, 5.5], [2, 8], color="#888888", linewidth=3, zorder=1)
    base_labels = ["A", "U", "G", "C"]
    for y, base in zip([3, 4.5, 6, 7.5], base_labels):
        ax.plot(5.5, y, "o", color="#E74C3C", markersize=6, zorder=2)
        ax.plot(6.3, y, "D", color="#F39C12", markersize=10, zorder=2)
        ax.text(6.3, y, base, ha="center", va="center", fontsize=6,
               fontweight="bold", color="black", zorder=3)

    modern_residues = [
        (3, "E", "#E74C3C", "base", "charge\nreversal"),
        (4.5, "P", "#F39C12", "2oh", "2'OH\nrecognition"),
        (6, "G", "#F39C12", "space", "space\nfor fit"),
        (7.5, "Y", "#9B59B6", "stack", "aromatic\nstacking"),
    ]

    for y, res, color, contact_type, label in modern_residues:
        if contact_type == "stack":
            ax.annotate("", xy=(6.6, y), xytext=(8.0, y),
                        arrowprops=dict(arrowstyle="<->", color=color, lw=1.5,
                                       linestyle="--"))
        elif contact_type == "2oh":
            ax.annotate("", xy=(5.9, y - 0.2), xytext=(8.0, y),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.5))
        elif contact_type == "base":
            ax.annotate("", xy=(6.6, y), xytext=(8.0, y),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.5))
        else:
            ax.annotate("", xy=(6.3, y), xytext=(8.0, y),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.5,
                                       linestyle=":"))

        ax.plot(8.3, y, "s", color=color, markersize=12, zorder=2)
        ax.text(8.3, y, res, ha="center", va="center", fontsize=8,
               fontweight="bold", color="white", zorder=3)
        ax.text(9.2, y, label, ha="left", va="center", fontsize=5.5,
               color=color, fontstyle="italic")

    ax.text(7.5, 1.5, "Diverse specific contacts\n(substrate-discriminating)",
           ha="center", fontsize=7, fontstyle="italic", color="#9B59B6")

    # arrow between
    ax.annotate("", xy=(5.0, 5.25), xytext=(4.5, 5.25),
                arrowprops=dict(arrowstyle="->", color="black", lw=2))
    ax.text(4.75, 5.8, "Specificity\nevolution", ha="center", fontsize=6.5,
           fontweight="bold")

    fig.tight_layout()
    save_fig(fig, "fig3_lysine_claw")


# ============================================================================
# FIGURE 1 — NAS-BENCH VALIDATION (priority 3)
# ============================================================================

def make_figure_1():
    """tool validation: modern specialists cleanly separated, LUCA domains shifted."""
    print("\n=== FIGURE 1: NAS-Bench Validation ===")

    modern = pd.read_csv(os.path.join(RESULTS, "nasbench_modern_controls.tsv"), sep="\t")
    luca = pd.read_csv(os.path.join(RESULTS, "nasbench_full_luca.tsv"), sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(7, 3))

    # ── panel A: modern SI comparison ───────────────────────────────────────

    ax = axes[0]
    rna_specs = modern[modern["category"] == "RNA_specialist"]
    dna_specs = modern[modern["category"] == "DNA_specialist"]

    jitter = 0.15
    for i, (group, color, label) in enumerate([
        (rna_specs, RNA_COLOR, "RNA specialists"),
        (dna_specs, DNA_COLOR, "DNA specialists"),
    ]):
        x_base = i
        np.random.seed(42)
        x_jitter = x_base + np.random.uniform(-jitter, jitter, len(group))
        ax.scatter(x_jitter, group["SI"], c=color, s=30, alpha=0.8,
                  edgecolors="black", linewidth=0.3, zorder=3)

        for j, (_, row) in enumerate(group.iterrows()):
            name = row["name"]
            if len(name) > 12:
                name = name[:10] + ".."
            ax.text(x_jitter[j] + 0.18, group["SI"].iloc[j], name,
                   fontsize=4, va="center", color=color, alpha=0.8)

    u_stat, p_val = stats.mannwhitneyu(rna_specs["SI"], dna_specs["SI"],
                                        alternative="two-sided")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["RNA\nspecialists", "DNA\nspecialists"], fontsize=7)
    ax.set_ylabel("Specificity Index (SI)")
    ax.set_ylim(-0.05, 1.1)
    ax.set_title("a", fontweight="bold", loc="left", fontsize=12)

    y_bracket = 1.03
    ax.plot([0, 0, 1, 1], [y_bracket - 0.02, y_bracket, y_bracket, y_bracket - 0.02],
            color="black", linewidth=0.8)
    ax.text(0.5, y_bracket + 0.01, f"p = {p_val:.4f}", ha="center", fontsize=6)

    # ── panel B: LUCA SI_RNA vs SI_DNA scatter ──────────────────────────────

    ax = axes[1]

    di_colors = []
    for _, row in luca.iterrows():
        di = row["DI"]
        if di < 0.10:
            di_colors.append(GENERALIST_COLOR)
        elif di < 0.25:
            di_colors.append(MODERATE_COLOR)
        else:
            di_colors.append(SPECIALIST_COLOR)

    ax.scatter(luca["dna_SI"], luca["rna_SI"], c=di_colors, s=40,
              edgecolors="black", linewidth=0.3, zorder=3)

    ax.plot([0, 1], [0, 1], "--", color="gray", linewidth=0.5, alpha=0.5, zorder=1)

    for _, row in luca.iterrows():
        if row["DI"] > 0.4 or row["DI"] < 0.02 or row["pfam_id"] == "PF01131":
            ax.text(row["dna_SI"] + 0.02, row["rna_SI"] + 0.02,
                   row["pfam_id"], fontsize=4.5, fontstyle="italic")

    ax.set_xlabel("SI on DNA structure")
    ax.set_ylabel("SI on RNA structure")
    ax.set_xlim(-0.05, 0.75)
    ax.set_ylim(-0.05, 0.85)
    ax.set_title("b", fontweight="bold", loc="left", fontsize=12)
    ax.set_aspect("equal")

    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=GENERALIST_COLOR,
               markersize=6, label="Generalist (DI<0.10)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=MODERATE_COLOR,
               markersize=6, label="Moderate (0.10-0.25)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=SPECIALIST_COLOR,
               markersize=6, label="Specialist (DI>0.25)"),
    ]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=5, frameon=False)

    # ── panel C: DI distribution ────────────────────────────────────────────

    ax = axes[2]

    luca_di = luca["DI"].dropna().values

    parts = ax.violinplot([luca_di], positions=[0], showmeans=True,
                          showmedians=True, widths=0.6)
    for pc in parts["bodies"]:
        pc.set_facecolor(GENERALIST_COLOR)
        pc.set_alpha(0.4)
    parts["cmeans"].set_color("black")
    parts["cmedians"].set_color(SPECIALIST_COLOR)

    np.random.seed(42)
    x_jitter = np.random.uniform(-0.15, 0.15, len(luca_di))
    ax.scatter(x_jitter, luca_di, c=GENERALIST_COLOR, s=15, alpha=0.6,
              edgecolors="black", linewidth=0.3, zorder=3)

    ax.axhline(y=0.10, color="gray", linestyle=":", linewidth=0.5, alpha=0.7)
    ax.axhline(y=0.25, color="gray", linestyle=":", linewidth=0.5, alpha=0.7)
    ax.text(0.4, 0.10, "DI=0.10", fontsize=5, va="bottom", color="gray")
    ax.text(0.4, 0.25, "DI=0.25", fontsize=5, va="bottom", color="gray")

    median_di = np.median(luca_di)
    ax.text(-0.5, median_di, f"median={median_di:.2f}", fontsize=5.5,
           ha="right", va="center", color=SPECIALIST_COLOR)

    ax.set_xticks([0])
    ax.set_xticklabels(["LUCA domains\n(n=18)"], fontsize=7)
    ax.set_ylabel("Discrimination Index (DI)")
    ax.set_ylim(-0.05, 0.75)
    ax.set_title("c", fontweight="bold", loc="left", fontsize=12)

    fig.tight_layout()
    save_fig(fig, "fig1_validation")


# ============================================================================
# FIGURE 4 — CROSS-TOOL DISAGREEMENT HEATMAP (priority 4)
# ============================================================================

def make_figure_4():
    """heatmap: three tools disagree on LUCA classifications."""
    print("\n=== FIGURE 4: Cross-Tool Disagreement ===")

    nasbench = pd.read_csv(os.path.join(RESULTS, "nasbench_full_luca.tsv"), sep="\t")
    prona = pd.read_csv(os.path.join(RESULTS, "prona2020_luca_domains.tsv"), sep="\t")
    drbp = pd.read_csv(os.path.join(RESULTS, "drbp_edp_luca_domains.tsv"), sep="\t")

    families = nasbench.sort_values("DI")["pfam_id"].tolist()

    # classify each tool
    nb_class = {}
    for _, row in nasbench.iterrows():
        di = row["DI"]
        delta = row["delta_SI"]
        if di < 0.10:
            nb_class[row["pfam_id"]] = "dual"
        elif delta > 0:
            nb_class[row["pfam_id"]] = "RNA"
        else:
            nb_class[row["pfam_id"]] = "DNA"

    prona_class = {}
    for _, row in prona.iterrows():
        if row["pfam_id"] not in families:
            continue
        cls = str(row.get("predicted_class", "")).strip()
        if cls == "dual":
            prona_class[row["pfam_id"]] = "dual"
        elif "RNA" in cls:
            prona_class[row["pfam_id"]] = "RNA"
        elif "DNA" in cls:
            prona_class[row["pfam_id"]] = "DNA"
        else:
            prona_class[row["pfam_id"]] = "other"

    drbp_class = {}
    for _, row in drbp.iterrows():
        if row["pfam_id"] not in families:
            continue
        s2 = str(row.get("stage2_pred", "")).strip()
        if s2 == "RNA":
            drbp_class[row["pfam_id"]] = "RNA"
        elif s2 == "DNA":
            drbp_class[row["pfam_id"]] = "DNA"
        else:
            s1 = str(row.get("stage1_pred", "")).strip()
            if "NAB" in s1:
                drbp_class[row["pfam_id"]] = "dual"
            else:
                drbp_class[row["pfam_id"]] = "other"

    class_to_num = {"RNA": 0, "DNA": 1, "dual": 2, "other": 3}
    class_colors = [RNA_COLOR, DNA_COLOR, DUAL_COLOR, "#CCCCCC"]

    matrix = np.full((len(families), 3), 3)
    tools = ["NAS-Bench", "ProNA2020", "DRBP-EDP"]

    for i, fam in enumerate(families):
        matrix[i, 0] = class_to_num.get(nb_class.get(fam, "other"), 3)
        matrix[i, 1] = class_to_num.get(prona_class.get(fam, "other"), 3)
        matrix[i, 2] = class_to_num.get(drbp_class.get(fam, "other"), 3)

    fig, ax = plt.subplots(figsize=(3.5, 5))

    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(class_colors)

    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5,
                   interpolation="nearest")

    for x in [-0.5, 0.5, 1.5, 2.5]:
        ax.axvline(x, color="white", linewidth=1)
    for y in np.arange(-0.5, len(families), 1):
        ax.axhline(y, color="white", linewidth=0.5)

    ax.set_xticks(range(3))
    ax.set_xticklabels(tools, fontsize=7, rotation=30, ha="right")
    ax.set_yticks(range(len(families)))

    di_map = nasbench.set_index("pfam_id")["DI"].to_dict()
    ylabels = [f"{fam} ({di_map.get(fam, 0):.2f})" for fam in families]
    ax.set_yticklabels(ylabels, fontsize=5.5)

    ax.set_title("Cross-tool classification\nof LUCA families", fontsize=9, fontweight="bold")

    legend_elements = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor=RNA_COLOR,
               markersize=8, label="RNA"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=DNA_COLOR,
               markersize=8, label="DNA"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=DUAL_COLOR,
               markersize=8, label="Dual/ambiguous"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#CCCCCC",
               markersize=8, label="Not classified"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=6,
             frameon=True, framealpha=0.9)

    n_agree = sum(1 for i in range(len(families)) if len(set(matrix[i, :])) == 1)
    ax.text(0.5, -0.08, f"Full agreement: {n_agree}/{len(families)} families",
           transform=ax.transAxes, ha="center", fontsize=7, fontstyle="italic")

    fig.tight_layout()
    save_fig(fig, "fig4_cross_tool")


# ============================================================================
# FIGURE 5 — KH DOMAIN PRONA2020 (priority 5)
# ============================================================================

def make_figure_5():
    """KH domain: ProNA2020 P_RNA vs P_DNA for all nodes."""
    print("\n=== FIGURE 5: KH Domain ProNA2020 ===")

    kh_prona = pd.read_csv(os.path.join(ASR, "kh/kh_prona_all_nodes.tsv"), sep="\t")

    fig, ax = plt.subplots(figsize=(5, 4))

    # column names: node_id, node_type, organism, domain_of_life, P_DNA, P_RNA
    tips = kh_prona[kh_prona["node_type"] == "tip"].copy()
    ancestors = kh_prona[kh_prona["node_type"] == "ancestor"].copy()

    domain_colors = {
        "Archaea": "#F39C12",
        "Bacteria": "#3498DB",
        "Eukaryota": "#27AE60",
    }
    for domain, color in domain_colors.items():
        subset = tips[tips["domain_of_life"] == domain]
        if len(subset) > 0:
            ax.scatter(subset["P_DNA"], subset["P_RNA"],
                      c=color, s=25, alpha=0.7, edgecolors="black",
                      linewidth=0.3, label=f"Tips: {domain} (n={len(subset)})",
                      zorder=3)

    ax.scatter(ancestors["P_DNA"], ancestors["P_RNA"],
              c="#999999", s=40, alpha=0.8, edgecolors="black",
              linewidth=0.5, marker="D", label=f"Ancestors (n={len(ancestors)})",
              zorder=4)

    root = ancestors[ancestors["node_id"].str.contains("Node1|root", case=False, na=False)]
    if len(root) == 0 and len(ancestors) > 0:
        root = ancestors.head(1)
    if len(root) > 0:
        r = root.iloc[0]
        ax.annotate(f"Root\nP_RNA={r['P_RNA']:.3f}\nP_DNA={r['P_DNA']:.3f}",
                    xy=(r["P_DNA"], r["P_RNA"]),
                    xytext=(r["P_DNA"] + 0.1, r["P_RNA"] - 0.1),
                    fontsize=6, fontstyle="italic",
                    arrowprops=dict(arrowstyle="->", color="black", lw=0.8))

    ax.plot([0, 1], [0, 1], "--", color="gray", linewidth=0.5, alpha=0.5)
    ax.text(0.05, 0.95, "RNA-biased", fontsize=7, color=RNA_COLOR, alpha=0.5)
    ax.text(0.85, 0.05, "DNA-biased", fontsize=7, color=DNA_COLOR, alpha=0.5, ha="right")

    ax.set_xlabel("ProNA2020 P(DNA-binding)")
    ax.set_ylabel("ProNA2020 P(RNA-binding)")
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=6, frameon=True, framealpha=0.8)
    ax.set_title("KH domain: sequence-level binding predictions", fontsize=9, fontweight="bold")

    fig.tight_layout()
    save_fig(fig, "fig5_kh_phylogeny")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("generating publication figures...")
    print(f"output: {FIGDIR}/")

    for name, func in [
        ("Figure 2", make_figure_2),
        ("Figure 3", make_figure_3),
        ("Figure 1", make_figure_1),
        ("Figure 4", make_figure_4),
        ("Figure 5", make_figure_5),
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
