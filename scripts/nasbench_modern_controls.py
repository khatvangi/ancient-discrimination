#!/usr/bin/env python3
"""
NAS-Bench control experiment: modern specialist proteins.
validates that SI cleanly separates RNA vs DNA specialists,
then contrasts with LUCA domain DI distribution.

category A: RNA-only specialists (expect SI > 0.5)
category B: DNA-only specialists (expect SI < 0.5)
category C: known dual-binders (expect low DI)

output:
  results/nasbench_modern_controls.tsv — per-protein SI values
  figures/fig1b_modern_vs_luca.png     — control panel for Figure 1
"""

import sys, os, time, gc
import numpy as np
import pandas as pd
from scipy import stats

# add scripts dir to path for nasbench imports
sys.path.insert(0, os.path.dirname(__file__))
from nasbench import compute_si, compute_si_split, fetch

BASE = '/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination'

# ═══════════════════════════════════════════════════════════════════════
# define control proteins
# ═══════════════════════════════════════════════════════════════════════

# category A: RNA-only specialists
# each tuple: (name, pdb_id, description)
RNA_SPECIALISTS = [
    ('U1A RRM',          '1URN', 'U1A spliceosomal protein + U1 snRNA hairpin'),
    ('MS2 coat',         '1AQ3', 'MS2 coat protein + RNA stem-loop'),
    ('PUF/Pumilio',      '3Q0P', 'human PUF domain + NRE RNA'),
    ('PAZ domain',       '1SI3', 'Argonaute PAZ + siRNA'),
    ('dsRBD (Xlrbpa)',   '1DI2', 'Xenopus dsRBD + dsRNA'),
    ('Pumilio repeat',   '3BSX', 'C. elegans FBF-2 + RNA'),
    ('Fox-1 RRM',        '2ERR', 'Fox-1 RRM + UGCAUGU RNA'),
    ('hnRNP A1 RRM',     '2UP1', 'hnRNP A1 RRM + telomeric RNA'),
    ('SRSF1 RRM',        '6HPJ', 'SRSF1 RRM + RNA ESE'),
    ('L1 stalk rRNA',    '4V4Q', 'ribosomal L1 + rRNA fragment'),
]

# category B: DNA-only specialists
DNA_SPECIALISTS = [
    ('Lambda repressor',  '1LMB', 'Lambda cI repressor + operator DNA'),
    ('p53 core',          '1TSR', 'p53 core domain + response element DNA'),
    ('Zif268 ZnF',        '1ZAA', 'Zif268 zinc finger + DNA'),
    ('EcoRI',             '1ERI', 'EcoRI restriction enzyme + cognate DNA'),
    ('BamHI',             '1BHM', 'BamHI restriction enzyme + DNA'),
    ('TBP',               '1CDW', 'TATA-binding protein + TATA box DNA'),
    ('Lac repressor',     '1EFA', 'Lac repressor + operator DNA'),
    ('Engrailed HD',      '1HDD', 'Engrailed homeodomain + DNA'),
    ('Cro repressor',     '3CRO', 'Cro repressor + operator DNA'),
    ('EcoRV',             '1RVA', 'EcoRV restriction enzyme + DNA'),
]

# category C: known dual-binders (paired RNA + DNA structures)
# each tuple: (name, rna_pdb, dna_pdb, description)
DUAL_BINDERS = [
    ('TDP-43 RRM1', '4BS2', '4IUF', 'TDP-43 RRM1: RNA (UG repeat) vs DNA (TG repeat)'),
    ('CspA/CSD',    '3PF5', '2HAX', 'CspA cold shock domain: RNA vs DNA'),
    ('Lin28 ZnF+CSD', '4A76', '4A4I', 'Lin28 zinc knuckle+CSD: pre-let-7 RNA vs ssDNA'),
]


def compute_si_safe(pdb_id, target_na=None):
    """
    download PDB and compute SI.
    if target_na specified ('RNA' or 'DNA'), use split to get that fraction.
    returns dict with SI info or None on failure.
    """
    try:
        pdb_path = fetch(pdb_id)
        if pdb_path is None:
            print(f'    FAIL: could not download {pdb_id}')
            return None

        # first try split to detect mixed structures
        split = compute_si_split(pdb_path, pdb_id)

        if target_na is not None:
            target_key = target_na.lower()  # 'rna' or 'dna'
            if split and split.get(target_key) is not None:
                result = split[target_key]
                result['method'] = f'split({target_na})'
                return result
            # fallback: try standard compute_si
            result = compute_si(pdb_path, pdb_id)
            if result and result.get('n_total', 0) > 0:
                result['method'] = 'pure'
                return result
            return None
        else:
            # no target specified, use standard compute_si
            result = compute_si(pdb_path, pdb_id)
            if result and result.get('n_total', 0) > 0:
                result['method'] = 'pure'
                return result
            return None

    except Exception as e:
        print(f'    ERROR on {pdb_id}: {e}')
        return None
    finally:
        gc.collect()


def run_category(proteins, category, target_na=None):
    """
    run SI computation on a list of (name, pdb, desc) tuples.
    returns list of result dicts.
    """
    results = []
    for name, pdb_id, desc in proteins:
        print(f'  [{category}] {name} ({pdb_id})...')
        r = compute_si_safe(pdb_id, target_na=target_na)
        if r is not None:
            results.append({
                'name': name,
                'pdb_id': pdb_id,
                'category': category,
                'target_na': target_na or 'auto',
                'description': desc,
                'SI': r.get('SI', np.nan),
                'bb_frac': r.get('bb_frac', np.nan),
                'base_frac': r.get('base_frac', np.nan),
                'oh_frac': r.get('oh_frac', np.nan),
                'n_contacts': r.get('n_total', 0),
                'na_type': r.get('na_type', '?'),
                'method': r.get('method', '?'),
            })
            print(f'    SI = {r.get("SI", "N/A"):.4f}, '
                  f'n = {r.get("n_total", 0)}, '
                  f'na_type = {r.get("na_type", "?")}')
        else:
            print(f'    FAILED — no contacts or download error')
            results.append({
                'name': name,
                'pdb_id': pdb_id,
                'category': category,
                'target_na': target_na or 'auto',
                'description': desc,
                'SI': np.nan,
                'bb_frac': np.nan,
                'base_frac': np.nan,
                'oh_frac': np.nan,
                'n_contacts': 0,
                'na_type': 'FAIL',
                'method': 'FAIL',
            })
        time.sleep(0.3)  # be gentle with RCSB
    return results


def run_dual_binders(duals):
    """
    run paired SI on dual-binder proteins.
    returns list of result dicts with DI computation.
    """
    results = []
    for name, rna_pdb, dna_pdb, desc in duals:
        print(f'  [DUAL] {name}: RNA={rna_pdb}, DNA={dna_pdb}')

        r_rna = compute_si_safe(rna_pdb, target_na='RNA')
        r_dna = compute_si_safe(dna_pdb, target_na='DNA')

        if r_rna is not None and r_dna is not None:
            rna_si = r_rna.get('SI', np.nan)
            dna_si = r_dna.get('SI', np.nan)
            delta_si = rna_si - dna_si
            di = abs(delta_si) / (rna_si + dna_si) if (rna_si + dna_si) > 0 else np.nan
            print(f'    RNA SI = {rna_si:.4f}, DNA SI = {dna_si:.4f}, '
                  f'ΔSI = {delta_si:+.4f}, DI = {di:.4f}')

            results.append({
                'name': name,
                'rna_pdb': rna_pdb,
                'dna_pdb': dna_pdb,
                'category': 'dual',
                'description': desc,
                'rna_SI': rna_si,
                'dna_SI': dna_si,
                'delta_SI': delta_si,
                'DI': di,
                'rna_contacts': r_rna.get('n_total', 0),
                'dna_contacts': r_dna.get('n_total', 0),
                'rna_method': r_rna.get('method', '?'),
                'dna_method': r_dna.get('method', '?'),
            })
        else:
            fail_detail = []
            if r_rna is None: fail_detail.append(f'RNA({rna_pdb})')
            if r_dna is None: fail_detail.append(f'DNA({dna_pdb})')
            print(f'    FAILED: {", ".join(fail_detail)}')
            results.append({
                'name': name,
                'rna_pdb': rna_pdb,
                'dna_pdb': dna_pdb,
                'category': 'dual',
                'description': desc,
                'rna_SI': r_rna.get('SI', np.nan) if r_rna else np.nan,
                'dna_SI': r_dna.get('SI', np.nan) if r_dna else np.nan,
                'delta_SI': np.nan,
                'DI': np.nan,
                'rna_contacts': r_rna.get('n_total', 0) if r_rna else 0,
                'dna_contacts': r_dna.get('n_total', 0) if r_dna else 0,
                'rna_method': r_rna.get('method', 'FAIL') if r_rna else 'FAIL',
                'dna_method': r_dna.get('method', 'FAIL') if r_dna else 'FAIL',
            })
        time.sleep(0.3)
    return results


# ═══════════════════════════════════════════════════════════════════════
# statistics and figure
# ═══════════════════════════════════════════════════════════════════════

def compute_statistics(rna_results, dna_results, dual_results, luca_df):
    """compute all statistical comparisons."""
    rna_si = [r['SI'] for r in rna_results if not np.isnan(r['SI'])]
    dna_si = [r['SI'] for r in dna_results if not np.isnan(r['SI'])]

    print('\n' + '=' * 70)
    print('STATISTICAL ANALYSIS')
    print('=' * 70)

    # 1. SI distributions
    print(f'\nRNA specialists (n={len(rna_si)}):')
    print(f'  mean SI = {np.mean(rna_si):.4f} ± {np.std(rna_si):.4f}')
    print(f'  range: [{min(rna_si):.4f}, {max(rna_si):.4f}]')

    print(f'\nDNA specialists (n={len(dna_si)}):')
    print(f'  mean SI = {np.mean(dna_si):.4f} ± {np.std(dna_si):.4f}')
    print(f'  range: [{min(dna_si):.4f}, {max(dna_si):.4f}]')

    # 2. Mann-Whitney U: RNA SI vs DNA SI
    u_stat, p_rna_dna = stats.mannwhitneyu(rna_si, dna_si, alternative='greater')
    print(f'\nMann-Whitney U (RNA SI > DNA SI):')
    print(f'  U = {u_stat:.1f}, p = {p_rna_dna:.2e}')
    print(f'  effect size (rank-biserial r) = {1 - 2*u_stat/(len(rna_si)*len(dna_si)):.4f}')

    # 3. dual binder DI
    dual_di = [r['DI'] for r in dual_results if not np.isnan(r.get('DI', np.nan))]
    if dual_di:
        print(f'\nDual-binder DI (n={len(dual_di)}):')
        for r in dual_results:
            if not np.isnan(r.get('DI', np.nan)):
                print(f'  {r["name"]}: DI = {r["DI"]:.4f} (ΔSI = {r["delta_SI"]:+.4f})')

    # 4. LUCA DI vs hypothetical modern specialist DI
    # compute pseudo-DI: pair each RNA specialist with each DNA specialist
    pseudo_di = []
    for rsi in rna_si:
        for dsi in dna_si:
            if (rsi + dsi) > 0:
                pseudo_di.append(abs(rsi - dsi) / (rsi + dsi))

    luca_di = luca_df['DI'].dropna().values

    print(f'\nPseudo-DI for modern specialist pairs (n={len(pseudo_di)}):')
    print(f'  mean = {np.mean(pseudo_di):.4f} ± {np.std(pseudo_di):.4f}')
    print(f'  range: [{min(pseudo_di):.4f}, {max(pseudo_di):.4f}]')

    print(f'\nLUCA DI (n={len(luca_di)}):')
    print(f'  mean = {np.mean(luca_di):.4f} ± {np.std(luca_di):.4f}')
    print(f'  median = {np.median(luca_di):.4f}')

    # 5. Mann-Whitney: modern pseudo-DI vs LUCA DI
    u_luca, p_luca = stats.mannwhitneyu(pseudo_di, luca_di, alternative='greater')
    print(f'\nMann-Whitney U (modern pseudo-DI > LUCA DI):')
    print(f'  U = {u_luca:.1f}, p = {p_luca:.2e}')

    # 6. also direct comparison: LUCA DI vs dual-binder DI
    if dual_di:
        u_dual, p_dual = stats.mannwhitneyu(luca_di, dual_di, alternative='two-sided')
        print(f'\nMann-Whitney U (LUCA DI vs dual DI, two-sided):')
        print(f'  U = {u_dual:.1f}, p = {p_dual:.4f}')
        print(f'  (LUCA DI resembles modern dual-binders: {p_dual > 0.05})')

    return {
        'rna_si': rna_si,
        'dna_si': dna_si,
        'pseudo_di': pseudo_di,
        'luca_di': luca_di,
        'dual_di': dual_di,
        'p_rna_dna': p_rna_dna,
        'p_luca_vs_modern': p_luca,
    }


def make_figure(stat_dict, rna_results, dna_results, dual_results, luca_df):
    """create the control comparison figure."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.linewidth': 0.8,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

    fig, axes = plt.subplots(1, 3, figsize=(14, 5),
                              gridspec_kw={'width_ratios': [1, 1, 1]})

    # ── panel A: modern SI distributions ──────────────────────────────
    ax = axes[0]

    rna_si = stat_dict['rna_si']
    dna_si = stat_dict['dna_si']

    # strip plot with jitter
    np.random.seed(42)
    jitter = 0.08

    # RNA specialists
    y_rna = np.random.normal(0, jitter, len(rna_si))
    ax.scatter(rna_si, y_rna + 1, c='#4878CF', s=55, alpha=0.85,
               edgecolors='white', linewidths=0.5, zorder=5)

    # DNA specialists
    y_dna = np.random.normal(0, jitter, len(dna_si))
    ax.scatter(dna_si, y_dna + 0, c='#C44E52', s=55, alpha=0.85,
               edgecolors='white', linewidths=0.5, zorder=5)

    # means with error bars
    ax.errorbar(np.mean(rna_si), 1, xerr=np.std(rna_si),
                fmt='D', color='#2d4fa5', markersize=7, capsize=4,
                capthick=1.2, linewidth=1.2, zorder=6)
    ax.errorbar(np.mean(dna_si), 0, xerr=np.std(dna_si),
                fmt='D', color='#a03035', markersize=7, capsize=4,
                capthick=1.2, linewidth=1.2, zorder=6)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(['DNA\nspecialists', 'RNA\nspecialists'], fontsize=9)
    ax.set_xlabel('Specificity Index (SI)')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.5, 1.5)
    ax.set_title('A. Modern specialists', fontweight='bold')

    # p-value annotation
    p = stat_dict['p_rna_dna']
    p_text = f'p = {p:.1e}' if p < 0.001 else f'p = {p:.4f}'
    ax.annotate(f'Mann-Whitney\n{p_text}', xy=(0.5, 0.5),
                xycoords='axes fraction', ha='center', fontsize=7,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#f7f7f7',
                          edgecolor='#cccccc'))

    # ── panel B: LUCA SI (RNA vs DNA, paired) ─────────────────────────
    ax = axes[1]

    # get valid LUCA pairs
    valid = luca_df.dropna(subset=['rna_SI', 'dna_SI', 'DI']).copy()

    # scatter: x = RNA SI, y = DNA SI
    for _, row in valid.iterrows():
        color = '#6ACC65' if row['DI'] < 0.10 else (
                '#B47CC7' if row['DI'] < 0.25 else '#C44E52')
        ax.scatter(row['rna_SI'], row['dna_SI'], c=color, s=55,
                   edgecolors='white', linewidths=0.5, zorder=5, alpha=0.85)

    # diagonal line (perfect generalism)
    ax.plot([0, 1], [0, 1], '--', color='#999999', linewidth=0.8, alpha=0.6)

    # add modern specialist means as reference markers
    ax.axvline(np.mean(rna_si), color='#4878CF', linestyle=':', linewidth=0.8,
               alpha=0.5, label=f'Modern RNA mean ({np.mean(rna_si):.2f})')
    ax.axhline(np.mean(dna_si), color='#C44E52', linestyle=':', linewidth=0.8,
               alpha=0.5, label=f'Modern DNA mean ({np.mean(dna_si):.2f})')

    ax.set_xlabel('RNA complex SI')
    ax.set_ylabel('DNA complex SI')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.set_title('B. LUCA domains (paired)', fontweight='bold')
    ax.set_aspect('equal')

    # legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#6ACC65',
               markersize=7, label='Generalist (DI<0.10)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#B47CC7',
               markersize=7, label='Moderate (0.10-0.25)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#C44E52',
               markersize=7, label='Specialist (DI>0.25)'),
        Line2D([0], [0], color='#999999', linestyle='--', linewidth=0.8,
               label='Perfect generalism'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=6.5,
              frameon=True, framealpha=0.9, edgecolor='#cccccc')

    # ── panel C: DI comparison ────────────────────────────────────────
    ax = axes[2]

    luca_di = stat_dict['luca_di']
    pseudo_di = stat_dict['pseudo_di']
    dual_di = stat_dict['dual_di']

    # violin / box comparison
    positions = [0, 1, 2]
    data = [luca_di, dual_di if dual_di else [np.nan], pseudo_di]
    labels = ['LUCA\ndomains', 'Modern\ndual-binders', 'Modern\nspecialists\n(cross-paired)']
    colors_box = ['#6ACC65', '#B47CC7', '#C44E52']

    # use strip plot + box
    for i, (d, c) in enumerate(zip(data, colors_box)):
        d_arr = np.array([x for x in d if not np.isnan(x)])
        if len(d_arr) > 0:
            # box
            bp = ax.boxplot([d_arr], positions=[i], widths=0.4,
                            patch_artist=True, showfliers=False,
                            medianprops=dict(color='black', linewidth=1.5))
            bp['boxes'][0].set_facecolor(c)
            bp['boxes'][0].set_alpha(0.3)
            bp['boxes'][0].set_edgecolor(c)

            # strip jitter (limit to 30 points for pseudo-DI to avoid clutter)
            if len(d_arr) > 30:
                d_plot = np.random.choice(d_arr, 30, replace=False)
            else:
                d_plot = d_arr
            jit = np.random.normal(0, 0.06, len(d_plot))
            ax.scatter(jit + i, d_plot, c=c, s=20, alpha=0.6,
                       edgecolors='white', linewidths=0.3, zorder=5)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel('Discrimination Index (DI)')
    ax.set_title('C. DI comparison', fontweight='bold')

    # threshold lines
    ax.axhline(0.10, color='#555555', linestyle='--', linewidth=0.6, alpha=0.5)
    ax.axhline(0.25, color='#555555', linestyle='--', linewidth=0.6, alpha=0.5)
    ax.text(2.6, 0.10, 'DI=0.10', fontsize=6, color='#555555', va='center')
    ax.text(2.6, 0.25, 'DI=0.25', fontsize=6, color='#555555', va='center')

    # p-value annotation
    p_luca = stat_dict['p_luca_vs_modern']
    p_text = f'p = {p_luca:.1e}' if p_luca < 0.001 else f'p = {p_luca:.4f}'
    ax.annotate(f'LUCA vs modern specialist\n{p_text}',
                xy=(1.5, max(pseudo_di) * 0.85 if pseudo_di else 0.5),
                ha='center', fontsize=7,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#f7f7f7',
                          edgecolor='#cccccc'))

    plt.suptitle('NAS-Bench validation: modern specialists vs LUCA domains',
                 fontweight='bold', fontsize=12, y=1.02)
    plt.tight_layout()

    figdir = f'{BASE}/figures'
    fig.savefig(f'{figdir}/fig1b_modern_vs_luca.png')
    fig.savefig(f'{figdir}/fig1b_modern_vs_luca.pdf')
    plt.close()
    print(f'\nfigure saved to {figdir}/fig1b_modern_vs_luca.{{png,pdf}}')


# ═══════════════════════════════════════════════════════════════════════
# main
# ═══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print('=' * 70)
    print('NAS-BENCH CONTROL EXPERIMENT: modern specialist proteins')
    print('=' * 70)

    # ── run category A: RNA specialists ───────────────────────────────
    print('\n--- CATEGORY A: RNA-only specialists ---')
    rna_results = run_category(RNA_SPECIALISTS, 'RNA_specialist', target_na='RNA')

    # ── run category B: DNA specialists ───────────────────────────────
    print('\n--- CATEGORY B: DNA-only specialists ---')
    dna_results = run_category(DNA_SPECIALISTS, 'DNA_specialist', target_na='DNA')

    # ── run category C: dual binders ──────────────────────────────────
    print('\n--- CATEGORY C: known dual-binders ---')
    dual_results = run_dual_binders(DUAL_BINDERS)

    # ── save results TSV ──────────────────────────────────────────────
    outpath = f'{BASE}/results/nasbench_modern_controls.tsv'

    # singles
    singles_df = pd.DataFrame(rna_results + dna_results)
    singles_df.to_csv(outpath, sep='\t', index=False)
    print(f'\nsingles results saved to {outpath}')

    # duals as separate file
    duals_df = pd.DataFrame(dual_results)
    dual_outpath = f'{BASE}/results/nasbench_modern_duals.tsv'
    duals_df.to_csv(dual_outpath, sep='\t', index=False)
    print(f'dual results saved to {dual_outpath}')

    # ── load LUCA data ────────────────────────────────────────────────
    luca_df = pd.read_csv(f'{BASE}/results/nasbench_full_luca.tsv', sep='\t')
    # filter to complete, deduplicated
    luca_df = luca_df[luca_df['generalism_mode'] != 'incomplete']
    luca_df = luca_df[luca_df['DI'].notna()]
    luca_df = luca_df[luca_df['pfam_id'] != 'PF10996']  # dedup

    # ── statistics ────────────────────────────────────────────────────
    stat_dict = compute_statistics(rna_results, dna_results,
                                   dual_results, luca_df)

    # ── figure ────────────────────────────────────────────────────────
    print('\ngenerating comparison figure...')
    make_figure(stat_dict, rna_results, dna_results, dual_results, luca_df)

    print('\ndone.')
