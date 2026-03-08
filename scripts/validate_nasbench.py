#!/usr/bin/env python3
"""
validate_nasbench.py — formal NAS-Bench metric validation
==========================================================

addresses reviewer concern: "no rigorous benchmark section that defines
expected null distributions, confidence intervals, family-level uncertainty,
robustness to PDB choice, or sensitivity to the exact DI cutoff."

outputs:
  results/nasbench_validation.tsv — per-family bootstrap CIs + sensitivity
  stdout — full validation report

formulas:
  SI  = (n_base + n_2oh) / (n_base + n_2oh + n_bb)
      where contacts are protein-NA atom pairs within 4.0 A cutoff:
        n_base  = contacts to nucleobase atoms
        n_2oh   = contacts to ribose 2'-OH (RNA-specific)
        n_bb    = contacts to phosphodiester backbone atoms (P, OP1, OP2, O5', O3')
      sugar-ring contacts (C1'-C5', O4') are excluded from the SI ratio
      as they are neither fully specific nor fully generic.

  DI  = |SI_rna - SI_dna| / (SI_rna + SI_dna)
      normalized absolute difference in specificity index between the
      RNA-bound and DNA-bound co-crystal structures of the same domain.
      range: [0, 1]. DI=0 means identical contact profiles. DI=1 means
      one SI is zero while the other is nonzero.

  classification thresholds (default):
      DI < 0.10  → generalist (does not distinguish RNA from DNA)
      DI 0.10-0.25 → moderate
      DI > 0.25  → specialist (discriminates RNA from DNA)

validation tests:
  1. bootstrap 95% CI on DI for each LUCA family (resample residue contacts)
  2. threshold sensitivity: fraction of families changing category at ±0.05
  3. negative controls: known RNA-specialists and DNA-specialists from
     modern controls data — verify high DI for specialists, low for generalists
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASR_DIR = os.path.join(BASE_DIR, "results", "asr")

np.random.seed(42)
N_BOOT = 1000


def load_contacts(filepath):
    """load per-residue contact file and return arrays of (n_bb, n_base, n_2oh).

    each row is one protein residue that contacts nucleic acid.
    we return the per-residue contact counts as the bootstrap unit.
    """
    df = pd.read_csv(filepath, sep='\t')
    # columns vary slightly between alba/kh files and pf* files
    required = ['n_bb', 'n_base', 'n_2oh']
    for col in required:
        if col not in df.columns:
            return None
    return df[required].values  # shape (n_residues, 3)


def compute_si_from_counts(bb, base, oh):
    """compute SI from aggregate contact counts."""
    spec = base + oh
    gen = bb
    if (spec + gen) == 0:
        return float('nan')
    return spec / (spec + gen)


def compute_di(si_rna, si_dna):
    """compute discrimination index."""
    if (si_rna + si_dna) == 0:
        return float('nan')
    return abs(si_rna - si_dna) / (si_rna + si_dna)


def bootstrap_di(rna_contacts, dna_contacts, n_boot=N_BOOT):
    """bootstrap DI by resampling residue-level contacts with replacement.

    for each iteration:
      1. resample RNA residues → recompute SI_rna
      2. resample DNA residues → recompute SI_dna
      3. compute DI from resampled SI values

    returns: (mean_di, ci_lo, ci_hi, all_di_values)
    """
    n_rna = len(rna_contacts)
    n_dna = len(dna_contacts)
    dis = np.empty(n_boot)

    for b in range(n_boot):
        # resample residues with replacement
        rna_idx = np.random.randint(0, n_rna, size=n_rna)
        dna_idx = np.random.randint(0, n_dna, size=n_dna)

        rna_samp = rna_contacts[rna_idx]
        dna_samp = dna_contacts[dna_idx]

        # aggregate resampled counts
        rna_bb = rna_samp[:, 0].sum()
        rna_base = rna_samp[:, 1].sum()
        rna_oh = rna_samp[:, 2].sum()

        dna_bb = dna_samp[:, 0].sum()
        dna_base = dna_samp[:, 1].sum()
        dna_oh = dna_samp[:, 2].sum()

        si_rna = compute_si_from_counts(rna_bb, rna_base, rna_oh)
        si_dna = compute_si_from_counts(dna_bb, dna_base, dna_oh)
        dis[b] = compute_di(si_rna, si_dna)

    # remove any nans
    dis = dis[~np.isnan(dis)]
    if len(dis) == 0:
        return float('nan'), float('nan'), float('nan'), np.array([])

    ci_lo = np.percentile(dis, 2.5)
    ci_hi = np.percentile(dis, 97.5)
    return np.mean(dis), ci_lo, ci_hi, dis


def find_contact_files(pfam_id):
    """find RNA and DNA contact files for a given pfam family.

    returns: (rna_file, dna_file) or (None, None) if not found.
    """
    pfam_lower = pfam_id.lower()
    family_dir = os.path.join(ASR_DIR, pfam_lower)
    if not os.path.isdir(family_dir):
        return None, None

    rna_file = None
    dna_file = None
    for f in os.listdir(family_dir):
        if f.endswith('_rna_contacts.tsv') or f.endswith('_rna_residue_contacts.tsv'):
            rna_file = os.path.join(family_dir, f)
        if f.endswith('_dna_contacts.tsv') or f.endswith('_dna_residue_contacts.tsv'):
            dna_file = os.path.join(family_dir, f)

    return rna_file, dna_file


def classify_di(di, lo_thresh=0.10, hi_thresh=0.25):
    """classify DI into generalist / moderate / specialist."""
    if np.isnan(di):
        return "unknown"
    if di < lo_thresh:
        return "generalist"
    elif di < hi_thresh:
        return "moderate"
    else:
        return "specialist"


def main():
    print("=" * 80)
    print("NAS-Bench Metric Validation Report")
    print("=" * 80)

    # ===================================================================
    # section 0: formal metric definitions
    # ===================================================================
    print("""
FORMAL METRIC DEFINITIONS
--------------------------
  SI  = (n_base + n_2oh) / (n_base + n_2oh + n_bb)
  DI  = |SI_rna - SI_dna| / (SI_rna + SI_dna)

  classification:
    DI < 0.10  → generalist
    DI 0.10-0.25 → moderate
    DI > 0.25  → specialist

  bootstrap: resample per-residue contacts (n=1000 iterations), recompute
  SI and DI from resampled aggregate counts. report 95% percentile CI.
""")

    # ===================================================================
    # section 1: load LUCA family data and bootstrap CIs
    # ===================================================================
    luca_path = os.path.join(BASE_DIR, "results", "nasbench_full_luca.tsv")
    luca = pd.read_csv(luca_path, sep='\t')
    # keep only complete rows (have both RNA and DNA SI)
    complete = luca[luca['generalism_mode'] != 'incomplete'].copy()

    print(f"SECTION 1: Bootstrap 95% CI on DI ({N_BOOT} iterations, residue resampling)")
    print("-" * 80)
    print(f"{'Pfam':<10} {'DI_obs':>7} {'DI_boot':>8} {'CI_lo':>7} {'CI_hi':>7} {'CI_width':>8} {'Cat_obs':<20} {'CI_spans_thresh'}")
    print("-" * 80)

    validation_rows = []

    for _, row in complete.iterrows():
        pfam_id = row['pfam_id']
        di_obs = row['DI']
        cat_obs = row['generalism_mode']

        rna_file, dna_file = find_contact_files(pfam_id)

        if rna_file is None or dna_file is None:
            # try alba / kh special directories
            if pfam_id == 'PF01918':
                rna_file = os.path.join(ASR_DIR, 'alba', 'alba_rna_residue_contacts.tsv')
                dna_file = os.path.join(ASR_DIR, 'alba', 'alba_dna_residue_contacts.tsv')
            elif pfam_id == 'PF00013':
                # kh family
                rna_file = os.path.join(ASR_DIR, 'kh', 'kh_rna_residue_contacts.tsv')
                dna_file = os.path.join(ASR_DIR, 'kh', 'kh_dna_residue_contacts.tsv')

        if rna_file is None or dna_file is None:
            print(f"{pfam_id:<10} {di_obs:>7.4f} {'N/A':>8} {'N/A':>7} {'N/A':>7} {'N/A':>8} {cat_obs:<20} no contact files")
            validation_rows.append({
                'pfam_id': pfam_id, 'DI_observed': di_obs,
                'DI_boot_mean': '', 'CI_lo': '', 'CI_hi': '', 'CI_width': '',
                'category_observed': cat_obs,
                'CI_spans_0.10': '', 'CI_spans_0.25': '',
                'n_rna_residues': '', 'n_dna_residues': '',
            })
            continue

        rna_contacts = load_contacts(rna_file)
        dna_contacts = load_contacts(dna_file)

        if rna_contacts is None or dna_contacts is None:
            print(f"{pfam_id:<10} {di_obs:>7.4f} {'ERR':>8} {'ERR':>7} {'ERR':>7} {'ERR':>8} {cat_obs:<20} contact parse error")
            continue

        if len(rna_contacts) < 2 or len(dna_contacts) < 2:
            print(f"{pfam_id:<10} {di_obs:>7.4f} {'<2res':>8} {'':>7} {'':>7} {'':>8} {cat_obs:<20} too few residues")
            continue

        di_mean, ci_lo, ci_hi, all_dis = bootstrap_di(rna_contacts, dna_contacts)
        ci_width = ci_hi - ci_lo

        # check if CI spans either threshold
        spans_010 = "YES" if (ci_lo < 0.10 < ci_hi) else "no"
        spans_025 = "YES" if (ci_lo < 0.25 < ci_hi) else "no"
        spans_note = ""
        if spans_010 == "YES" or spans_025 == "YES":
            spans_note = f"spans {'0.10' if spans_010 == 'YES' else ''}{'+' if spans_010 == 'YES' and spans_025 == 'YES' else ''}{'0.25' if spans_025 == 'YES' else ''}"
        else:
            spans_note = "stable"

        print(f"{pfam_id:<10} {di_obs:>7.4f} {di_mean:>8.4f} {ci_lo:>7.4f} {ci_hi:>7.4f} {ci_width:>8.4f} {cat_obs:<20} {spans_note}")

        validation_rows.append({
            'pfam_id': pfam_id, 'DI_observed': round(di_obs, 4),
            'DI_boot_mean': round(di_mean, 4),
            'CI_lo': round(ci_lo, 4), 'CI_hi': round(ci_hi, 4),
            'CI_width': round(ci_width, 4),
            'category_observed': cat_obs,
            'CI_spans_0.10': spans_010, 'CI_spans_0.25': spans_025,
            'n_rna_residues': len(rna_contacts),
            'n_dna_residues': len(dna_contacts),
        })

    # ===================================================================
    # section 2: threshold sensitivity analysis
    # ===================================================================
    print(f"\n{'=' * 80}")
    print("SECTION 2: Threshold Sensitivity Analysis")
    print("=" * 80)
    print("""
  test: what fraction of families change DI category if thresholds shift ±0.05?
  default thresholds: generalist < 0.10, moderate 0.10-0.25, specialist > 0.25
""")

    # get observed DI values
    di_values = complete['DI'].values
    n_families = len(di_values)

    # define threshold variants
    threshold_variants = [
        ("strict",  0.05, 0.20),
        ("default", 0.10, 0.25),
        ("lenient", 0.15, 0.30),
    ]

    print(f"{'Variant':<10} {'Lo_thr':>7} {'Hi_thr':>7} {'n_gen':>6} {'n_mod':>6} {'n_spec':>7} {'%_gen':>6} {'%_mod':>6} {'%_spec':>7}")
    print("-" * 70)

    # store default classification
    default_cats = [classify_di(d, 0.10, 0.25) for d in di_values]

    for name, lo, hi in threshold_variants:
        cats = [classify_di(d, lo, hi) for d in di_values]
        n_gen = cats.count("generalist")
        n_mod = cats.count("moderate")
        n_spec = cats.count("specialist")
        print(f"{name:<10} {lo:>7.2f} {hi:>7.2f} {n_gen:>6} {n_mod:>6} {n_spec:>7} {100*n_gen/n_families:>5.0f}% {100*n_mod/n_families:>5.0f}% {100*n_spec/n_families:>6.0f}%")

    # count families that change category
    changed_strict = sum(1 for d, c in zip(di_values, default_cats) if classify_di(d, 0.05, 0.20) != c)
    changed_lenient = sum(1 for d, c in zip(di_values, default_cats) if classify_di(d, 0.15, 0.30) != c)

    print(f"\n  families changing category (default → strict):  {changed_strict}/{n_families} ({100*changed_strict/n_families:.0f}%)")
    print(f"  families changing category (default → lenient): {changed_lenient}/{n_families} ({100*changed_lenient/n_families:.0f}%)")

    # find which specific families are threshold-sensitive
    print(f"\n  threshold-sensitive families (near boundary):")
    for _, row in complete.iterrows():
        d = row['DI']
        pfam = row['pfam_id']
        # near 0.10 boundary
        if 0.05 <= d <= 0.15:
            print(f"    {pfam} DI={d:.4f} — near generalist/moderate boundary (0.10)")
        # near 0.25 boundary
        elif 0.20 <= d <= 0.30:
            print(f"    {pfam} DI={d:.4f} — near moderate/specialist boundary (0.25)")

    # ===================================================================
    # section 3: negative-control validation with modern specialists
    # ===================================================================
    print(f"\n{'=' * 80}")
    print("SECTION 3: Negative-Control Validation (Modern Specialists)")
    print("=" * 80)
    print("""
  rationale: if SI and DI are valid metrics, then:
    - known RNA-specialists should have high SI on RNA (>0.80)
    - known DNA-specialists should have varied SI on DNA (sequence-specific
      binders contact bases; structural binders contact backbone)
    - LUCA generalists should have low DI (<0.10), i.e., similar contact
      profiles regardless of whether they bind RNA or DNA
""")

    controls_path = os.path.join(BASE_DIR, "results", "nasbench_modern_controls.tsv")
    controls = pd.read_csv(controls_path, sep='\t')

    # exclude failed entries
    controls = controls[controls['method'] != 'FAIL'].copy()
    controls['SI'] = pd.to_numeric(controls['SI'], errors='coerce')

    print(f"{'Name':<20} {'PDB':>5} {'Category':<18} {'NA':>4} {'SI':>7} {'BB%':>6} {'Base%':>7} {'2OH%':>6}")
    print("-" * 80)
    for _, r in controls.iterrows():
        si = r['SI']
        if pd.isna(si):
            continue
        print(f"{r['name']:<20} {r['pdb_id']:>5} {r['category']:<18} {r['na_type']:>4} {si:>7.3f} {r['bb_frac']:>6.1%} {r['base_frac']:>7.1%} {r['oh_frac']:>6.1%}")

    # compute summary statistics by category
    print(f"\n  summary by category:")
    rna_specs = controls[controls['category'] == 'RNA_specialist']
    dna_specs = controls[controls['category'] == 'DNA_specialist']

    if len(rna_specs) > 0:
        rna_si = rna_specs['SI'].dropna()
        print(f"    RNA-specialists (n={len(rna_si)}): mean SI = {rna_si.mean():.3f}, "
              f"range [{rna_si.min():.3f}, {rna_si.max():.3f}]")
    if len(dna_specs) > 0:
        dna_si = dna_specs['SI'].dropna()
        print(f"    DNA-specialists (n={len(dna_si)}): mean SI = {dna_si.mean():.3f}, "
              f"range [{dna_si.min():.3f}, {dna_si.max():.3f}]")

    # compute expected DI for specialist pairs
    # if we paired the top RNA-specialist with each DNA-specialist:
    print(f"\n  expected DI for specialist cross-comparisons:")
    if len(rna_specs) > 0 and len(dna_specs) > 0:
        # pick canonical RNA specialist (U1A = 1URN with highest SI)
        best_rna = rna_specs.loc[rna_specs['SI'].idxmax()]
        print(f"    reference RNA-specialist: {best_rna['name']} ({best_rna['pdb_id']}) SI={best_rna['SI']:.3f}")

        for _, d in dna_specs.iterrows():
            di = compute_di(best_rna['SI'], d['SI'])
            cat = classify_di(di)
            print(f"    vs {d['name']:<16} ({d['pdb_id']}) SI={d['SI']:.3f} → DI={di:.3f} [{cat}]")

    # compare to LUCA families
    print(f"\n  comparison with LUCA ancient domains:")
    luca_complete = complete[complete['DI'].notna()]
    luca_di = luca_complete['DI'].values
    print(f"    LUCA families (n={len(luca_di)}): mean DI = {np.mean(luca_di):.3f}, "
          f"median DI = {np.median(luca_di):.3f}")
    n_gen = sum(1 for d in luca_di if d < 0.10)
    n_mod = sum(1 for d in luca_di if 0.10 <= d < 0.25)
    n_spec = sum(1 for d in luca_di if d >= 0.25)
    print(f"    generalist: {n_gen}/{len(luca_di)} ({100*n_gen/len(luca_di):.0f}%)")
    print(f"    moderate:   {n_mod}/{len(luca_di)} ({100*n_mod/len(luca_di):.0f}%)")
    print(f"    specialist: {n_spec}/{len(luca_di)} ({100*n_spec/len(luca_di):.0f}%)")

    # ===================================================================
    # section 4: save validation results
    # ===================================================================
    out_path = os.path.join(BASE_DIR, "results", "nasbench_validation.tsv")
    out_df = pd.DataFrame(validation_rows)
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f"\n{'=' * 80}")
    print(f"validation results saved to: {out_path}")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
