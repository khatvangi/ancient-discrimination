#!/usr/bin/env python3
"""
nasbench_full_luca.py — scale NAS-Bench to ALL LUCA NA-binding domains.

step 1: query RCSB PDB for each of 87 LUCA domains to find RNA + DNA co-crystals
step 2: for domains with BOTH, find best-resolution isolated structures
step 3: compute SI (or SI_split for mixed structures)
step 4: classify generalism mode and output comprehensive TSV

outputs:
  results/luca_pdb_census.tsv — which domains have RNA/DNA co-crystals
  results/nasbench_full_luca.tsv — SI values + discrimination index for all paired domains
"""

import sys, os, json, time, gc, warnings, tempfile, urllib.request, csv
import numpy as np
from collections import defaultdict

warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nasbench import compute_si, compute_si_split, fetch

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ribosome keywords to filter out whole-ribosome structures
RIBOSOME_KEYWORDS = ["70s", "80s", "ribosome", "30s ribosomal subunit", "50s ribosomal subunit",
                     "40s ribosomal subunit", "60s ribosomal subunit"]


def load_luca_na_domains():
    """load all LUCA+preLUCA domains that are RNA or DNA binding from phase1 census."""
    domains = []
    census_path = os.path.join(BASE_DIR, "results", "phase1_census.tsv")
    with open(census_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['age_class'] not in ('LUCA', 'preLUCA'):
                continue
            is_rna = row['rna_union'] == 'True'
            is_dna = row['dna_binding'] == 'True'
            if is_rna or is_dna:
                domains.append({
                    'pfam_id': row['pfam_id'],
                    'age': row['age_class'],
                    'rna_annotated': is_rna,
                    'dna_annotated': is_dna,
                })
    return domains


def query_pdb_count(pfam_id, na_type, max_results=5):
    """query RCSB for structures with this Pfam + NA type.
    returns list of PDB IDs sorted by resolution (best first).
    """
    attr = ("rcsb_entry_info.polymer_entity_count_RNA" if na_type == "RNA"
            else "rcsb_entry_info.polymer_entity_count_DNA")

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                        "operator": "exact_match",
                        "value": pfam_id
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": attr,
                        "operator": "greater",
                        "value": 0
                    }
                },
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_results},
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            "scoring_strategy": "combined",
        }
    }

    try:
        data = json.dumps(query).encode("utf-8")
        req = urllib.request.Request(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            data=data,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=30) as resp:
            result = json.loads(resp.read().decode())

        total = result.get("total_count", 0)
        entries = [r["identifier"] for r in result.get("result_set", [])]
        return total, entries

    except Exception as e:
        return 0, []


def check_pdb_title(pdb_id):
    """fetch PDB title to check if it's a ribosome structure."""
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        title = data.get("struct", {}).get("title", "").lower()
        return title
    except:
        return ""


def is_ribosome(title):
    """check if a PDB title indicates a whole ribosome structure."""
    t = title.lower()
    for kw in RIBOSOME_KEYWORDS:
        if kw in t:
            return True
    return False


def find_best_structure(pfam_id, na_type, max_candidates=8):
    """find best-resolution non-ribosome structure for this Pfam + NA type.
    returns (pdb_id, title) or (None, None).
    """
    total, candidates = query_pdb_count(pfam_id, na_type, max_results=max_candidates)
    time.sleep(0.35)  # rate limit

    for pdb_id in candidates:
        title = check_pdb_title(pdb_id)
        time.sleep(0.15)

        if is_ribosome(title):
            continue

        return pdb_id, title, total

    # if all candidates are ribosome, return the first one anyway but flag it
    if candidates:
        title = check_pdb_title(candidates[0])
        return candidates[0], title + " [ribosome]", total

    return None, None, total


def compute_si_safe(pdb_id, na_type_target=None):
    """compute SI for a PDB structure, using split if mixed.
    if na_type_target is specified, try to get SI for just that NA type.
    returns dict with SI result or None.
    """
    try:
        path = fetch(pdb_id)

        # first try split to check if mixed
        split = compute_si_split(path, pdb_id)
        has_rna = split['rna'] is not None
        has_dna = split['dna'] is not None

        if na_type_target == "RNA":
            if has_rna and not has_dna:
                # pure RNA — use standard compute_si
                r = compute_si(path, pdb_id)
                r['_method'] = 'pure'
                return r
            elif has_rna:
                # mixed — use RNA portion from split
                r = split['rna']
                r['_method'] = f'split(RNA:{",".join(split["rna_chains"])},DNA:{",".join(split["dna_chains"])})'
                return r
            else:
                return None

        elif na_type_target == "DNA":
            if has_dna and not has_rna:
                r = compute_si(path, pdb_id)
                r['_method'] = 'pure'
                return r
            elif has_dna:
                r = split['dna']
                r['_method'] = f'split(RNA:{",".join(split["rna_chains"])},DNA:{",".join(split["dna_chains"])})'
                return r
            else:
                return None

        else:
            # no target — return combined
            r = compute_si(path, pdb_id)
            r['_method'] = 'combined'
            return r

    except Exception as e:
        return None


def classify_generalism(rna_si, dna_si, di):
    """classify the generalism mode based on SI values and DI."""
    if di < 0.10:
        # generalist — check which mode
        if rna_si < 0.5 and dna_si < 0.5:
            return "backbone_generalist"
        elif rna_si > 0.8 and dna_si > 0.8:
            return "base_generalist"
        else:
            return "generalist"
    elif di < 0.25:
        return "moderate"
    else:
        return "specialist"


def main():
    # step 1: load LUCA domains
    domains = load_luca_na_domains()
    print(f"loaded {len(domains)} LUCA+preLUCA NA-binding domains")

    # step 2: PDB census — find which have RNA and/or DNA co-crystals
    print(f"\n{'='*80}")
    print("STEP 1: PDB co-crystal census")
    print(f"{'='*80}")

    census = []
    for i, dom in enumerate(domains):
        pfam_id = dom['pfam_id']
        print(f"  [{i+1}/{len(domains)}] {pfam_id}...", end="", flush=True)

        n_rna, rna_pdbs = query_pdb_count(pfam_id, "RNA", max_results=3)
        time.sleep(0.3)
        n_dna, dna_pdbs = query_pdb_count(pfam_id, "DNA", max_results=3)
        time.sleep(0.3)

        has_both = n_rna > 0 and n_dna > 0
        tag = "BOTH" if has_both else ("RNA" if n_rna > 0 else ("DNA" if n_dna > 0 else "NONE"))
        print(f" RNA={n_rna} DNA={n_dna} → {tag}")

        dom['n_rna_pdb'] = n_rna
        dom['n_dna_pdb'] = n_dna
        dom['has_both'] = has_both
        dom['rna_pdbs'] = rna_pdbs
        dom['dna_pdbs'] = dna_pdbs
        census.append(dom)

    # save census
    census_path = os.path.join(BASE_DIR, "results", "luca_pdb_census.tsv")
    with open(census_path, 'w') as f:
        f.write("pfam_id\tage\trna_annotated\tdna_annotated\tn_rna_pdb\tn_dna_pdb\thas_both\n")
        for dom in census:
            f.write(f"{dom['pfam_id']}\t{dom['age']}\t{dom['rna_annotated']}\t{dom['dna_annotated']}\t{dom['n_rna_pdb']}\t{dom['n_dna_pdb']}\t{dom['has_both']}\n")
    print(f"\ncensus written to {census_path}")

    both_domains = [d for d in census if d['has_both']]
    print(f"\ndomains with BOTH RNA+DNA co-crystals: {len(both_domains)}/{len(domains)}")

    # step 3: for BOTH domains, find best structures and compute SI
    print(f"\n{'='*80}")
    print(f"STEP 2: NAS-Bench SI computation for {len(both_domains)} dual domains")
    print(f"{'='*80}")

    results = []
    for i, dom in enumerate(both_domains):
        pfam_id = dom['pfam_id']
        print(f"\n  [{i+1}/{len(both_domains)}] {pfam_id}")

        # find best RNA structure
        rna_pdb, rna_title, rna_total = find_best_structure(pfam_id, "RNA")
        if rna_pdb:
            print(f"    RNA: {rna_pdb} ({rna_title[:60]}...)" if len(rna_title or '') > 60 else f"    RNA: {rna_pdb} ({rna_title})")
        else:
            print(f"    RNA: no valid structure found")

        # find best DNA structure
        dna_pdb, dna_title, dna_total = find_best_structure(pfam_id, "DNA")
        if dna_pdb:
            print(f"    DNA: {dna_pdb} ({dna_title[:60]}...)" if len(dna_title or '') > 60 else f"    DNA: {dna_pdb} ({dna_title})")
        else:
            print(f"    DNA: no valid structure found")

        # compute SI for RNA
        rna_si_result = None
        if rna_pdb:
            rna_si_result = compute_si_safe(rna_pdb, "RNA")
            if rna_si_result:
                print(f"    RNA SI: {rna_si_result['SI']:.3f} (bb={rna_si_result['bb_frac']:.0%} base={rna_si_result['base_frac']:.0%} 2oh={rna_si_result['oh_frac']:.0%} n={rna_si_result['n_total']}) [{rna_si_result['_method']}]")
            else:
                print(f"    RNA SI: FAILED")
            gc.collect()

        # compute SI for DNA
        dna_si_result = None
        if dna_pdb:
            dna_si_result = compute_si_safe(dna_pdb, "DNA")
            if dna_si_result:
                print(f"    DNA SI: {dna_si_result['SI']:.3f} (bb={dna_si_result['bb_frac']:.0%} base={dna_si_result['base_frac']:.0%} 2oh={dna_si_result['oh_frac']:.0%} n={dna_si_result['n_total']}) [{dna_si_result['_method']}]")
            else:
                print(f"    DNA SI: FAILED")
            gc.collect()

        # if same PDB was selected for both, use split for the missing type
        if rna_pdb and dna_pdb and rna_pdb == dna_pdb:
            print(f"    NOTE: same PDB for both — using split analysis")
            # already handled by compute_si_safe with na_type_target

        # compute paired metrics
        row = {
            'pfam_id': pfam_id,
            'age': dom['age'],
            'rna_annotated': dom['rna_annotated'],
            'dna_annotated': dom['dna_annotated'],
            'n_rna_pdb': dom['n_rna_pdb'],
            'n_dna_pdb': dom['n_dna_pdb'],
            'rna_pdb': rna_pdb or '',
            'dna_pdb': dna_pdb or '',
        }

        if rna_si_result and dna_si_result:
            rsi = rna_si_result['SI']
            dsi = dna_si_result['SI']
            delta = rsi - dsi
            di = abs(delta) / (rsi + dsi) if (rsi + dsi) > 0 else 0
            mode = classify_generalism(rsi, dsi, di)

            row.update({
                'rna_SI': rsi,
                'rna_bb_frac': rna_si_result['bb_frac'],
                'rna_base_frac': rna_si_result['base_frac'],
                'rna_2oh_frac': rna_si_result['oh_frac'],
                'rna_contacts': rna_si_result['n_total'],
                'rna_method': rna_si_result['_method'],
                'dna_SI': dsi,
                'dna_bb_frac': dna_si_result['bb_frac'],
                'dna_base_frac': dna_si_result['base_frac'],
                'dna_2oh_frac': dna_si_result['oh_frac'],
                'dna_contacts': dna_si_result['n_total'],
                'dna_method': dna_si_result['_method'],
                'delta_SI': delta,
                'DI': di,
                'generalism_mode': mode,
            })
            print(f"    → ΔSI={delta:+.3f} DI={di:.3f} → {mode}")
        else:
            row.update({
                'rna_SI': rna_si_result['SI'] if rna_si_result else '',
                'rna_bb_frac': rna_si_result['bb_frac'] if rna_si_result else '',
                'rna_base_frac': rna_si_result['base_frac'] if rna_si_result else '',
                'rna_2oh_frac': rna_si_result['oh_frac'] if rna_si_result else '',
                'rna_contacts': rna_si_result['n_total'] if rna_si_result else '',
                'rna_method': rna_si_result.get('_method', '') if rna_si_result else '',
                'dna_SI': dna_si_result['SI'] if dna_si_result else '',
                'dna_bb_frac': dna_si_result['bb_frac'] if dna_si_result else '',
                'dna_base_frac': dna_si_result['base_frac'] if dna_si_result else '',
                'dna_2oh_frac': dna_si_result['oh_frac'] if dna_si_result else '',
                'dna_contacts': dna_si_result['n_total'] if dna_si_result else '',
                'dna_method': dna_si_result.get('_method', '') if dna_si_result else '',
                'delta_SI': '',
                'DI': '',
                'generalism_mode': 'incomplete',
            })
            print(f"    → incomplete (missing {'RNA' if not rna_si_result else 'DNA'} SI)")

        results.append(row)

    # step 4: output
    out_path = os.path.join(BASE_DIR, "results", "nasbench_full_luca.tsv")
    cols = ['pfam_id', 'age', 'rna_annotated', 'dna_annotated', 'n_rna_pdb', 'n_dna_pdb',
            'rna_pdb', 'rna_SI', 'rna_bb_frac', 'rna_base_frac', 'rna_2oh_frac', 'rna_contacts', 'rna_method',
            'dna_pdb', 'dna_SI', 'dna_bb_frac', 'dna_base_frac', 'dna_2oh_frac', 'dna_contacts', 'dna_method',
            'delta_SI', 'DI', 'generalism_mode']
    with open(out_path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for row in results:
            vals = []
            for c in cols:
                v = row.get(c, '')
                if isinstance(v, float):
                    vals.append(f"{v:.4f}")
                else:
                    vals.append(str(v))
            f.write('\t'.join(vals) + '\n')
    print(f"\nfull results written to {out_path}")

    # summary
    complete = [r for r in results if r['generalism_mode'] != 'incomplete']
    print(f"\n{'='*80}")
    print(f"SUMMARY: {len(complete)}/{len(results)} domains with complete paired SI")
    print(f"{'='*80}")

    if complete:
        generalists = [r for r in complete if 'generalist' in r['generalism_mode']]
        bb_gen = [r for r in complete if r['generalism_mode'] == 'backbone_generalist']
        base_gen = [r for r in complete if r['generalism_mode'] == 'base_generalist']
        moderate = [r for r in complete if r['generalism_mode'] == 'moderate']
        specialist = [r for r in complete if r['generalism_mode'] == 'specialist']

        print(f"\n  generalist (DI < 0.10):    {len(generalists)}/{len(complete)} ({100*len(generalists)/len(complete):.0f}%)")
        print(f"    backbone generalist:      {len(bb_gen)}")
        print(f"    base generalist:          {len(base_gen)}")
        print(f"    mixed generalist:         {len(generalists) - len(bb_gen) - len(base_gen)}")
        print(f"  moderate (DI 0.10-0.25):   {len(moderate)}/{len(complete)} ({100*len(moderate)/len(complete):.0f}%)")
        print(f"  specialist (DI > 0.25):    {len(specialist)}/{len(complete)} ({100*len(specialist)/len(complete):.0f}%)")

        dis = [r['DI'] for r in complete]
        print(f"\n  mean DI:   {np.mean(dis):.3f}")
        print(f"  median DI: {np.median(dis):.3f}")

        # table sorted by DI
        print(f"\n{'Pfam':<10} {'Age':<8} {'RNA_PDB':>7} {'RNA_SI':>7} {'DNA_PDB':>7} {'DNA_SI':>7} {'ΔSI':>7} {'DI':>6} {'Mode':<20}")
        print("-" * 95)
        for r in sorted(complete, key=lambda x: x['DI']):
            print(f"{r['pfam_id']:<10} {r['age']:<8} {r['rna_pdb']:>7} {r['rna_SI']:>7.3f} {r['dna_pdb']:>7} {r['dna_SI']:>7.3f} {r['delta_SI']:>+7.3f} {r['DI']:>6.3f} {r['generalism_mode']:<20}")


if __name__ == '__main__':
    main()
