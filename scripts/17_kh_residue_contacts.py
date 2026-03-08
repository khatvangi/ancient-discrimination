#!/usr/bin/env python3
"""
KH domain per-residue contact analysis.

computes per-residue NA contact profiles for two KH domain co-crystal structures:
  5WWW — KH1 of MEX-3C + RNA (UniProt Q5U5Q3, residues 221-306)
  2P2R — KH3 of PCBP2 + DNA (UniProt Q15366, residues 285-359)

outputs:
  results/asr/kh/kh_rna_residue_contacts.tsv — per-residue contacts from 5WWW
  results/asr/kh/kh_dna_residue_contacts.tsv — per-residue contacts from 2P2R
  results/asr/kh/kh_specificity_map.tsv     — aligned comparison + classification
  results/asr/kh/kh_critical_columns.tsv    — PDB-to-UniProt mapped contact positions
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from nasbench import compute_residue_contacts, fetch, AA_3TO1
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner
import urllib.request
import csv

OUTDIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'asr', 'kh')

# ========================
# structures to analyze
# ========================
STRUCTURES = {
    '5WWW': {
        'desc': 'KH1 of MEX-3C + RNA',
        'na_type': 'RNA',
        'uniprot': 'Q5U5Q3',
        'prot_chain': 'A',
        'pdb_range': (221, 306),  # from DBREF
    },
    '2P2R': {
        'desc': 'KH3 of PCBP2 + DNA',
        'na_type': 'DNA',
        'uniprot': 'Q15366',
        'prot_chain': 'A',
        'pdb_range': (285, 359),  # from DBREF
    },
}


def extract_pdb_sequence(pdb_path, pdb_id, chain_id):
    """extract amino acid sequence from ATOM records of a PDB chain.

    returns list of (resnum, resname, aa1) tuples and the plain sequence string.
    """
    from Bio.PDB import PDBParser
    s = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    m = s[0]

    residues = []
    for r in m[chain_id]:
        rn = r.get_resname().strip().upper()
        if rn not in AA_3TO1:
            continue
        if r.get_id()[0] not in (' ', 'A'):
            continue
        aa1 = AA_3TO1[rn]
        residues.append((r.get_id()[1], rn, aa1))

    seq = ''.join(aa for _, _, aa in residues)
    return residues, seq


def fetchUniProtSequence(uniprot_id):
    """download full UniProt sequence in FASTA format."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = urllib.request.urlopen(url)
    lines = response.read().decode().strip().split('\n')
    # skip header line starting with '>'
    seq = ''.join(l.strip() for l in lines if not l.startswith('>'))
    return seq


def alignSequences(seq1, seq2):
    """align two sequences using Needleman-Wunsch (BioPython PairwiseAligner).

    returns the best alignment as two aligned strings.
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(seq1, seq2)
    best = alignments[0]

    # extract aligned strings
    aligned = str(best).split('\n')
    # biopython PairwiseAligner format: line 0 = seq1 aligned, line 2 = seq2 aligned
    aln_seq1 = aligned[0]
    aln_seq2 = aligned[2]

    return aln_seq1, aln_seq2


def buildPdbToUniprotMap(pdb_residues, pdb_seq, uniprot_seq):
    """map PDB residue numbers to UniProt positions via LOCAL sequence alignment.

    uses local alignment because PDB sequences are fragments of the full UniProt
    sequence, often with expression tags (His-tags etc) that aren't in UniProt.

    pdb_residues: list of (resnum, resname, aa1) from PDB
    pdb_seq: string sequence from PDB
    uniprot_seq: full UniProt sequence

    returns dict: pdb_resnum -> uniprot_position (1-based), or None for unmapped
    """
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    alignments = aligner.align(pdb_seq, uniprot_seq)
    best = alignments[0]

    # the alignment object has .aligned property:
    # best.aligned[0] = list of (start, end) intervals in pdb_seq
    # best.aligned[1] = list of (start, end) intervals in uniprot_seq
    pdb_blocks = best.aligned[0]
    uni_blocks = best.aligned[1]

    # build a mapping from pdb_seq index -> uniprot_seq index (both 0-based)
    seq_idx_map = {}
    for (ps, pe), (us, ue) in zip(pdb_blocks, uni_blocks):
        for offset in range(pe - ps):
            seq_idx_map[ps + offset] = us + offset

    # now map pdb residue numbers through: resnum -> pdb_seq index -> uniprot index -> uniprot position
    mapping = {}
    for idx, (resnum, resname, aa1) in enumerate(pdb_residues):
        if idx in seq_idx_map:
            # convert 0-based uniprot index to 1-based position
            mapping[resnum] = seq_idx_map[idx] + 1
        else:
            mapping[resnum] = None

    return mapping


def saveResidueContacts(residue_data, filepath):
    """save per-residue contact table as TSV."""
    fieldnames = ['chain', 'resnum', 'resname', 'aa1',
                  'n_bb', 'n_base', 'n_2oh', 'n_sr', 'n_total',
                  'dominant', 'na_type']

    with open(filepath, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t',
                           extrasaction='ignore')
        w.writeheader()
        for rd in residue_data:
            w.writerow(rd)

    print(f"  saved {len(residue_data)} residues to {filepath}")


# ========================================================================
# MAIN
# ========================================================================
def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # ---------------------------------------------------------------
    # step 3: compute per-residue contacts for both structures
    # ---------------------------------------------------------------
    print("=" * 60)
    print("step 3: compute per-residue contacts")
    print("=" * 60)

    all_residue_data = {}
    all_summaries = {}

    for pdb_id, info in STRUCTURES.items():
        print(f"\n--- {pdb_id}: {info['desc']} ---")
        path = fetch(pdb_id)
        residues, summary = compute_residue_contacts(path, pdb_id)

        all_residue_data[pdb_id] = residues
        all_summaries[pdb_id] = summary

        print(f"  NA type: {summary['na_type']}")
        print(f"  contacting residues: {summary['n_pres']}")
        print(f"  total contacts: {summary['n_total']}")
        print(f"  SI = {summary['SI']:.3f}")
        print(f"  bb={summary['n_bb']} base={summary['n_base']} 2oh={summary['n_2oh']} sr={summary['n_sr']}")

        # print 2'OH contacting residues if any
        oh_residues = [r for r in residues if r['n_2oh'] > 0]
        if oh_residues:
            print(f"  residues contacting 2'-OH ({len(oh_residues)}):")
            for r in oh_residues:
                print(f"    {r['chain']}:{r['resnum']} {r['resname']}({r['aa1']}) "
                      f"— 2oh={r['n_2oh']} bb={r['n_bb']} base={r['n_base']} sr={r['n_sr']}")

    # save TSVs
    rna_out = os.path.join(OUTDIR, 'kh_rna_residue_contacts.tsv')
    dna_out = os.path.join(OUTDIR, 'kh_dna_residue_contacts.tsv')
    saveResidueContacts(all_residue_data['5WWW'], rna_out)
    saveResidueContacts(all_residue_data['2P2R'], dna_out)

    # ---------------------------------------------------------------
    # step 4: align KH domain sequences and classify specificity
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("step 4: identify specificity residues")
    print("=" * 60)

    # extract PDB sequences
    pdb_seqs = {}
    pdb_res_lists = {}
    for pdb_id, info in STRUCTURES.items():
        path = fetch(pdb_id)
        res_list, seq = extract_pdb_sequence(path, pdb_id, info['prot_chain'])
        pdb_seqs[pdb_id] = seq
        pdb_res_lists[pdb_id] = res_list
        print(f"\n{pdb_id} chain {info['prot_chain']}: {len(res_list)} residues")
        print(f"  range: {res_list[0][0]}-{res_list[-1][0]}")
        print(f"  seq: {seq[:60]}...")

    # align the two KH domain sequences to find equivalent positions
    print("\naligning KH domains (5WWW vs 2P2R)...")
    rna_seq = pdb_seqs['5WWW']
    dna_seq = pdb_seqs['2P2R']
    aln_rna, aln_dna = alignSequences(rna_seq, dna_seq)

    print(f"  alignment length: {len(aln_rna)}")
    # show alignment
    for start in range(0, len(aln_rna), 60):
        end = min(start + 60, len(aln_rna))
        mid = ''
        for k in range(start, end):
            if aln_rna[k] == '-' or aln_dna[k] == '-':
                mid += ' '
            elif aln_rna[k] == aln_dna[k]:
                mid += '|'
            else:
                mid += '.'
        print(f"  5WWW: {aln_rna[start:end]}")
        print(f"        {mid}")
        print(f"  2P2R: {aln_dna[start:end]}")
        print()

    # build contact lookup dicts keyed by resnum
    rna_contact_map = {r['resnum']: r for r in all_residue_data['5WWW']}
    dna_contact_map = {r['resnum']: r for r in all_residue_data['2P2R']}

    # walk through alignment and classify each position
    rna_idx, dna_idx = 0, 0
    rna_residues = pdb_res_lists['5WWW']
    dna_residues = pdb_res_lists['2P2R']

    specificity_rows = []

    for col in range(len(aln_rna)):
        has_rna = aln_rna[col] != '-'
        has_dna = aln_dna[col] != '-'

        rna_resnum = rna_residues[rna_idx][0] if has_rna and rna_idx < len(rna_residues) else None
        rna_resname = rna_residues[rna_idx][1] if has_rna and rna_idx < len(rna_residues) else ''
        dna_resnum = dna_residues[dna_idx][0] if has_dna and dna_idx < len(dna_residues) else None
        dna_resname = dna_residues[dna_idx][1] if has_dna and dna_idx < len(dna_residues) else ''

        # get contact profiles if residue is in contact with NA
        rna_c = rna_contact_map.get(rna_resnum) if rna_resnum else None
        dna_c = dna_contact_map.get(dna_resnum) if dna_resnum else None

        # format contact profile strings
        def profileStr(c):
            if c is None:
                return "no_contact"
            return f"bb={c['n_bb']},base={c['n_base']},2oh={c['n_2oh']},sr={c['n_sr']}"

        # classify the position
        classification = "non_contact"
        if rna_c and dna_c:
            # both structures have contact at this aligned position
            if rna_c['n_2oh'] > 0:
                classification = "rna_specific"  # contacts 2'-OH in RNA structure
            elif rna_c['n_bb'] > 0 and dna_c['n_bb'] > 0:
                classification = "generalist"  # backbone contacts in both
            elif rna_c['n_base'] > 0 and dna_c['n_base'] > 0:
                classification = "generalist"  # base contacts in both
            else:
                classification = "generalist"  # contact in both, mixed categories
        elif rna_c and not dna_c:
            if rna_c['n_2oh'] > 0:
                classification = "rna_specific"
            else:
                classification = "rna_only"  # contact only in RNA structure
        elif dna_c and not rna_c:
            classification = "dna_only"  # contact only in DNA structure
        # else: non_contact (neither structure has contact here)

        # refine: if dna_c has base contacts but equivalent rna_c position has none,
        # that's potentially DNA-specific
        if classification == "dna_only" and dna_c and dna_c['n_base'] > 0:
            classification = "dna_specific"

        row = {
            'aln_col': col + 1,  # 1-based
            'rna_resnum': rna_resnum if rna_resnum else '',
            'rna_resname': rna_resname,
            'rna_contact_profile': profileStr(rna_c),
            'dna_resnum': dna_resnum if dna_resnum else '',
            'dna_resname': dna_resname,
            'dna_contact_profile': profileStr(dna_c),
            'classification': classification,
        }
        specificity_rows.append(row)

        if has_rna:
            rna_idx += 1
        if has_dna:
            dna_idx += 1

    # save specificity map
    spec_out = os.path.join(OUTDIR, 'kh_specificity_map.tsv')
    fieldnames = ['aln_col', 'rna_resnum', 'rna_resname', 'rna_contact_profile',
                  'dna_resnum', 'dna_resname', 'dna_contact_profile', 'classification']
    with open(spec_out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for row in specificity_rows:
            w.writerow(row)
    print(f"saved specificity map: {spec_out}")

    # summarize classifications
    from collections import Counter
    class_counts = Counter(r['classification'] for r in specificity_rows)
    print(f"\nposition classifications:")
    for cls, cnt in sorted(class_counts.items()):
        print(f"  {cls}: {cnt}")

    # highlight specificity residues
    print("\n--- RNA-specific residues (2'-OH contacts) ---")
    for row in specificity_rows:
        if row['classification'] == 'rna_specific':
            print(f"  aln_col={row['aln_col']} 5WWW:{row['rna_resnum']} {row['rna_resname']} "
                  f"profile={row['rna_contact_profile']}")

    print("\n--- DNA-specific residues ---")
    for row in specificity_rows:
        if row['classification'] == 'dna_specific':
            print(f"  aln_col={row['aln_col']} 2P2R:{row['dna_resnum']} {row['dna_resname']} "
                  f"profile={row['dna_contact_profile']}")

    print("\n--- generalist residues (contact in both) ---")
    for row in specificity_rows:
        if row['classification'] == 'generalist':
            print(f"  aln_col={row['aln_col']} "
                  f"5WWW:{row['rna_resnum']} {row['rna_resname']}({row['rna_contact_profile']}) | "
                  f"2P2R:{row['dna_resnum']} {row['dna_resname']}({row['dna_contact_profile']})")

    # ---------------------------------------------------------------
    # step 5: map PDB positions to UniProt
    # ---------------------------------------------------------------
    print("\n" + "=" * 60)
    print("step 5: map PDB contacts to UniProt positions")
    print("=" * 60)

    critical_rows = []

    for pdb_id, info in STRUCTURES.items():
        print(f"\n--- {pdb_id}: {info['desc']} ---")

        # fetch UniProt sequence
        uniprot_id = info['uniprot']
        print(f"  fetching UniProt {uniprot_id}...")
        uniprot_seq = fetchUniProtSequence(uniprot_id)
        print(f"  UniProt sequence length: {len(uniprot_seq)}")

        # align PDB sequence to UniProt
        pdb_res = pdb_res_lists[pdb_id]
        pdb_seq = pdb_seqs[pdb_id]

        print(f"  aligning PDB seq ({len(pdb_seq)} aa) to UniProt ({len(uniprot_seq)} aa)...")
        pdb_to_uniprot = buildPdbToUniprotMap(pdb_res, pdb_seq, uniprot_seq)

        # show mapping range
        mapped_positions = [v for v in pdb_to_uniprot.values() if v is not None]
        if mapped_positions:
            print(f"  mapped range: UniProt {min(mapped_positions)}-{max(mapped_positions)}")

        # for each contacting residue, determine its classification from specificity map
        contact_data = all_residue_data[pdb_id]
        # build lookup from specificity map: resnum -> classification
        resnum_to_class = {}
        for row in specificity_rows:
            key_field = 'rna_resnum' if pdb_id == '5WWW' else 'dna_resnum'
            if row[key_field]:
                resnum = int(row[key_field]) if row[key_field] else None
                if resnum is not None:
                    resnum_to_class[resnum] = row['classification']

        for rd in contact_data:
            pdb_resnum = rd['resnum']
            uni_pos = pdb_to_uniprot.get(pdb_resnum)

            # determine dominant contact type
            contact_type = rd['dominant']

            # get classification from specificity map
            classification = resnum_to_class.get(pdb_resnum, 'unknown')

            critical_rows.append({
                'pdb_source': pdb_id,
                'pdb_resnum': pdb_resnum,
                'pdb_resname': rd['resname'],
                'aa1': rd['aa1'],
                'uniprot_id': info['uniprot'],
                'uniprot_position': uni_pos if uni_pos else '',
                'n_bb': rd['n_bb'],
                'n_base': rd['n_base'],
                'n_2oh': rd['n_2oh'],
                'n_sr': rd['n_sr'],
                'contact_type': contact_type,
                'classification': classification,
            })

            if uni_pos:
                print(f"  {pdb_resnum} {rd['resname']}({rd['aa1']}) -> UniProt {uni_pos} "
                      f"[{contact_type}] {classification}")
            else:
                print(f"  {pdb_resnum} {rd['resname']}({rd['aa1']}) -> UniProt ??? "
                      f"[{contact_type}] {classification}")

    # save critical columns
    crit_out = os.path.join(OUTDIR, 'kh_critical_columns.tsv')
    fieldnames = ['pdb_source', 'pdb_resnum', 'pdb_resname', 'aa1',
                  'uniprot_id', 'uniprot_position',
                  'n_bb', 'n_base', 'n_2oh', 'n_sr',
                  'contact_type', 'classification']
    with open(crit_out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        w.writeheader()
        for row in critical_rows:
            w.writerow(row)
    print(f"\nsaved critical columns: {crit_out}")

    # final summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\n5WWW (RNA): {all_summaries['5WWW']['n_pres']} contacting residues, SI={all_summaries['5WWW']['SI']:.3f}")
    print(f"2P2R (DNA): {all_summaries['2P2R']['n_pres']} contacting residues, SI={all_summaries['2P2R']['SI']:.3f}")

    rna_spec = sum(1 for r in specificity_rows if r['classification'] == 'rna_specific')
    dna_spec = sum(1 for r in specificity_rows if r['classification'] == 'dna_specific')
    generalist = sum(1 for r in specificity_rows if r['classification'] == 'generalist')
    print(f"\nspecificity classification across aligned positions:")
    print(f"  rna_specific (2'-OH contacts): {rna_spec}")
    print(f"  dna_specific: {dna_spec}")
    print(f"  generalist: {generalist}")

    print(f"\noutput files:")
    print(f"  {rna_out}")
    print(f"  {dna_out}")
    print(f"  {spec_out}")
    print(f"  {crit_out}")


if __name__ == '__main__':
    main()
