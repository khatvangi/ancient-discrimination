#!/usr/bin/env python3
"""step 1-2: build a taxonomically diverse Alba domain (PF01918) dataset.

strategy:
  1. parse the Pfam seed alignment (76 sequences) — these are already domain-level
  2. extract UniProt accessions from the #=GS AC annotations
  3. batch-query UniProt for taxonomy of all 76 seed sequences
  4. classify into Archaea/Bacteria/Eukaryota
  5. ensure mandatory PDB structure proteins are included:
     - 3IAB (Pop6, S. cerevisiae — RNA-binding RNase P/MRP subunit)
     - 3U6Y (Ape10b2, Aeropyrum pernix — DNA-binding)
  6. CD-HIT at 95% identity (not 90%, Alba is small — preserve diversity)
  7. select 30-40 balanced across domains of life

usage:
  python asr_alba_step1_build_dataset.py
"""

import os
import sys
import re
import time
import json
import requests
import subprocess
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.join(SCRIPT_DIR, '..')
ALBA_DIR = os.path.join(PROJECT_DIR, 'results', 'asr', 'alba')
SEED_STO = os.path.join(ALBA_DIR, 'PF01918_seed.sto')

# output files
RAW_FASTA = os.path.join(ALBA_DIR, 'alba_raw_sequences.fasta')
CDHIT_FASTA = os.path.join(ALBA_DIR, 'alba_cdhit95.fasta')
FINAL_FASTA = os.path.join(ALBA_DIR, 'alba_sequences.fasta')
METADATA_TSV = os.path.join(ALBA_DIR, 'alba_sequence_metadata.tsv')
TAXONOMY_CACHE = os.path.join(ALBA_DIR, 'alba_taxonomy_cache.json')

# mandatory PDB proteins — must be in the final set
# 3IAB: Pop6 (S. cerevisiae RNase P/MRP subunit, RNA-binding)
#   Pop6 = POP6_YEAST = P53218 — already in seed alignment!
# 3U6Y: Ape10b2 = Alba2 from Aeropyrum pernix (DNA-binding)
#   ALBA2_AERPE = Q9YAX2 — already in seed alignment!
# also add:
# 1H0X/1H0Y: Sso10b (Sulfolobus solfataricus) = Alba1 — DNA-binding
#   Q97ZF4 (ALBA2_SACS2) is Sulfolobus solfataricus in seed, but may be a different paralog
# Rpp25 (RPP25_MOUSE, Q91WE3) — RNA-binding RNase P subunit, already in seed
MANDATORY_PDB = {
    'P53218': {'pdb': '3IAB', 'name': 'Pop6 (yeast)', 'na_type': 'RNA',
               'note': 'RNase P/MRP subunit'},
    'Q9YAX2': {'pdb': '3U6Y', 'name': 'Alba2 (A. pernix)', 'na_type': 'DNA',
               'note': 'crenarchaeal DNA-binding'},
    'Q91WE3': {'pdb': '3Q5H', 'name': 'Rpp25 (mouse)', 'na_type': 'RNA',
               'note': 'RNase P subunit'},
    'Q9YAW1': {'pdb': '2BKY', 'name': 'Alba1 (A. pernix)', 'na_type': 'DNA+RNA',
               'note': 'dual-binding Alba'},
}


def parseSeedAlignment(sto_path):
    """parse Stockholm alignment.
    returns:
      sequences: dict {seq_id: ungapped_sequence}
      accessions: dict {seq_id: uniprot_accession}
    """
    sequences = {}
    accessions = {}

    with open(sto_path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            # parse #=GS lines for accession mapping
            if line.startswith('#=GS'):
                # format: #=GS SEQ_ID  AC  ACCESSION.version
                match = re.match(r'#=GS\s+(\S+)\s+AC\s+(\S+)', line)
                if match:
                    seq_id = match.group(1)
                    acc_raw = match.group(2)
                    # remove version suffix (e.g. P38199.1 -> P38199)
                    acc = acc_raw.split('.')[0]
                    accessions[seq_id] = acc
                continue

            if line.startswith('#') or line.startswith('/'):
                continue

            # sequence lines: seq_id followed by aligned sequence
            parts = line.split()
            if len(parts) == 2:
                seq_id, aligned_seq = parts
                ungapped = aligned_seq.replace('.', '').replace('-', '')
                if seq_id in sequences:
                    sequences[seq_id] += ungapped
                else:
                    sequences[seq_id] = ungapped

    return sequences, accessions


def batchQueryTaxonomy(accession_list, batch_size=100):
    """query UniProt REST API in batches for taxonomy information.
    returns dict {accession: {organism, domain_of_life, lineage, phylum}}
    """
    taxonomy = {}

    for i in range(0, len(accession_list), batch_size):
        batch = accession_list[i:i+batch_size]
        # build query: accession:(P12345 OR P67890 OR ...)
        acc_query = ' OR '.join(batch)
        query = f'accession:({acc_query})'

        params = {
            'query': query,
            'format': 'json',
            'size': batch_size,
            'fields': 'accession,organism_name,organism_id,lineage',
        }

        try:
            resp = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params=params, timeout=60
            )
            if resp.status_code == 200:
                data = resp.json()
                for entry in data.get('results', []):
                    acc = entry.get('primaryAccession', '')
                    org = entry.get('organism', {}).get('scientificName', '')
                    org_id = entry.get('organism', {}).get('taxonId', 0)
                    lineage = entry.get('organism', {}).get('lineage', [])

                    domain_of_life = 'Unknown'
                    for taxon in lineage:
                        if taxon in ('Archaea', 'Bacteria', 'Eukaryota'):
                            domain_of_life = taxon
                            break

                    # extract phylum (first lineage element after domain)
                    phylum = ''
                    if lineage:
                        found_domain = False
                        for taxon in lineage:
                            if taxon in ('Archaea', 'Bacteria', 'Eukaryota'):
                                found_domain = True
                                continue
                            if found_domain:
                                phylum = taxon
                                break

                    taxonomy[acc] = {
                        'organism': org,
                        'organism_id': org_id,
                        'domain_of_life': domain_of_life,
                        'phylum': phylum,
                        'lineage': lineage,
                    }
            else:
                print(f"    warning: HTTP {resp.status_code} for batch starting at {i}")
        except Exception as e:
            print(f"    warning: API error at batch {i}: {e}")

        done = min(i + batch_size, len(accession_list))
        print(f"    queried {done}/{len(accession_list)} accessions, got {len(taxonomy)} results")
        time.sleep(1.0)

    return taxonomy


def main():
    os.makedirs(ALBA_DIR, exist_ok=True)

    if not os.path.exists(SEED_STO):
        print(f"ERROR: seed alignment not found at {SEED_STO}")
        sys.exit(1)

    # ================================================================
    print("=" * 70)
    print("STEP 1: parse Pfam seed alignment for Alba (PF01918)")
    print("=" * 70)
    sequences, accessions = parseSeedAlignment(SEED_STO)
    print(f"  parsed {len(sequences)} sequences")
    print(f"  accession mappings: {len(accessions)}")

    # get unique accessions
    unique_accs = list(set(accessions.values()))
    print(f"  unique UniProt accessions: {len(unique_accs)}")

    # length stats
    lengths = [len(s) for s in sequences.values()]
    print(f"  sequence lengths: min={min(lengths)}, max={max(lengths)}, "
          f"median={sorted(lengths)[len(lengths)//2]}")

    # check mandatory PDB proteins are in the seed
    seed_accs = set(accessions.values())
    print(f"\n  mandatory PDB protein check:")
    for acc, info in MANDATORY_PDB.items():
        in_seed = acc in seed_accs
        status = "IN SEED" if in_seed else "NOT in seed — will add"
        print(f"    {acc} ({info['name']}, PDB: {info['pdb']}): {status}")

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 2: query UniProt for taxonomy")
    print("=" * 70)

    if os.path.exists(TAXONOMY_CACHE):
        print("  loading cached taxonomy...")
        with open(TAXONOMY_CACHE) as f:
            taxonomy = json.load(f)
        print(f"  loaded {len(taxonomy)} cached entries")
        missing = [a for a in unique_accs if a not in taxonomy]
        if missing:
            print(f"  querying {len(missing)} missing accessions...")
            new_tax = batchQueryTaxonomy(missing)
            taxonomy.update(new_tax)
            with open(TAXONOMY_CACHE, 'w') as f:
                json.dump(taxonomy, f, indent=1)
    else:
        print(f"  querying taxonomy for {len(unique_accs)} accessions...")
        taxonomy = batchQueryTaxonomy(unique_accs)
        with open(TAXONOMY_CACHE, 'w') as f:
            json.dump(taxonomy, f, indent=1)
        print(f"  cached {len(taxonomy)} entries")

    # classify seed sequences by domain of life
    by_domain = defaultdict(list)
    by_phylum = defaultdict(lambda: defaultdict(list))
    no_taxonomy = []

    for seq_id, seq in sequences.items():
        acc = accessions.get(seq_id, '')
        if acc and acc in taxonomy:
            tax = taxonomy[acc]
            dol = tax['domain_of_life']
            phylum = tax.get('phylum', 'Unknown')
            by_domain[dol].append((seq_id, acc, seq, tax))
            by_phylum[dol][phylum].append((seq_id, acc, seq, tax))
        else:
            no_taxonomy.append((seq_id, acc, seq))

    print(f"\n  seed taxonomy breakdown:")
    for d in ['Archaea', 'Bacteria', 'Eukaryota', 'Unknown']:
        n = len(by_domain.get(d, []))
        if n > 0:
            print(f"    {d}: {n} sequences")
            phyla = by_phylum.get(d, {})
            for p in sorted(phyla.keys(), key=lambda x: len(phyla[x]), reverse=True):
                print(f"      {p}: {len(phyla[p])}")
    if no_taxonomy:
        print(f"    no taxonomy: {len(no_taxonomy)}")

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 3: build raw FASTA (seed sequences with annotated headers)")
    print("=" * 70)

    raw_sequences = {}

    for seq_id, seq in sequences.items():
        acc = accessions.get(seq_id, '')
        tax = taxonomy.get(acc, {})
        dol = tax.get('domain_of_life', 'Unknown')
        org = tax.get('organism', 'Unknown')
        org_clean = org.replace(' ', '_').replace('(', '').replace(')', '').replace(',', '').replace("'", '')[:40]

        # check if this is a mandatory PDB accession
        has_pdb = ''
        for mand_acc, mand_info in MANDATORY_PDB.items():
            if acc == mand_acc:
                has_pdb = mand_info['pdb']
                break

        new_id = f"{org_clean}|{acc}|{dol}|{has_pdb if has_pdb else 'no'}"
        raw_sequences[new_id] = {
            'sequence': seq,
            'accession': acc,
            'original_id': seq_id,
            'organism': org,
            'domain_of_life': dol,
            'has_pdb': has_pdb,
        }

    # write raw FASTA
    with open(RAW_FASTA, 'w') as f:
        for seq_id, meta in raw_sequences.items():
            f.write(f'>{seq_id}\n')
            seq = meta['sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

    print(f"  wrote {len(raw_sequences)} sequences to {RAW_FASTA}")

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 4: CD-HIT clustering at 95% identity")
    print("  (using 95% not 90% because Alba is a small family — preserve diversity)")
    print("=" * 70)

    cmd = [
        'conda', 'run', '-n', 'phylo_asr',
        'cd-hit',
        '-i', RAW_FASTA,
        '-o', CDHIT_FASTA,
        '-c', '0.95',
        '-n', '5',
        '-M', '4000',
        '-T', '8',
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  CD-HIT stderr: {result.stderr[:500]}")

    # count output
    cdhit_seqs = {}
    current_id = None
    current_seq = []
    with open(CDHIT_FASTA) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    cdhit_seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            cdhit_seqs[current_id] = ''.join(current_seq)

    print(f"  input: {len(raw_sequences)} sequences")
    print(f"  output: {len(cdhit_seqs)} representatives at 95% identity")

    # classify by domain of life
    cdhit_by_domain = defaultdict(list)
    cdhit_by_phylum = defaultdict(lambda: defaultdict(list))
    for seq_id in cdhit_seqs:
        parts = seq_id.split('|')
        dol = parts[2] if len(parts) >= 3 else 'Unknown'
        acc = parts[1] if len(parts) >= 2 else ''
        tax = taxonomy.get(acc, {})
        phylum = tax.get('phylum', 'Unknown')
        cdhit_by_domain[dol].append(seq_id)
        cdhit_by_phylum[dol][phylum].append(seq_id)

    print(f"\n  CD-HIT taxonomy breakdown:")
    for d in ['Archaea', 'Bacteria', 'Eukaryota', 'Unknown']:
        n = len(cdhit_by_domain.get(d, []))
        if n > 0:
            print(f"    {d}: {n}")
            phyla = cdhit_by_phylum.get(d, {})
            for p in sorted(phyla.keys(), key=lambda x: len(phyla[x]), reverse=True):
                print(f"      {p}: {len(phyla[p])}")

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 5: select balanced set (target 30-40 sequences)")
    print("=" * 70)

    # with 76 seed seqs after 95% CD-HIT, we expect most to survive
    # if total <= 45, keep them all — Alba is small enough
    # if too many, subsample with balance

    target_max = 45
    total_available = len(cdhit_seqs)

    if total_available <= target_max:
        print(f"  {total_available} <= {target_max}: using ALL CD-HIT representatives")
        selected = dict(cdhit_seqs)
    else:
        print(f"  {total_available} > {target_max}: subsampling with taxonomic balance")

        n_archaea = len(cdhit_by_domain.get('Archaea', []))
        n_eukaryota = len(cdhit_by_domain.get('Eukaryota', []))

        # for Alba: archaea are critical (the DNA-binding lineage),
        # keep all archaea, subsample eukaryota if needed
        alloc = {
            'Archaea': n_archaea,  # keep all archaea
            'Eukaryota': min(n_eukaryota, target_max - n_archaea),
        }
        # bacteria probably don't have Alba, but handle just in case
        n_bacteria = len(cdhit_by_domain.get('Bacteria', []))
        if n_bacteria > 0:
            alloc['Bacteria'] = min(n_bacteria, 5)
            alloc['Eukaryota'] = min(n_eukaryota, target_max - n_archaea - alloc['Bacteria'])

        print(f"  allocation: {dict(alloc)}")

        selected = {}

        # always include mandatory PDB proteins
        for seq_id, seq in cdhit_seqs.items():
            parts = seq_id.split('|')
            has_pdb = parts[3] if len(parts) >= 4 else 'no'
            if has_pdb != 'no':
                selected[seq_id] = seq
                dol = parts[2] if len(parts) >= 3 else 'Unknown'
                if dol in alloc:
                    alloc[dol] = max(0, alloc[dol] - 1)

        # for each domain, select with phylum diversity
        for domain, target_count in alloc.items():
            phyla = cdhit_by_phylum.get(domain, {})
            candidates = cdhit_by_domain.get(domain, [])
            candidates = [c for c in candidates if c not in selected]

            if len(candidates) <= target_count:
                for c in candidates:
                    selected[c] = cdhit_seqs[c]
            else:
                # round-robin across phyla for diversity
                added = 0
                phylum_lists = {p: list(ids) for p, ids in phyla.items()}
                phylum_order = sorted(phylum_lists.keys(),
                                      key=lambda p: len(phylum_lists[p]),
                                      reverse=True)

                while added < target_count:
                    progress = False
                    for p in phylum_order:
                        if added >= target_count:
                            break
                        ids = phylum_lists[p]
                        for cand_id in ids:
                            if cand_id not in selected:
                                selected[cand_id] = cdhit_seqs[cand_id]
                                ids.remove(cand_id)
                                added += 1
                                progress = True
                                break
                    if not progress:
                        break

    # verify mandatory proteins are included
    print(f"\n  mandatory PDB protein verification:")
    for acc, info in MANDATORY_PDB.items():
        found = any(f"|{acc}|" in sid for sid in selected)
        if found:
            print(f"    {acc} ({info['name']}): INCLUDED")
        else:
            print(f"    {acc} ({info['name']}): MISSING — adding from raw")
            for raw_id, raw_meta in raw_sequences.items():
                if raw_meta['accession'] == acc:
                    selected[raw_id] = raw_meta['sequence']
                    print(f"      added: {raw_id}")
                    break

    # final count
    final_by_domain = defaultdict(int)
    for seq_id in selected:
        parts = seq_id.split('|')
        dol = parts[2] if len(parts) >= 3 else 'Unknown'
        final_by_domain[dol] += 1

    print(f"\n  final set: {len(selected)} sequences")
    for d in ['Archaea', 'Bacteria', 'Eukaryota', 'Unknown']:
        if final_by_domain[d] > 0:
            print(f"    {d}: {final_by_domain[d]}")

    # write final FASTA
    with open(FINAL_FASTA, 'w') as f:
        for seq_id, seq in selected.items():
            f.write(f'>{seq_id}\n')
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

    print(f"\n  wrote {len(selected)} sequences to {FINAL_FASTA}")

    # write metadata TSV
    with open(METADATA_TSV, 'w') as f:
        f.write('seq_id\taccession\torganism\tdomain_of_life\tphylum\t'
                'seq_length\thas_pdb\tna_type\tnote\n')
        for seq_id in selected:
            parts = seq_id.split('|')
            acc = parts[1] if len(parts) >= 2 else ''
            dol = parts[2] if len(parts) >= 3 else 'Unknown'
            has_pdb = parts[3] if len(parts) >= 4 else 'no'
            tax = taxonomy.get(acc, {})
            org = tax.get('organism', parts[0].replace('_', ' ') if parts else 'Unknown')
            phylum = tax.get('phylum', 'Unknown')
            seq_len = len(selected[seq_id])

            # annotate known NA-binding type
            na_type = ''
            note = ''
            if acc in MANDATORY_PDB:
                na_type = MANDATORY_PDB[acc]['na_type']
                note = MANDATORY_PDB[acc]['note']

            f.write(f"{seq_id}\t{acc}\t{org}\t{dol}\t{phylum}\t"
                    f"{seq_len}\t{has_pdb}\t{na_type}\t{note}\n")

    print(f"  wrote metadata to {METADATA_TSV}")

    # ================================================================
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  seed alignment:    {len(sequences)} sequences")
    print(f"  after CD-HIT 95%: {len(cdhit_seqs)} representatives")
    print(f"  final selected:    {len(selected)} sequences")
    print(f"    Archaea:    {final_by_domain.get('Archaea', 0)}")
    print(f"    Bacteria:   {final_by_domain.get('Bacteria', 0)}")
    print(f"    Eukaryota:  {final_by_domain.get('Eukaryota', 0)}")
    print(f"    Unknown:    {final_by_domain.get('Unknown', 0)}")
    print()
    print("  key Alba biology:")
    print("    - archaea: ancestrally RNA-binding chromatin protein")
    print("    - crenarchaea (Sulfolobus, Aeropyrum): derived DNA-binding")
    print("    - eukaryota: Rpp20/Rpp25/Pop6/Pop7 = RNA-binding (RNase P/MRP)")
    print()
    print(f"  output files:")
    print(f"    {FINAL_FASTA}")
    print(f"    {METADATA_TSV}")
    print()
    print("  NEXT: run MAFFT alignment (step 3)")


if __name__ == '__main__':
    main()
