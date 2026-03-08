#!/usr/bin/env python3
"""step 1-2 v2: build a taxonomically diverse KH domain (PF00013) dataset.

strategy:
  1. parse the Pfam seed alignment (735 sequences) — these are already domain-level
  2. extract UniProt accessions from the #=GS AC annotations
  3. batch-query UniProt for taxonomy of all 735 seed sequences
  4. classify into Archaea/Bacteria/Eukaryota
  5. add mandatory PDB structure proteins (extract KH domains)
  6. CD-HIT at 90% identity to remove redundancy
  7. select ~60-80 balanced across domains of life

usage:
  python asr_kh_step1v2_build_dataset.py
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
KH_DIR = os.path.join(PROJECT_DIR, 'results', 'asr', 'kh')
SEED_STO = os.path.join(KH_DIR, 'PF00013_seed.sto')

# output files
SEED_FASTA = os.path.join(KH_DIR, 'PF00013_seed.fasta')
RAW_FASTA = os.path.join(KH_DIR, 'kh_raw_sequences.fasta')
CDHIT_FASTA = os.path.join(KH_DIR, 'kh_cdhit90.fasta')
FINAL_FASTA = os.path.join(KH_DIR, 'kh_sequences.fasta')
METADATA_TSV = os.path.join(KH_DIR, 'kh_sequence_metadata.tsv')
TAXONOMY_CACHE = os.path.join(KH_DIR, 'seed_taxonomy_cache.json')

# mandatory PDB proteins — must be in the final set
MANDATORY_PDB = {
    'Q5U5Q3': {'pdb': '5WWW', 'name': 'MEX3C', 'na_type': 'RNA'},
    'Q15366': {'pdb': '2P2R', 'name': 'PCBP2', 'na_type': 'DNA'},
    'Q9UNW9': {'pdb': '1EC6', 'name': 'Nova-2', 'na_type': 'RNA'},
    'Q15365': {'pdb': '1ZTG', 'name': 'PCBP1', 'na_type': 'DNA'},
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

            # sequence lines
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
    returns dict {accession: {organism, domain_of_life, lineage}}
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

                    # extract phylum (second element in lineage after domain)
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

        # progress
        done = min(i + batch_size, len(accession_list))
        print(f"    queried {done}/{len(accession_list)} accessions, got {len(taxonomy)} results")
        time.sleep(1.0)  # be polite

    return taxonomy


def fetchProteinWithDomains(accession):
    """fetch a single protein from UniProt with its domain annotations."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            entry = resp.json()
            seq = entry.get('sequence', {}).get('value', '')
            org = entry.get('organism', {}).get('scientificName', '')
            lineage = entry.get('organism', {}).get('lineage', [])
            entry_name = entry.get('uniProtkbId', '')

            domain_of_life = 'Unknown'
            for taxon in lineage:
                if taxon in ('Archaea', 'Bacteria', 'Eukaryota'):
                    domain_of_life = taxon
                    break

            prot_name = ''
            if 'proteinDescription' in entry:
                pd = entry['proteinDescription']
                if 'recommendedName' in pd:
                    prot_name = pd['recommendedName'].get('fullName', {}).get('value', '')

            # extract KH domain positions
            kh_domains = []
            for feature in entry.get('features', []):
                if feature.get('type') == 'Domain':
                    desc = feature.get('description', '')
                    if 'KH' in desc.upper() or 'K HOMOLOGY' in desc.upper() or 'K-HOMOLOGY' in desc.upper():
                        loc = feature.get('location', {})
                        start = loc.get('start', {}).get('value', 0)
                        end = loc.get('end', {}).get('value', 0)
                        if start and end:
                            kh_domains.append((start, end, desc))

            return {
                'accession': accession,
                'entry_name': entry_name,
                'sequence': seq,
                'organism': org,
                'domain_of_life': domain_of_life,
                'protein_name': prot_name,
                'kh_domains': kh_domains,
            }
        else:
            return None
    except Exception as e:
        print(f"    warning: error fetching {accession}: {e}")
        return None


def main():
    os.makedirs(KH_DIR, exist_ok=True)

    if not os.path.exists(SEED_STO):
        print(f"ERROR: seed alignment not found at {SEED_STO}")
        sys.exit(1)

    # ================================================================
    print("=" * 70)
    print("STEP 1: parse Pfam seed alignment")
    print("=" * 70)
    sequences, accessions = parseSeedAlignment(SEED_STO)
    print(f"  parsed {len(sequences)} sequences")
    print(f"  accession mappings: {len(accessions)}")

    # get unique accessions
    unique_accs = list(set(accessions.values()))
    print(f"  unique UniProt accessions: {len(unique_accs)}")

    # length stats
    lengths = [len(s) for s in sequences.values()]
    print(f"  sequence lengths: min={min(lengths)}, max={max(lengths)}, median={sorted(lengths)[len(lengths)//2]}")

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 2: query UniProt for taxonomy of all seed sequences")
    print("=" * 70)

    # check for cached taxonomy
    if os.path.exists(TAXONOMY_CACHE):
        print("  loading cached taxonomy...")
        with open(TAXONOMY_CACHE) as f:
            taxonomy = json.load(f)
        print(f"  loaded {len(taxonomy)} cached entries")
        # check if we need to query more
        missing = [a for a in unique_accs if a not in taxonomy]
        if missing:
            print(f"  querying {len(missing)} missing accessions...")
            new_tax = batchQueryTaxonomy(missing)
            taxonomy.update(new_tax)
            with open(TAXONOMY_CACHE, 'w') as f:
                json.dump(taxonomy, f, indent=1)
    else:
        print(f"  querying taxonomy for {len(unique_accs)} accessions (this takes ~1 min)...")
        taxonomy = batchQueryTaxonomy(unique_accs)
        with open(TAXONOMY_CACHE, 'w') as f:
            json.dump(taxonomy, f, indent=1)
        print(f"  cached {len(taxonomy)} entries to {TAXONOMY_CACHE}")

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
            # show phyla
            phyla = by_phylum.get(d, {})
            for p in sorted(phyla.keys(), key=lambda x: len(phyla[x]), reverse=True)[:5]:
                print(f"      {p}: {len(phyla[p])}")
            if len(phyla) > 5:
                print(f"      ... and {len(phyla) - 5} more phyla")

    if no_taxonomy:
        print(f"    no taxonomy: {len(no_taxonomy)}")

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 3: fetch mandatory PDB proteins and extract KH domains")
    print("=" * 70)

    mandatory_kh_seqs = {}
    for acc, info in MANDATORY_PDB.items():
        print(f"  {acc} ({info['name']}, PDB: {info['pdb']})...", end=' ', flush=True)
        prot = fetchProteinWithDomains(acc)
        if prot:
            if prot['kh_domains']:
                # take the first KH domain that's in the right size range
                for start, end, desc in prot['kh_domains']:
                    domain_seq = prot['sequence'][start-1:end]
                    if 40 <= len(domain_seq) <= 100:
                        key = f"{prot['entry_name']}/{start}-{end}"
                        mandatory_kh_seqs[key] = {
                            'seq_id': key,
                            'accession': acc,
                            'sequence': domain_seq,
                            'organism': prot['organism'],
                            'domain_of_life': prot['domain_of_life'],
                            'protein_name': prot['protein_name'],
                            'has_pdb': info['pdb'],
                            'na_type': info['na_type'],
                        }
                        print(f"KH {start}-{end} ({len(domain_seq)} aa)")
                        break
                else:
                    print(f"no suitable KH domain found in {len(prot['kh_domains'])} annotations")
            else:
                print("no KH domain annotations")
        else:
            print("FAILED to fetch")
        time.sleep(0.3)

    print(f"  mandatory KH domain sequences: {len(mandatory_kh_seqs)}")

    # check if mandatory proteins are already in seed (by accession)
    seed_accs = set(accessions.values())
    for acc in MANDATORY_PDB:
        if acc in seed_accs:
            print(f"  note: {acc} already in seed alignment")

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 4: build raw FASTA (seed + mandatory)")
    print("=" * 70)

    # combine seed sequences with mandatory proteins
    raw_sequences = {}

    # add seed sequences with taxonomy-annotated headers
    for seq_id, seq in sequences.items():
        acc = accessions.get(seq_id, '')
        tax = taxonomy.get(acc, {})
        dol = tax.get('domain_of_life', 'Unknown')
        org = tax.get('organism', 'Unknown')
        org_clean = org.replace(' ', '_').replace('(', '').replace(')', '').replace(',', '')[:30]

        # check if this is one of the mandatory PDB accessions
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

    # add mandatory KH domain seqs if their accession is NOT already in seed
    for key, meta in mandatory_kh_seqs.items():
        if meta['accession'] not in seed_accs:
            org_clean = meta['organism'].replace(' ', '_').replace('(', '').replace(')', '')[:30]
            new_id = f"{org_clean}|{meta['accession']}|{meta['domain_of_life']}|{meta['has_pdb']}"
            raw_sequences[new_id] = {
                'sequence': meta['sequence'],
                'accession': meta['accession'],
                'original_id': key,
                'organism': meta['organism'],
                'domain_of_life': meta['domain_of_life'],
                'has_pdb': meta['has_pdb'],
            }
            print(f"  added mandatory: {meta['accession']} ({meta['protein_name']}, PDB: {meta['has_pdb']})")

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
    print("STEP 5: CD-HIT clustering at 90% identity")
    print("=" * 70)

    cmd = [
        'conda', 'run', '-n', 'phylo_asr',
        'cd-hit',
        '-i', RAW_FASTA,
        '-o', CDHIT_FASTA,
        '-c', '0.9',
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
    print(f"  output: {len(cdhit_seqs)} representatives at 90% identity")

    # classify CD-HIT representatives
    cdhit_by_domain = defaultdict(list)
    cdhit_by_phylum = defaultdict(lambda: defaultdict(list))
    for seq_id in cdhit_seqs:
        # parse domain of life from header
        parts = seq_id.split('|')
        if len(parts) >= 3:
            dol = parts[2]
        else:
            dol = 'Unknown'
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

    # ================================================================
    print()
    print("=" * 70)
    print("STEP 6: select balanced set (~70 sequences)")
    print("=" * 70)

    # strategy: take all CD-HIT reps (they're already non-redundant at 90%)
    # if too many, subsample maintaining domain of life balance
    # if too few, we're done

    target_n = 70
    total_available = len(cdhit_seqs)

    if total_available <= target_n:
        print(f"  {total_available} <= {target_n}: using ALL CD-HIT representatives")
        selected = cdhit_seqs
    else:
        print(f"  {total_available} > {target_n}: subsampling with taxonomic balance")

        # allocate slots: try ~20 each for Bac/Arc/Euk, fill rest proportionally
        n_archaea = len(cdhit_by_domain.get('Archaea', []))
        n_bacteria = len(cdhit_by_domain.get('Bacteria', []))
        n_eukaryota = len(cdhit_by_domain.get('Eukaryota', []))

        # proportional allocation with minimum 15 per domain if available
        alloc = {}
        for d, n_avail in [('Archaea', n_archaea), ('Bacteria', n_bacteria), ('Eukaryota', n_eukaryota)]:
            alloc[d] = min(n_avail, max(15, int(target_n * n_avail / total_available)))

        # adjust to hit target
        total_alloc = sum(alloc.values())
        if total_alloc < target_n:
            surplus = target_n - total_alloc
            for d in sorted(alloc.keys(), key=lambda x: len(cdhit_by_domain.get(x, [])), reverse=True):
                can_add = len(cdhit_by_domain.get(d, [])) - alloc[d]
                add_now = min(surplus, can_add)
                alloc[d] += add_now
                surplus -= add_now
                if surplus <= 0:
                    break

        print(f"  allocation: {dict(alloc)}")

        selected = {}

        # always include mandatory PDB proteins
        mandatory_accs = set(MANDATORY_PDB.keys())
        for seq_id, seq in cdhit_seqs.items():
            parts = seq_id.split('|')
            acc = parts[1] if len(parts) >= 2 else ''
            has_pdb = parts[3] if len(parts) >= 4 else 'no'
            if has_pdb != 'no':
                selected[seq_id] = seq
                # deduct from the domain's allocation
                dol = parts[2] if len(parts) >= 3 else 'Unknown'
                if dol in alloc:
                    alloc[dol] = max(0, alloc[dol] - 1)

        # for each domain, select with phylum diversity
        for domain, target_count in alloc.items():
            phyla = cdhit_by_phylum.get(domain, {})
            candidates = cdhit_by_domain.get(domain, [])

            # remove already-selected
            candidates = [c for c in candidates if c not in selected]

            if len(candidates) <= target_count:
                for c in candidates:
                    selected[c] = cdhit_seqs[c]
            else:
                # maximize phylum diversity: round-robin across phyla
                added = 0
                phylum_lists = {p: list(ids) for p, ids in phyla.items()}
                phylum_order = sorted(phylum_lists.keys(), key=lambda p: len(phylum_lists[p]), reverse=True)

                while added < target_count:
                    progress = False
                    for p in phylum_order:
                        if added >= target_count:
                            break
                        ids = phylum_lists[p]
                        # find one not yet selected
                        for cand_id in ids:
                            if cand_id not in selected:
                                selected[cand_id] = cdhit_seqs[cand_id]
                                ids.remove(cand_id)
                                added += 1
                                progress = True
                                break
                    if not progress:
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

    # check mandatory proteins are included
    for acc, info in MANDATORY_PDB.items():
        found = any(f"|{acc}|" in sid for sid in selected)
        if found:
            print(f"    mandatory {acc} ({info['name']}): INCLUDED")
        else:
            print(f"    mandatory {acc} ({info['name']}): MISSING (will try to add)")
            # find it in raw_sequences and add
            for raw_id, raw_meta in raw_sequences.items():
                if raw_meta['accession'] == acc:
                    selected[raw_id] = raw_meta['sequence']
                    print(f"      added from raw: {raw_id}")
                    break

    # write final FASTA
    with open(FINAL_FASTA, 'w') as f:
        for seq_id, seq in selected.items():
            f.write(f'>{seq_id}\n')
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

    print(f"\n  wrote {len(selected)} sequences to {FINAL_FASTA}")

    # write metadata TSV
    with open(METADATA_TSV, 'w') as f:
        f.write('seq_id\taccession\torganism\tdomain_of_life\tphylum\tseq_length\thas_pdb\n')
        for seq_id in selected:
            parts = seq_id.split('|')
            acc = parts[1] if len(parts) >= 2 else ''
            dol = parts[2] if len(parts) >= 3 else 'Unknown'
            has_pdb = parts[3] if len(parts) >= 4 else 'no'
            tax = taxonomy.get(acc, {})
            org = tax.get('organism', parts[0].replace('_', ' ') if parts else 'Unknown')
            phylum = tax.get('phylum', 'Unknown')
            seq_len = len(selected[seq_id])
            f.write(f"{seq_id}\t{acc}\t{org}\t{dol}\t{phylum}\t{seq_len}\t{has_pdb}\n")

    print(f"  wrote metadata to {METADATA_TSV}")

    # ================================================================
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  seed alignment:   {len(sequences)} sequences")
    print(f"  after CD-HIT 90%: {len(cdhit_seqs)} representatives")
    print(f"  final selected:   {len(selected)} sequences")
    print(f"    Archaea:    {final_by_domain.get('Archaea', 0)}")
    print(f"    Bacteria:   {final_by_domain.get('Bacteria', 0)}")
    print(f"    Eukaryota:  {final_by_domain.get('Eukaryota', 0)}")
    print(f"    Unknown:    {final_by_domain.get('Unknown', 0)}")
    print(f"\n  output files:")
    print(f"    {FINAL_FASTA}")
    print(f"    {METADATA_TSV}")


if __name__ == '__main__':
    main()
