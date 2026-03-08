#!/usr/bin/env python3
"""step 1-2: build a taxonomically diverse KH domain (PF00013) dataset for phylogenetic analysis.

converts the Pfam seed alignment (Stockholm) to FASTA,
queries UniProt for additional sequences ensuring taxonomic diversity,
runs CD-HIT to cluster at 90% identity,
and produces a curated set of ~60-80 sequences.

mandatory PDB structure organisms included:
  - 5WWW: MEX3C (Q5U5Q3, Homo sapiens) — RNA-bound KH
  - 2P2R: PCBP2 (Q15366, Homo sapiens) — DNA-bound KH
  - 1EC6: Nova-2 (Q9UNW9, Homo sapiens) — RNA specialist
  - 1ZTG: PCBP1 (Q15365, Homo sapiens) — DNA-bound KH

usage:
  python asr_kh_step1_build_dataset.py
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

# mandatory PDB proteins — these must be in the final set
# these are full-length proteins; we'll extract KH domains from them
MANDATORY_PDB = {
    'Q5U5Q3': {'pdb': '5WWW', 'name': 'MEX3C', 'na_type': 'RNA', 'organism': 'Homo sapiens'},
    'Q15366': {'pdb': '2P2R', 'name': 'PCBP2', 'na_type': 'DNA', 'organism': 'Homo sapiens'},
    'Q9UNW9': {'pdb': '1EC6', 'name': 'Nova-2', 'na_type': 'RNA', 'organism': 'Homo sapiens'},
    'Q15365': {'pdb': '1ZTG', 'name': 'PCBP1', 'na_type': 'DNA', 'organism': 'Homo sapiens'},
}

# taxonomic groups to sample from, with approximate desired counts
# format: (query_term, label, domain_of_life, desired_count)
TAXA_QUERIES = [
    # archaea — maximize phylum diversity
    ('taxonomy_id:2157 AND taxonomy_name:euryarchaeota', 'Euryarchaeota', 'Archaea', 5),
    ('taxonomy_id:2157 AND taxonomy_name:crenarchaeota', 'Crenarchaeota', 'Archaea', 3),
    ('taxonomy_id:2157 AND taxonomy_name:asgard', 'Asgardarchaeota', 'Archaea', 3),
    ('taxonomy_id:2157 AND NOT taxonomy_name:euryarchaeota AND NOT taxonomy_name:crenarchaeota AND NOT taxonomy_name:asgard', 'Other_Archaea', 'Archaea', 4),
    # bacteria — maximize phylum diversity
    ('taxonomy_id:1224', 'Proteobacteria', 'Bacteria', 4),
    ('taxonomy_id:1239', 'Firmicutes', 'Bacteria', 3),
    ('taxonomy_id:1760', 'Actinobacteria', 'Bacteria', 2),
    ('taxonomy_id:1117', 'Cyanobacteria', 'Bacteria', 2),
    ('taxonomy_id:200783', 'Aquificae', 'Bacteria', 2),
    ('taxonomy_id:188787', 'Thermotogae', 'Bacteria', 2),
    ('taxonomy_id:976', 'Bacteroidetes', 'Bacteria', 2),
    ('taxonomy_id:203691', 'Spirochaetes', 'Bacteria', 1),
    ('taxonomy_id:1297', 'Deinococcus-Thermus', 'Bacteria', 1),
    # eukaryota — include known specialists
    ('taxonomy_id:7742', 'Vertebrata', 'Eukaryota', 6),
    ('taxonomy_id:6656', 'Arthropoda', 'Eukaryota', 3),
    ('taxonomy_id:4751', 'Fungi', 'Eukaryota', 4),
    ('taxonomy_id:3193', 'Embryophyta', 'Eukaryota', 3),
    ('taxonomy_id:2759 AND NOT taxonomy_id:7742 AND NOT taxonomy_id:6656 AND NOT taxonomy_id:4751 AND NOT taxonomy_id:3193', 'Other_Eukaryota', 'Eukaryota', 3),
]

# well-known KH proteins to specifically search for
NAMED_PROTEINS = [
    ('Nova-1', 'NOVA1'),
    ('Nova-2', 'NOVA2'),
    ('FMRP', 'FMR1'),
    ('hnRNP K', 'HNRNPK'),
    ('PCBP2', 'PCBP2'),
    ('Vigilin', 'HDLBP'),
    ('PNPase', 'PNPT1'),
]


def parseStoToFasta(sto_path, fasta_path):
    """convert Stockholm alignment to FASTA (ungapped sequences).
    returns dict of {seq_id: sequence}
    """
    sequences = {}
    with open(sto_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('/'):
                continue
            parts = line.split()
            if len(parts) == 2:
                seq_id, aligned_seq = parts
                # remove gaps
                ungapped = aligned_seq.replace('.', '').replace('-', '')
                if seq_id in sequences:
                    sequences[seq_id] += ungapped
                else:
                    sequences[seq_id] = ungapped

    # write FASTA
    with open(fasta_path, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f'>{seq_id}\n')
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

    print(f"converted {len(sequences)} sequences from Stockholm to FASTA")
    print(f"  output: {fasta_path}")

    # report length stats
    lengths = [len(s) for s in sequences.values()]
    print(f"  length range: {min(lengths)}-{max(lengths)} aa, median: {sorted(lengths)[len(lengths)//2]} aa")

    return sequences


def queryUniprotForKH(query_filter, max_results=10, reviewed_only=True):
    """query UniProt for PF00013-containing proteins matching a taxonomic filter.
    returns list of dicts with accession, name, organism, taxonomy, sequence, domain_positions.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    # build query: PF00013 domain + taxonomic filter
    query_parts = ['(xref:pfam-PF00013)']
    if query_filter:
        query_parts.append(f'({query_filter})')
    if reviewed_only:
        query_parts.append('(reviewed:true)')

    query = ' AND '.join(query_parts)

    params = {
        'query': query,
        'format': 'json',
        'size': max_results,
        'fields': 'accession,id,organism_name,organism_id,lineage,sequence,protein_name,ft_domain',
    }

    try:
        resp = requests.get(base_url, params=params, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            results = []
            for entry in data.get('results', []):
                acc = entry.get('primaryAccession', '')
                entry_name = entry.get('uniProtkbId', '')
                org = entry.get('organism', {}).get('scientificName', '')
                org_id = entry.get('organism', {}).get('taxonId', 0)

                # get lineage for domain of life classification
                lineage = entry.get('organism', {}).get('lineage', [])
                domain_of_life = 'Unknown'
                for taxon in lineage:
                    if taxon in ('Archaea', 'Bacteria', 'Eukaryota'):
                        domain_of_life = taxon
                        break

                seq_data = entry.get('sequence', {})
                seq = seq_data.get('value', '')
                seq_len = seq_data.get('length', 0)

                # get protein name
                prot_name = ''
                if 'proteinDescription' in entry:
                    pd = entry['proteinDescription']
                    if 'recommendedName' in pd:
                        prot_name = pd['recommendedName'].get('fullName', {}).get('value', '')
                    elif 'submissionNames' in pd and pd['submissionNames']:
                        prot_name = pd['submissionNames'][0].get('fullName', {}).get('value', '')

                # extract KH domain positions from features
                kh_domains = []
                for feature in entry.get('features', []):
                    if feature.get('type') == 'Domain':
                        desc = feature.get('description', '')
                        if 'KH' in desc.upper() or 'K HOMOLOGY' in desc.upper():
                            loc = feature.get('location', {})
                            start = loc.get('start', {}).get('value', 0)
                            end = loc.get('end', {}).get('value', 0)
                            if start and end:
                                kh_domains.append((start, end, desc))

                results.append({
                    'accession': acc,
                    'entry_name': entry_name,
                    'organism': org,
                    'organism_id': org_id,
                    'domain_of_life': domain_of_life,
                    'lineage': lineage,
                    'sequence': seq,
                    'seq_length': seq_len,
                    'protein_name': prot_name,
                    'kh_domains': kh_domains,
                })
            return results
        else:
            print(f"  warning: HTTP {resp.status_code} for query: {query_filter}")
            return []
    except Exception as e:
        print(f"  warning: API error for {query_filter}: {e}")
        return []


def extractKHDomain(protein, min_len=50, max_len=100):
    """extract KH domain sequence from a full-length protein.
    if domain boundaries are annotated, use them.
    otherwise return None (will use full sequence through alignment later).
    """
    if protein['kh_domains']:
        # take the first KH domain
        start, end, desc = protein['kh_domains'][0]
        domain_seq = protein['sequence'][start-1:end]  # 1-based to 0-based
        if min_len <= len(domain_seq) <= max_len:
            return domain_seq, start, end
        # if annotated domain is too short/long, try all of them
        for s, e, d in protein['kh_domains']:
            ds = protein['sequence'][s-1:e]
            if min_len <= len(ds) <= max_len:
                return ds, s, e

    return None, None, None


def fetchMandatoryProteins():
    """fetch the mandatory PDB-structure proteins from UniProt."""
    mandatory_seqs = {}
    for acc, info in MANDATORY_PDB.items():
        print(f"  fetching mandatory {acc} ({info['name']})...", end=' ', flush=True)
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                entry = resp.json()
                seq = entry.get('sequence', {}).get('value', '')
                org = entry.get('organism', {}).get('scientificName', '')
                lineage = entry.get('organism', {}).get('lineage', [])
                domain_of_life = 'Unknown'
                for taxon in lineage:
                    if taxon in ('Archaea', 'Bacteria', 'Eukaryota'):
                        domain_of_life = taxon
                        break

                # extract KH domains
                kh_domains = []
                for feature in entry.get('features', []):
                    if feature.get('type') == 'Domain':
                        desc = feature.get('description', '')
                        if 'KH' in desc.upper() or 'K HOMOLOGY' in desc.upper():
                            loc = feature.get('location', {})
                            start = loc.get('start', {}).get('value', 0)
                            end = loc.get('end', {}).get('value', 0)
                            if start and end:
                                kh_domains.append((start, end, desc))

                prot_name = ''
                if 'proteinDescription' in entry:
                    pd = entry['proteinDescription']
                    if 'recommendedName' in pd:
                        prot_name = pd['recommendedName'].get('fullName', {}).get('value', '')

                mandatory_seqs[acc] = {
                    'accession': acc,
                    'entry_name': entry.get('uniProtkbId', ''),
                    'organism': org,
                    'organism_id': entry.get('organism', {}).get('taxonId', 0),
                    'domain_of_life': domain_of_life,
                    'lineage': lineage,
                    'sequence': seq,
                    'seq_length': len(seq),
                    'protein_name': prot_name,
                    'kh_domains': kh_domains,
                    'has_pdb': info['pdb'],
                    'na_type': info['na_type'],
                }
                print(f"OK ({len(seq)} aa, {len(kh_domains)} KH domains annotated)")
            else:
                print(f"FAILED (HTTP {resp.status_code})")
        except Exception as e:
            print(f"FAILED ({e})")
        time.sleep(0.3)

    return mandatory_seqs


def buildDiverseDataset():
    """build a taxonomically diverse dataset of KH domain sequences."""

    print("=" * 70)
    print("STEP 1: convert Pfam seed alignment")
    print("=" * 70)
    seed_seqs = parseStoToFasta(SEED_STO, SEED_FASTA)

    print()
    print("=" * 70)
    print("STEP 2: fetch mandatory PDB proteins")
    print("=" * 70)
    mandatory = fetchMandatoryProteins()

    print()
    print("=" * 70)
    print("STEP 3: query UniProt for taxonomically diverse KH proteins")
    print("=" * 70)

    all_proteins = {}  # accession -> protein dict

    # add mandatory proteins first
    for acc, prot in mandatory.items():
        all_proteins[acc] = prot

    # query each taxonomic group
    total_by_domain = defaultdict(int)
    for query_filter, label, domain, desired_n in TAXA_QUERIES:
        print(f"\n  [{domain}] {label} (want {desired_n})...", end=' ', flush=True)

        results = queryUniprotForKH(query_filter, max_results=desired_n * 2, reviewed_only=True)

        if not results:
            # fallback: try without reviewed restriction
            results = queryUniprotForKH(query_filter, max_results=desired_n * 2, reviewed_only=False)
            if results:
                print(f"(unreviewed) ", end='')

        if not results:
            print("NONE FOUND")
            continue

        # pick up to desired_n, preferring diverse organisms
        added = 0
        seen_genera = set()
        for prot in results:
            if added >= desired_n:
                break
            acc = prot['accession']
            if acc in all_proteins:
                continue

            # try to get genus for diversity
            genus = prot['organism'].split()[0] if prot['organism'] else 'unknown'
            if genus in seen_genera and added > 0:
                continue  # skip duplicate genera
            seen_genera.add(genus)

            prot['source_group'] = label
            all_proteins[acc] = prot
            total_by_domain[domain] += 1
            added += 1

        print(f"got {added}")
        time.sleep(0.5)  # be polite to API

    print(f"\n  total proteins collected: {len(all_proteins)}")
    for d in ['Archaea', 'Bacteria', 'Eukaryota']:
        print(f"    {d}: {total_by_domain[d]}")

    # now extract KH domain sequences
    # for proteins with annotated KH domains, extract the first one
    # for proteins without, we'll use the full sequence and let MAFFT handle it
    # actually — for phylogenetic analysis of the KH domain specifically,
    # we need domain-level sequences. let's extract just the KH domain.

    print()
    print("=" * 70)
    print("STEP 4: extract KH domain regions")
    print("=" * 70)

    domain_sequences = {}  # key -> (seq, metadata)
    no_annotation = []

    for acc, prot in all_proteins.items():
        kh_doms = prot.get('kh_domains', [])
        has_pdb = prot.get('has_pdb', '')

        if kh_doms:
            # extract each KH domain (proteins can have multiple)
            # for multi-KH proteins, take the first one (the one most likely
            # in the PDB structure) unless this is a mandatory PDB protein
            for idx, (start, end, desc) in enumerate(kh_doms):
                domain_seq = prot['sequence'][start-1:end]
                if len(domain_seq) < 40 or len(domain_seq) > 120:
                    continue  # skip oddly sized annotations

                key = f"{acc}_KH{idx+1}"
                org_short = prot['organism'].replace(' ', '_')[:20]
                domain_sequences[key] = {
                    'sequence': domain_seq,
                    'accession': acc,
                    'kh_index': idx + 1,
                    'domain_start': start,
                    'domain_end': end,
                    'organism': prot['organism'],
                    'domain_of_life': prot.get('domain_of_life', 'Unknown'),
                    'protein_name': prot.get('protein_name', ''),
                    'entry_name': prot.get('entry_name', ''),
                    'has_pdb': has_pdb,
                    'source_group': prot.get('source_group', 'mandatory'),
                }

                # for non-mandatory, just take the first KH domain
                if not has_pdb:
                    break
        else:
            # no KH annotation — will need to handle via hmmsearch or alignment
            no_annotation.append(acc)

    print(f"  extracted {len(domain_sequences)} KH domain sequences from annotated proteins")
    print(f"  {len(no_annotation)} proteins lack KH domain annotation")

    # for proteins without KH annotation, we'll do a simple approach:
    # use the full sequence and let MAFFT trim to the domain region.
    # but for a cleaner tree, let's search for the KH domain using hmmer on the seed.
    # for now, let's include the full sequence if it's short enough (<200 aa),
    # which would be a single-domain protein
    for acc in no_annotation:
        prot = all_proteins[acc]
        if prot['seq_length'] <= 200:
            key = f"{acc}_full"
            domain_sequences[key] = {
                'sequence': prot['sequence'],
                'accession': acc,
                'kh_index': 0,
                'domain_start': 1,
                'domain_end': prot['seq_length'],
                'organism': prot['organism'],
                'domain_of_life': prot.get('domain_of_life', 'Unknown'),
                'protein_name': prot.get('protein_name', ''),
                'entry_name': prot.get('entry_name', ''),
                'has_pdb': prot.get('has_pdb', ''),
                'source_group': prot.get('source_group', ''),
            }
        else:
            print(f"  skipping {acc} ({prot['protein_name']}, {prot['seq_length']} aa) — no KH annotation and too long")

    print(f"  total domain sequences for alignment: {len(domain_sequences)}")

    # write raw FASTA
    print()
    print("=" * 70)
    print("STEP 5: write raw FASTA")
    print("=" * 70)

    with open(RAW_FASTA, 'w') as f:
        for key, meta in domain_sequences.items():
            # header: >ORGANISM|ACCESSION_KHn|DOMAIN_OF_LIFE|HAS_PDB
            org_clean = meta['organism'].replace(' ', '_').replace('(', '').replace(')', '')
            has_pdb_str = meta['has_pdb'] if meta['has_pdb'] else 'no'
            header = f">{org_clean}|{meta['accession']}_KH{meta['kh_index']}|{meta['domain_of_life']}|{has_pdb_str}"
            f.write(header + '\n')
            seq = meta['sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

    print(f"  wrote {len(domain_sequences)} sequences to {RAW_FASTA}")

    # write metadata TSV
    with open(METADATA_TSV, 'w') as f:
        f.write('key\taccession\tkh_index\torganism\tdomain_of_life\tprotein_name\t'
                'entry_name\tdomain_start\tdomain_end\tseq_length\thas_pdb\tsource_group\n')
        for key, meta in domain_sequences.items():
            f.write(f"{key}\t{meta['accession']}\t{meta['kh_index']}\t{meta['organism']}\t"
                    f"{meta['domain_of_life']}\t{meta['protein_name']}\t{meta['entry_name']}\t"
                    f"{meta['domain_start']}\t{meta['domain_end']}\t{len(meta['sequence'])}\t"
                    f"{meta['has_pdb']}\t{meta['source_group']}\n")

    print(f"  wrote metadata to {METADATA_TSV}")

    return domain_sequences


def runCDHit(input_fasta, output_fasta, identity=0.9):
    """run CD-HIT to cluster sequences at given identity threshold."""
    print()
    print("=" * 70)
    print(f"STEP 6: CD-HIT clustering at {identity*100:.0f}% identity")
    print("=" * 70)

    cmd = [
        'conda', 'run', '-n', 'phylo_asr',
        'cd-hit',
        '-i', input_fasta,
        '-o', output_fasta,
        '-c', str(identity),
        '-n', '5',
        '-M', '4000',  # memory limit
        '-T', '8',     # threads
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  CD-HIT failed: {result.stderr}")
        return 0

    # count output sequences
    n_out = 0
    with open(output_fasta) as f:
        for line in f:
            if line.startswith('>'):
                n_out += 1

    print(f"  input: {sum(1 for line in open(input_fasta) if line.startswith('>'))} sequences")
    print(f"  output: {n_out} representative sequences at {identity*100:.0f}% identity")

    return n_out


def selectFinalSet(cdhit_fasta, domain_sequences, target_n=70):
    """from CD-HIT output, select a balanced set ensuring taxonomic diversity.
    mandatory PDB proteins always included.
    """
    print()
    print("=" * 70)
    print(f"STEP 7: select final balanced set (target ~{target_n} sequences)")
    print("=" * 70)

    # read CD-HIT representatives
    cdhit_ids = set()
    with open(cdhit_fasta) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].split()[0]
                cdhit_ids.add(seq_id)

    print(f"  CD-HIT representatives: {len(cdhit_ids)}")

    # classify by domain of life
    by_domain = defaultdict(list)
    mandatory_keys = set()

    for key, meta in domain_sequences.items():
        # reconstruct the FASTA header to match CD-HIT output
        org_clean = meta['organism'].replace(' ', '_').replace('(', '').replace(')', '')
        has_pdb_str = meta['has_pdb'] if meta['has_pdb'] else 'no'
        fasta_id = f"{org_clean}|{meta['accession']}_KH{meta['kh_index']}|{meta['domain_of_life']}|{has_pdb_str}"

        if fasta_id in cdhit_ids:
            by_domain[meta['domain_of_life']].append((key, meta, fasta_id))

        # mark mandatory
        if meta['has_pdb']:
            mandatory_keys.add(key)

    for d in ['Archaea', 'Bacteria', 'Eukaryota', 'Unknown']:
        n = len(by_domain.get(d, []))
        if n > 0:
            print(f"    {d}: {n} CD-HIT representatives")

    # selection strategy:
    # 1. always include mandatory PDB proteins
    # 2. balance across domains of life
    # 3. prefer reviewed/Swiss-Prot over TrEMBL

    selected = {}

    # add all mandatory first
    for key, meta in domain_sequences.items():
        if meta['has_pdb']:
            org_clean = meta['organism'].replace(' ', '_').replace('(', '').replace(')', '')
            has_pdb_str = meta['has_pdb']
            fasta_id = f"{org_clean}|{meta['accession']}_KH{meta['kh_index']}|{meta['domain_of_life']}|{has_pdb_str}"
            selected[fasta_id] = meta
            print(f"    mandatory: {meta['accession']} ({meta['protein_name']}, {meta['has_pdb']})")

    # allocate remaining slots proportionally
    remaining = target_n - len(selected)
    targets = {
        'Archaea': min(remaining // 3, len(by_domain.get('Archaea', []))),
        'Bacteria': min(remaining // 3, len(by_domain.get('Bacteria', []))),
        'Eukaryota': min(remaining // 3, len(by_domain.get('Eukaryota', []))),
    }
    # fill any surplus from the domain with most available
    surplus = remaining - sum(targets.values())
    for d in sorted(by_domain.keys(), key=lambda x: len(by_domain.get(x, [])), reverse=True):
        if d not in targets:
            continue
        can_add = len(by_domain[d]) - targets[d]
        add_now = min(surplus, can_add)
        targets[d] += add_now
        surplus -= add_now
        if surplus <= 0:
            break

    for domain, items in by_domain.items():
        if domain not in targets:
            continue
        target_for_domain = targets[domain]
        added = 0
        for key, meta, fasta_id in items:
            if fasta_id in selected:
                added += 1  # already counted (mandatory)
                continue
            if added >= target_for_domain:
                break
            selected[fasta_id] = meta
            added += 1

    print(f"\n  selected {len(selected)} sequences total:")
    domain_counts = defaultdict(int)
    for fasta_id, meta in selected.items():
        domain_counts[meta['domain_of_life']] += 1
    for d in ['Archaea', 'Bacteria', 'Eukaryota', 'Unknown']:
        if domain_counts[d] > 0:
            print(f"    {d}: {domain_counts[d]}")

    # write final FASTA
    with open(FINAL_FASTA, 'w') as f:
        for fasta_id, meta in selected.items():
            f.write(f'>{fasta_id}\n')
            seq = meta['sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

    print(f"\n  wrote final set to {FINAL_FASTA}")

    return selected


def main():
    os.makedirs(KH_DIR, exist_ok=True)

    # check seed file exists
    if not os.path.exists(SEED_STO):
        print(f"ERROR: seed alignment not found at {SEED_STO}")
        print("run: wget https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00013/?annotation=alignment:seed")
        sys.exit(1)

    # build the dataset
    domain_sequences = buildDiverseDataset()

    # run CD-HIT
    n_clustered = runCDHit(RAW_FASTA, CDHIT_FASTA, identity=0.9)

    if n_clustered == 0:
        print("CD-HIT failed or produced no output. using raw sequences.")
        # copy raw to cdhit output
        import shutil
        shutil.copy(RAW_FASTA, CDHIT_FASTA)

    # select final balanced set
    selected = selectFinalSet(CDHIT_FASTA, domain_sequences, target_n=70)

    print()
    print("=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"  seed FASTA:    {SEED_FASTA}")
    print(f"  raw FASTA:     {RAW_FASTA}")
    print(f"  CD-HIT FASTA:  {CDHIT_FASTA}")
    print(f"  final FASTA:   {FINAL_FASTA}")
    print(f"  metadata:      {METADATA_TSV}")
    print(f"  final count:   {len(selected)} sequences")


if __name__ == '__main__':
    main()
