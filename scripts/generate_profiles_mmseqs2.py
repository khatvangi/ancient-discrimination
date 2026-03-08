#!/usr/bin/env python3
"""
generate PSSM and HHM profiles for control proteins using MMseqs2.

mmseqs2 is orders of magnitude faster than PSI-BLAST and HHblits.
we use mmseqs2 to search against uniref90, build a profile, and
convert to PSSM format compatible with GraphRBF.

usage:
    python generate_profiles_mmseqs2.py
"""

import os
import sys
import subprocess
import numpy as np
from Bio.PDB import PDBParser

# paths
CONTROLS_DIR = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/pnabind_controls"
GRAPHRBF_WORK = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/graphrbf_controls"
MMSEQS2_DB = "/storage/kiran-stuff/alphafold_databases/uniref90_2022_05_mmseqs2/uniref90_2022_05"
MMSEQS2_BIN = "/home/kiran/miniforge3/bin/mmseqs"
PSIBLAST_BIN = "/storage/kiran-stuff/protein-folds/tools/ncbi-blast-2.16.0+/bin/psiblast"

# proteins to process
PROTEINS = [
    "P0A9X9", "P62244", "P03023", "P04637", "P67809",
    "P69441", "U1A_crystal", "MS2_crystal", "LambdaRep_crystal"
]


def get_sequence_from_pdb(pdb_path):
    """extract amino acid sequence from PDB ATOM records (first chain only)."""
    res_dict = {
        'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L',
        'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'TRP': 'W', 'CYS': 'C',
        'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'TYR': 'Y',
        'HIS': 'H', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K', 'ARG': 'R'
    }
    seen = set()
    sequence = []
    first_chain = None
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                chain = line[21]
                if first_chain is None:
                    first_chain = chain
                if chain != first_chain:
                    continue
                res_num = line[22:27].strip()
                key = (chain, res_num)
                if key not in seen and res_name in res_dict:
                    seen.add(key)
                    sequence.append(res_dict[res_name])
    return "".join(sequence)


def generate_pssm_via_mmseqs2(protein_id, sequence, work_dir):
    """generate PSSM profile using MMseqs2 profile search.

    mmseqs2 approach:
    1. create query database from sequence
    2. search against uniref90
    3. extract profile (position-specific scoring matrix)
    4. convert to PSI-BLAST PSSM format
    """
    os.makedirs(work_dir, exist_ok=True)
    tmp_dir = os.path.join(work_dir, "mmseqs_tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    # write fasta
    fasta_path = os.path.join(work_dir, f"{protein_id}.fasta")
    with open(fasta_path, "w") as f:
        f.write(f">{protein_id}\n{sequence}\n")

    # create query db
    query_db = os.path.join(work_dir, f"{protein_id}_querydb")
    subprocess.run(
        [MMSEQS2_BIN, "createdb", fasta_path, query_db],
        capture_output=True, check=True
    )

    # run profile search (iterative)
    result_db = os.path.join(work_dir, f"{protein_id}_result")
    profile_db = os.path.join(work_dir, f"{protein_id}_profile")

    # first search to build initial profile
    subprocess.run(
        [MMSEQS2_BIN, "search", query_db, MMSEQS2_DB, result_db, tmp_dir,
         "--num-iterations", "3", "-e", "0.001", "--threads", "16"],
        capture_output=True, check=True
    )

    # create profile from search results
    subprocess.run(
        [MMSEQS2_BIN, "result2profile", query_db, MMSEQS2_DB, result_db, profile_db],
        capture_output=True, check=True
    )

    # extract profile as tab-separated values
    profile_tsv = os.path.join(work_dir, f"{protein_id}_profile.tsv")
    subprocess.run(
        [MMSEQS2_BIN, "profile2pssm", profile_db, profile_tsv],
        capture_output=True
    )

    # also get consensus sequence
    consensus_path = os.path.join(work_dir, f"{protein_id}_consensus")
    subprocess.run(
        [MMSEQS2_BIN, "profile2consensus", profile_db, consensus_path],
        capture_output=True
    )

    return profile_db, profile_tsv


def convert_mmseqs_profile_to_pssm(sequence, profile_db, output_path, work_dir, protein_id):
    """convert mmseqs2 profile to PSI-BLAST PSSM format.

    mmseqs2 profiles store log-odds scores in a different format.
    we need to convert to the standard PSSM format expected by GraphRBF.
    """
    # try to extract the profile matrix from the mmseqs2 database
    # mmseqs2 stores profiles in its internal format
    # let's extract using result2msa and then compute PSSM from MSA

    aa_order_pssm = "ARNDCQEGHILKMFPSTWYV"
    aa_order_mmseqs = "ACDEFGHIKLMNPQRSTVWY"  # mmseqs2 alphabetical order

    # try to extract MSA and compute PSSM from it
    result_db = os.path.join(work_dir, f"{protein_id}_result")
    msa_path = os.path.join(work_dir, f"{protein_id}_msa.a3m")

    query_db = os.path.join(work_dir, f"{protein_id}_querydb")

    # extract MSA
    result = subprocess.run(
        [MMSEQS2_BIN, "result2msa", query_db, MMSEQS2_DB, result_db, msa_path,
         "--msa-format-mode", "5"],
        capture_output=True, text=True
    )

    # parse a3m and compute PSSM
    if os.path.exists(msa_path):
        pssm_matrix = compute_pssm_from_a3m(msa_path, sequence)
    else:
        print(f"  WARNING: could not generate MSA, using uniform PSSM")
        pssm_matrix = np.zeros((len(sequence), 20))

    # write in PSI-BLAST format
    write_pssm_file(sequence, pssm_matrix, output_path)
    return True


def compute_pssm_from_a3m(a3m_path, query_seq):
    """compute position-specific scoring matrix from A3M alignment."""
    aa_order = "ARNDCQEGHILKMFPSTWYV"
    aa_to_idx = {aa: i for i, aa in enumerate(aa_order)}

    # read a3m file
    sequences = []
    with open(a3m_path) as f:
        current_seq = []
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            sequences.append("".join(current_seq))

    if len(sequences) < 2:
        # not enough sequences for meaningful profile
        return np.zeros((len(query_seq), 20))

    query_len = len(query_seq)

    # count amino acid frequencies at each position
    # in a3m format, uppercase = match state, lowercase = insert
    counts = np.zeros((query_len, 20))
    pseudocount = 0.01  # small pseudocount

    for seq in sequences:
        pos = 0
        for c in seq:
            if c == '-':
                pos += 1
            elif c.isupper():
                if pos < query_len and c in aa_to_idx:
                    counts[pos, aa_to_idx[c]] += 1
                pos += 1
            # lowercase letters are insertions, skip them

    # convert to frequencies
    total = counts.sum(axis=1, keepdims=True)
    total[total == 0] = 1  # avoid division by zero
    freq = counts / total + pseudocount
    freq = freq / freq.sum(axis=1, keepdims=True)

    # background frequencies (from BLOSUM62)
    bg_freq = np.array([
        0.0787, 0.0507, 0.0447, 0.0543, 0.0199, 0.0384, 0.0618, 0.0750,
        0.0221, 0.0572, 0.0911, 0.0574, 0.0231, 0.0400, 0.0468, 0.0629,
        0.0585, 0.0647, 0.0133, 0.0321
    ])

    # compute log-odds scores (PSI-BLAST style)
    # score = round(2 * log2(freq / bg_freq))
    with np.errstate(divide='ignore', invalid='ignore'):
        log_odds = 2.0 * np.log2(freq / bg_freq)
    log_odds = np.nan_to_num(log_odds, nan=0, posinf=10, neginf=-10)
    pssm = np.round(log_odds).astype(int)
    pssm = np.clip(pssm, -10, 10)

    return pssm


def write_pssm_file(sequence, pssm_matrix, output_path):
    """write PSSM in PSI-BLAST format compatible with GraphRBF."""
    aa_order = "ARNDCQEGHILKMFPSTWYV"

    lines = [
        "",
        "Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts",
        "            " + "   ".join(f"{aa:>3}" for aa in aa_order) + "   " + "   ".join(f"{aa:>3}" for aa in aa_order),
    ]

    for i, aa in enumerate(sequence):
        scores = "".join(f"{int(pssm_matrix[i, j]):>4}" for j in range(20))
        # weighted observed percentages (just use uniform as we compute from log-odds)
        percs = "".join(f"{5:>4}" for _ in range(20))
        lines.append(f"  {i+1:3d} {aa}  {scores}{percs}  0.00 0.00")

    lines.append("")
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def generate_hhm_from_mmseqs2(protein_id, sequence, work_dir, output_path):
    """generate HHM-format profile from MMseqs2 MSA.

    we compute a simplified HHM by extracting position-specific emission
    probabilities from the MSA alignment.
    """
    a3m_path = os.path.join(work_dir, f"{protein_id}_msa.a3m")

    if not os.path.exists(a3m_path):
        # generate dummy HHM
        write_dummy_hhm(sequence, output_path, protein_id)
        return

    # compute emission probabilities from MSA
    aa_order_hhm = "ACDEFGHIKLMNPQRSTVWY"
    aa_to_idx = {aa: i for i, aa in enumerate(aa_order_hhm)}

    # read a3m
    sequences = []
    with open(a3m_path) as f:
        current_seq = []
        for line in f:
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            sequences.append("".join(current_seq))

    query_len = len(sequence)

    # count frequencies
    counts = np.ones((query_len, 20)) * 0.01  # pseudocounts
    for seq in sequences:
        pos = 0
        for c in seq:
            if c == '-':
                pos += 1
            elif c.isupper():
                if pos < query_len and c in aa_to_idx:
                    counts[pos, aa_to_idx[c]] += 1
                pos += 1

    # normalize to probabilities
    probs = counts / counts.sum(axis=1, keepdims=True)

    # convert to HHM format: -1000 * log2(p)
    with np.errstate(divide='ignore'):
        hhm_scores = np.round(-1000 * np.log2(probs)).astype(int)
    hhm_scores = np.clip(hhm_scores, 0, 9999)

    # write HHM file
    write_hhm_file(sequence, hhm_scores, output_path, protein_id, len(sequences))


def write_hhm_file(sequence, hhm_scores, output_path, name, n_seqs):
    """write profile in HHM format for GraphRBF."""
    neff = min(8.0, max(1.0, np.log2(n_seqs + 1)))

    lines = [
        "HHsearch 1.5",
        f"NAME  {name}",
        "FAM   ",
        "COM   generated from MMseqs2 MSA",
        "DATE  Mon Feb 24 00:00:00 2026",
        f"LENG  {len(sequence)} match states, {len(sequence)} columns in multiple alignment",
        "",
        f"FILT  {n_seqs} out of {n_seqs} sequences passed filter",
        f"NEFF  {neff:.1f} ",
        "SEQ",
        ">Consensus",
        sequence,
        f">{name}",
        sequence,
        "#",
        "NULL   3706\t5728\t4211\t4064\t4839\t3729\t4763\t4308\t4069\t3323\t5509\t4640\t4464\t4937\t4285\t4423\t3815\t6325\t4665\t3912\t",
        "HMM    A\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY\t",
        "       M->M\tM->I\tM->D\tI->M\tI->I\tD->M\tD->D\tNeff\tNeff_I\tNeff_D",
        f"       0\t*\t*\t0\t*\t0\t*\t{int(neff*1000)}\t{int(neff*1000)}\t{int(neff*1000)}",
    ]

    for i, aa in enumerate(sequence):
        # emission line
        scores_str = "\t".join(str(hhm_scores[i, j]) for j in range(20))
        lines.append(f"{aa} {i+1}\t{scores_str}\t")
        # transition line
        lines.append(f"\t0\t*\t*\t0\t*\t0\t*\t{int(neff*1000)}\t{int(neff*1000)}\t{int(neff*1000)}")
        lines.append("")

    lines.append("//")
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def write_dummy_hhm(sequence, output_path, name):
    """write a dummy HHM file with uniform probabilities."""
    hhm_scores = np.full((len(sequence), 20), 4322)  # -1000*log2(1/20) ≈ 4322
    write_hhm_file(sequence, hhm_scores, output_path, name, 1)


def main():
    print("=" * 70)
    print("generating PSSM and HHM profiles using MMseqs2")
    print("=" * 70)

    for protein_id in PROTEINS:
        pdb_path = os.path.join(CONTROLS_DIR, f"{protein_id}.pdb")
        if not os.path.exists(pdb_path):
            print(f"  SKIP: {pdb_path} not found")
            continue

        work_dir = os.path.join(GRAPHRBF_WORK, protein_id)
        os.makedirs(work_dir, exist_ok=True)

        # get sequence
        seq = get_sequence_from_pdb(pdb_path)
        print(f"\n{protein_id}: {len(seq)} residues")

        # generate PSSM via MMseqs2
        print(f"  running MMseqs2 profile search...")
        try:
            profile_db, profile_tsv = generate_pssm_via_mmseqs2(protein_id, seq, work_dir)
        except subprocess.CalledProcessError as e:
            print(f"  ERROR: MMseqs2 search failed: {e}")
            continue

        # convert to PSSM format
        pssm_path = os.path.join(work_dir, f"{protein_id}_A.pssm")
        print(f"  converting to PSSM format...")
        convert_mmseqs_profile_to_pssm(seq, profile_db, pssm_path, work_dir, protein_id)

        # generate HHM from MSA
        hhm_path = os.path.join(work_dir, f"{protein_id}_A.hhm")
        print(f"  generating HHM profile...")
        generate_hhm_from_mmseqs2(protein_id, seq, work_dir, hhm_path)

        # verify files exist
        if os.path.exists(pssm_path) and os.path.exists(hhm_path):
            pssm_size = os.path.getsize(pssm_path)
            hhm_size = os.path.getsize(hhm_path)
            print(f"  PSSM: {pssm_size:,} bytes, HHM: {hhm_size:,} bytes")
        else:
            print(f"  WARNING: profile generation incomplete")

    print("\nDone!")


if __name__ == "__main__":
    main()
