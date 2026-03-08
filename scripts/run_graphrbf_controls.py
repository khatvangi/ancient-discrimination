#!/usr/bin/env python3
"""
run graphrbf on control proteins for RNA vs DNA binding discrimination.

uses graphrbf (Wan et al. 2024, GigaScience) which has separate models for
DNA-binding and RNA-binding residue prediction. we compare the fraction and
mean probability of DNA-binding vs RNA-binding residues per protein.

note: this version uses UNIFORM (non-informative) PSSM and HHM features
because we don't have the BLAST/HHblits databases formatted yet.
the model was trained with real profiles, so absolute scores will be degraded,
but the RELATIVE comparison between DNA and RNA scores may still be informative.

usage:
    conda activate graphrbf
    python run_graphrbf_controls.py
"""

import os
import sys
import shutil
import subprocess
import pickle
import numpy as np
import pandas as pd
import torch
from itertools import repeat, product

# paths
GRAPHRBF_DIR = "/storage/kiran-stuff/GraphRBF/GraphRBF-main"
CONTROLS_DIR = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/pnabind_controls"
WORK_DIR = "/storage/kiran-stuff/Vocabulary-Proteins/ancient-discrimination/structures/graphrbf_controls"
DSSP_BIN = shutil.which("mkdssp") or "/home/kiran/miniforge3/bin/mkdssp"

sys.path.insert(0, GRAPHRBF_DIR)
os.chdir(GRAPHRBF_DIR)

from GN_model_guassian_posemb import GraphRBF as GraphRBFModel

# define FocalLoss in __main__ so torch.load can unpickle the saved model
import torch.nn as nn
class FocalLoss(nn.Module):
    def __init__(self, alpha=0.25, gamma=2, logits=False, sampling='mean'):
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.logits = logits
        self.sampling = sampling

    def forward(self, y_pred, y_true):
        alpha = self.alpha
        alpha_ = (1 - self.alpha)
        if self.logits:
            y_pred = torch.sigmoid(y_pred)
        pt_positive = torch.where(y_true == 1, y_pred, torch.ones_like(y_pred))
        pt_negative = torch.where(y_true == 0, y_pred, torch.zeros_like(y_pred))
        pt_positive = torch.clamp(pt_positive, 1e-3, .999)
        pt_negative = torch.clamp(pt_negative, 1e-3, .999)
        pos_ = (1 - pt_positive) ** self.gamma
        neg_ = pt_negative ** self.gamma
        pos_loss = -alpha * pos_ * torch.log(pt_positive)
        neg_loss = -alpha_ * neg_ * torch.log(1 - pt_negative)
        loss = pos_loss + neg_loss
        if self.sampling == "mean":
            return loss.mean()
        elif self.sampling == "sum":
            return loss.sum()
        elif self.sampling is None:
            return loss

# control proteins with known binding preferences
CONTROLS = {
    "P0A9X9": {"name": "CspA E.coli", "known": "dual RNA/DNA"},
    "P62244": {"name": "RPS15A human", "known": "RNA"},
    "P03023": {"name": "LacI E.coli", "known": "DNA"},
    "P04637": {"name": "TP53 human", "known": "DNA"},
    "P67809": {"name": "YBX1 human", "known": "dual RNA/DNA"},
    "P69441": {"name": "ADK E.coli", "known": "non-binder"},
    "U1A_crystal": {"name": "U1A snRNP", "known": "RNA"},
    "MS2_crystal": {"name": "MS2 coat protein", "known": "RNA"},
    "LambdaRep_crystal": {"name": "Lambda repressor", "known": "DNA"},
}


def get_sequence_from_pdb(pdb_path):
    """extract amino acid sequence from PDB ATOM records."""
    res_dict = {
        'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'ILE': 'I', 'LEU': 'L',
        'PHE': 'F', 'PRO': 'P', 'MET': 'M', 'TRP': 'W', 'CYS': 'C',
        'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'TYR': 'Y',
        'HIS': 'H', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K', 'ARG': 'R'
    }
    seen = set()
    sequence = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = line[22:27].strip()
                key = (chain, res_num)
                if key not in seen and res_name in res_dict:
                    seen.add(key)
                    sequence.append(res_dict[res_name])
    return "".join(sequence)


def generate_dummy_pssm(sequence, output_path):
    """generate a uniform PSSM file (all zeros) in PSI-BLAST format."""
    aa_order = "ARNDCQEGHILKMFPSTWYV"
    lines = [
        "",
        "Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts",
        "            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V",
    ]
    for i, aa in enumerate(sequence):
        # all zeros for PSSM scores
        scores = "   0" * 20
        # uniform percentages (5% each)
        percs = "   5" * 20
        lines.append(f"  {i+1:3d} {aa}  {scores}{percs}  0.00 0.00")
    lines.append("")
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def generate_dummy_hhm(sequence, output_path, name="protein"):
    """generate a minimal HHM file with uniform transition probabilities."""
    lines = [
        "HHsearch 1.5",
        f"NAME  {name}",
        "FAM   ",
        "COM   dummy",
        "DATE  Thu Jan 01 00:00:00 2026",
        f"LENG  {len(sequence)} match states, {len(sequence)} columns in multiple alignment",
        "",
        "FILT  1 out of 1 sequences passed filter",
        "NEFF  1.0 ",
        "SEQ",
        ">Consensus",
        sequence,
        f">{name}",
        sequence,
        "#",
    ]
    # header for HMM section
    lines.append("NULL   3706	5728	4211	4064	4839	3729	4763	4308	4069	3323	5509	4640	4464	4937	4285	4423	3815	6325	4665	3912	")
    lines.append("HMM    A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	")
    lines.append("       M->M	M->I	M->D	I->M	I->I	D->M	D->D	Neff	Neff_I	Neff_D")
    lines.append("       0	*	*	0	*	0	*	*	*	*")

    for i, aa in enumerate(sequence):
        # emission probabilities - uniform (all same value ~3000)
        emissions = "\t".join(["3000"] * 20)
        lines.append(f"{aa} {i+1}	{emissions}")
        # transition probabilities
        lines.append("\t0\t*\t*\t0\t*\t0\t*\t1000\t1000\t1000")
        lines.append("")

    lines.append("//")
    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def run_dssp(pdb_path, dssp_path):
    """run DSSP on a PDB file. handles both old (-i/-o) and new (positional) syntax."""
    # try new syntax first (mkdssp v4+): mkdssp input output --output-format dssp
    result = subprocess.run(
        [DSSP_BIN, "--output-format", "dssp", pdb_path, dssp_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        # try old syntax: mkdssp -i input -o output
        result = subprocess.run(
            [DSSP_BIN, "-i", pdb_path, "-o", dssp_path],
            capture_output=True, text=True
        )
    if result.returncode != 0:
        print(f"  DSSP warning: {result.stderr[:200]}")
    return os.path.exists(dssp_path)


def clean_pdb_for_graphrbf(input_pdb, output_pdb):
    """clean PDB: keep only standard amino acid ATOM records, single chain."""
    valid_res = {
        'GLY', 'ALA', 'VAL', 'ILE', 'LEU', 'PHE', 'PRO', 'MET', 'TRP', 'CYS',
        'SER', 'THR', 'ASN', 'GLN', 'TYR', 'HIS', 'ASP', 'GLU', 'LYS', 'ARG'
    }
    lines = []
    first_chain = None
    with open(input_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                chain = line[21]
                if first_chain is None:
                    first_chain = chain
                if chain == first_chain and res_name in valid_res:
                    lines.append(line)
    lines.append("TER\n")
    lines.append("END\n")

    with open(output_pdb, "w") as f:
        f.writelines(lines)
    return first_chain or "A"


def prepare_protein(protein_id, pdb_src, work_dir):
    """prepare all input files for one protein."""
    os.makedirs(work_dir, exist_ok=True)

    # clean PDB and get chain
    clean_pdb = os.path.join(work_dir, f"{protein_id}_A.pdb")
    chain = clean_pdb_for_graphrbf(pdb_src, clean_pdb)

    # also create the 'protein.pdb' that GraphRBF expects internally
    protein_pdb = os.path.join(work_dir, "protein.pdb")
    shutil.copy(clean_pdb, protein_pdb)

    # get sequence
    seq = get_sequence_from_pdb(clean_pdb)
    if len(seq) == 0:
        print(f"  ERROR: no standard residues found in {pdb_src}")
        return None

    # generate dummy PSSM
    pssm_path = os.path.join(work_dir, f"{protein_id}_A.pssm")
    generate_dummy_pssm(seq, pssm_path)

    # generate dummy HHM
    hhm_path = os.path.join(work_dir, f"{protein_id}_A.hhm")
    generate_dummy_hhm(seq, hhm_path, name=f"{protein_id}_A")

    # run DSSP
    dssp_path = os.path.join(work_dir, "protein.dssp")
    if not run_dssp(protein_pdb, dssp_path):
        print(f"  ERROR: DSSP failed for {protein_id}")
        return None

    return {"seq": seq, "chain": "A", "n_residues": len(seq)}


def norm_DSSP_v4(query_path, protein):
    """parse DSSP output with extended SS codes (handles P from DSSP v4).
    'P' (polyproline) is mapped to the same encoding as ' ' (coil)."""
    maxASA = {'G': 188, 'A': 198, 'V': 220, 'I': 233, 'L': 304, 'F': 272,
              'P': 203, 'M': 262, 'W': 317, 'C': 201, 'S': 234, 'T': 215,
              'N': 254, 'Q': 259, 'Y': 304, 'H': 258, 'D': 236, 'E': 262,
              'K': 317, 'R': 319}
    map_ss_8 = {
        ' ': [1, 0, 0, 0, 0, 0, 0, 0], 'S': [0, 1, 0, 0, 0, 0, 0, 0],
        'T': [0, 0, 1, 0, 0, 0, 0, 0], 'H': [0, 0, 0, 1, 0, 0, 0, 0],
        'G': [0, 0, 0, 0, 1, 0, 0, 0], 'I': [0, 0, 0, 0, 0, 1, 0, 0],
        'E': [0, 0, 0, 0, 0, 0, 1, 0], 'B': [0, 0, 0, 0, 0, 0, 0, 1],
        # DSSP v4 extensions: map to closest standard category
        'P': [1, 0, 0, 0, 0, 0, 0, 0],  # polyproline -> coil
        'C': [1, 0, 0, 0, 0, 0, 0, 0],  # coil explicit -> coil
    }
    file_path = '{}/{}.dssp'.format(query_path, protein)
    with open(file_path, 'r') as f:
        text = f.readlines()

    start_line = 0
    for i in range(len(text)):
        if text[i].startswith('  #  RESIDUE AA STRUCTURE'):
            start_line = i + 1
            break

    dssp = {}
    for i in range(start_line, len(text)):
        line = text[i]
        if len(line) < 115:
            continue
        if line[13] not in maxASA.keys() or line[9] == ' ':
            continue
        ss_code = line[16]
        if ss_code not in map_ss_8:
            ss_code = ' '  # fallback to coil
        res_id = f'{line[11]}{line[5:11].strip()}'
        res_dssp = np.zeros([14])
        res_dssp[:8] = map_ss_8[ss_code]
        try:
            res_dssp[8] = min(float(line[35:38]) / maxASA[line[13]], 1)
        except (ValueError, ZeroDivisionError):
            res_dssp[8] = 0.5
        try:
            res_dssp[9] = (float(line[85:91]) + 1) / 2
            res_dssp[10] = min(1, float(line[91:97]) / 180)
            res_dssp[11] = min(1, (float(line[97:103]) + 180) / 360)
            res_dssp[12] = min(1, (float(line[103:109]) + 180) / 360)
            res_dssp[13] = min(1, (float(line[109:115]) + 180) / 360)
        except (ValueError, IndexError):
            pass
        dssp[res_id] = res_dssp.reshape((1, -1))

    return dssp


def run_graphrbf_prediction(protein_id, work_dir):
    """run graphrbf prediction for both DNA and RNA ligands.
    returns dict with per-residue scores for each ligand type."""

    # import graphrbf functions
    from GraphRBF import (
        PDBFeature, norm_pssm, norm_hhm,
        PDBResidueFeature, seq_Dataset, predict
    )

    protein = "protein"

    # extract features
    onehot = PDBFeature(protein, work_dir, work_dir)
    pssm = norm_pssm(work_dir, protein_id, "A")
    hhm = norm_hhm(work_dir, protein_id, "A")
    dssp = norm_DSSP_v4(work_dir, protein)

    PDBResidueFeature(work_dir, protein, onehot, pssm, hhm, dssp)

    results = {}
    for ligand_name, model_subdir in [("DNA", "PDNA"), ("RNA", "PRNA")]:
        dist = 20
        query_data = seq_Dataset(work_dir, protein, dist)
        model = GraphRBFModel(
            gnn_steps=1, x_ind=92, edge_ind=2, x_hs=128, e_hs=128, u_hs=128,
            dropratio=0.5, bias=True, r_list=[10], dist=dist, max_nn=40
        )
        model_path = os.path.join(GRAPHRBF_DIR, "models", model_subdir, "model", "model0.pth")
        threshold, pred_score, pred_binary = predict(model, model_path, query_data)

        results[ligand_name] = {
            "threshold": float(threshold),
            "scores": pred_score.tolist(),
            "binary": pred_binary.tolist(),
            "residue_ids": query_data.res_id_list,
            "sequence": query_data.sequence,
            "mean_score": float(np.mean(pred_score)),
            "frac_binding": float(np.mean(pred_binary)),
            "n_binding": int(np.sum(pred_binary)),
        }

    return results


def compute_discrimination_index(results):
    """compute a discrimination index from DNA vs RNA predictions.

    DI = (mean_DNA_score - mean_RNA_score) / (mean_DNA_score + mean_RNA_score)
    DI > 0 means more DNA-binding character
    DI < 0 means more RNA-binding character
    DI ~ 0 means dual or non-specific
    """
    dna_mean = results["DNA"]["mean_score"]
    rna_mean = results["RNA"]["mean_score"]

    if dna_mean + rna_mean == 0:
        return 0.0

    di = (dna_mean - rna_mean) / (dna_mean + rna_mean)
    return di


def main():
    os.makedirs(WORK_DIR, exist_ok=True)
    os.chdir(GRAPHRBF_DIR)

    all_results = {}
    summary_rows = []

    print("=" * 80)
    print("GraphRBF DNA vs RNA Binding Discrimination Test")
    print("NOTE: using uniform (non-informative) PSSM/HHM features")
    print("=" * 80)

    for protein_id, info in CONTROLS.items():
        print(f"\n--- {protein_id}: {info['name']} (known: {info['known']}) ---")

        pdb_src = os.path.join(CONTROLS_DIR, f"{protein_id}.pdb")
        if not os.path.exists(pdb_src):
            print(f"  SKIP: {pdb_src} not found")
            continue

        protein_work_dir = os.path.join(WORK_DIR, protein_id)

        # step 1: prepare
        print("  preparing input files...")
        prep = prepare_protein(protein_id, pdb_src, protein_work_dir)
        if prep is None:
            continue

        print(f"  sequence length: {prep['n_residues']} residues")

        # step 2: run prediction
        print("  running GraphRBF prediction (DNA + RNA models)...")
        try:
            results = run_graphrbf_prediction(protein_id, protein_work_dir)
        except Exception as e:
            print(f"  ERROR: prediction failed: {e}")
            import traceback
            traceback.print_exc()
            continue

        # step 3: compute discrimination index
        di = compute_discrimination_index(results)

        all_results[protein_id] = results

        row = {
            "protein": protein_id,
            "name": info["name"],
            "known_binding": info["known"],
            "n_residues": prep["n_residues"],
            "DNA_mean_score": results["DNA"]["mean_score"],
            "RNA_mean_score": results["RNA"]["mean_score"],
            "DNA_frac_binding": results["DNA"]["frac_binding"],
            "RNA_frac_binding": results["RNA"]["frac_binding"],
            "DNA_n_binding": results["DNA"]["n_binding"],
            "RNA_n_binding": results["RNA"]["n_binding"],
            "discrimination_index": di,
        }
        summary_rows.append(row)

        print(f"  DNA: mean_score={results['DNA']['mean_score']:.3f}, "
              f"frac_binding={results['DNA']['frac_binding']:.3f} "
              f"({results['DNA']['n_binding']}/{prep['n_residues']})")
        print(f"  RNA: mean_score={results['RNA']['mean_score']:.3f}, "
              f"frac_binding={results['RNA']['frac_binding']:.3f} "
              f"({results['RNA']['n_binding']}/{prep['n_residues']})")
        print(f"  Discrimination Index (DI): {di:+.3f}")
        if di > 0.05:
            print(f"  -> PREDICTED: DNA-binding")
        elif di < -0.05:
            print(f"  -> PREDICTED: RNA-binding")
        else:
            print(f"  -> PREDICTED: dual/non-specific")

    # summary table
    if summary_rows:
        df = pd.DataFrame(summary_rows)
        output_csv = os.path.join(WORK_DIR, "graphrbf_discrimination_results.csv")
        df.to_csv(output_csv, index=False, float_format="%.4f")

        print("\n" + "=" * 80)
        print("SUMMARY TABLE")
        print("=" * 80)
        print(f"{'Protein':<20} {'Known':<15} {'DNA_mean':>10} {'RNA_mean':>10} {'DI':>8} {'Predicted':<15}")
        print("-" * 80)
        for _, row in df.iterrows():
            di = row["discrimination_index"]
            if di > 0.05:
                pred = "DNA"
            elif di < -0.05:
                pred = "RNA"
            else:
                pred = "dual/non-spec"
            print(f"{row['protein']:<20} {row['known_binding']:<15} "
                  f"{row['DNA_mean_score']:>10.3f} {row['RNA_mean_score']:>10.3f} "
                  f"{di:>+8.3f} {pred:<15}")

        print(f"\nResults saved to: {output_csv}")

        # count correct predictions
        correct = 0
        total = 0
        for _, row in df.iterrows():
            known = row["known_binding"]
            di = row["discrimination_index"]
            total += 1
            if known == "DNA" and di > 0.05:
                correct += 1
            elif known == "RNA" and di < -0.05:
                correct += 1
            elif "dual" in known and abs(di) <= 0.05:
                correct += 1
            elif known == "non-binder":
                correct += 1  # we don't penalize non-binders

        print(f"\nAccuracy: {correct}/{total} = {correct/total*100:.0f}%")
        print("(note: with dummy PSSM/HHM, accuracy will be lower than with real profiles)")


if __name__ == "__main__":
    main()
