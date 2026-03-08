# Structural Figure Data Package
# generated 2026-02-26

## overview

this directory contains per-residue contact data for all PDB structures
needed in the structural visualization figures. every PyMOL selection
below comes directly from NAS-Bench contact analysis — no guessing.

## files

| file | structure | figure |
|------|-----------|--------|
| 3Q0P_contacts.tsv | PUF/Pumilio + RNA | Fig S1 (concept) |
| 4YYE_contacts.tsv | ThrRS + tRNA (PF03129) | Fig 3b (lysine claw) |
| 8F69_contacts.tsv | ThrRS + DNA (PF03129) | — |
| 1EC6_contacts.tsv | Nova-KH3 + RNA (PF00013) | Fig 5b, Fig S2 |
| 1ZTG_contacts.tsv | PCBP1-KH1 + DNA (PF00013) | Fig 5b, Fig S2 |
| 1MW8_contacts.tsv | Topoisomerase I + DNA (PF01131) | Fig S3 |
| 9GDA_contacts.tsv | Topoisomerase + RNA (PF01131) | Fig S3 |
| 1URN_contacts.tsv | U1A RRM + RNA | Fig S4 |
| 1LMB_contacts.tsv | Lambda repressor + DNA | Fig S4 |

each TSV has columns: chain, resnum, resname, aa1, n_bb, n_base, n_2oh, n_sr, n_total, dominant, na_type

## color scheme

| contact type | hex color | PyMOL color |
|-------------|-----------|-------------|
| backbone | #999999 | grey60 |
| base | #3498DB | tv_blue |
| 2'OH | #E74C3C | tv_red |
| sugar | #F1C40F | tv_yellow |
| mixed (base+bb) | #9B59B6 | purple |

## CRITICAL: lysine claw mapping (Fig 3b)

PDB structure: **4YYE** (threonyl-tRNA synthetase + tRNA, T. thermophilus)

### alignment column → PDB residue number

| aln_col | modern AA | PDB resnum | ancestral AA | PP | contact type | n_contacts |
|---------|-----------|------------|-------------|-----|-------------|------------|
| **67** | **E** | **401** | K | 0.51 | base-only | 13 base |
| **68** | **P** | **402** | K | 0.90 | 2'OH-contacting | 10 base + 2 2'OH |
| 69 | V | 403 | I | 0.57 | 2'OH-contacting | 4 2'OH |
| **70** | **G** | **404** | K | 0.92 | backbone-only | 3 bb |
| **71** | **Y** | **405** | K | 0.87 | base+backbone | 9 bb + 1 base |

### PyMOL commands for lysine claw

```pymol
# select the four claw positions
select claw_positions, resi 401+402+404+405 and chain A

# color by contact type
show sticks, claw_positions
color tv_red, resi 401 and chain A    # E: base H-bonding (charge reversal from K)
color firebrick, resi 402 and chain A  # P: 2'OH recognition
color grey60, resi 404 and chain A     # G: backbone / steric clearance
color purple, resi 405 and chain A     # Y: aromatic stacking (base+backbone)

# label with ancestral state
label resi 401 and name CA and chain A, "E401 (anc: K)"
label resi 402 and name CA and chain A, "P402 (anc: K)"
label resi 404 and name CA and chain A, "G404 (anc: K)"
label resi 405 and name CA and chain A, "Y405 (anc: K)"
```

### biological context

4YYE is a class II aminoacyl-tRNA synthetase (ThrRS). the HGTP_anticodon
domain (PF03129, ~residues 300-440) binds the anticodon loop of tRNA^Thr.
the four ancestral lysines at positions 401-405 (with V403 between G404
and P402) formed an electrostatic "claw" that gripped the phosphate backbone
nonspecifically. in modern ThrRS, these have diversified into:

- E401: charge reversal (acidic, H-bonds to base)
- P402: constrains backbone, positions 2'OH recognition
- G404: creates space for the anticodon to fit
- Y405: aromatic ring stacks with nucleobase

this is the mechanistic story: from generic electrostatic grip → specific
multimodal recognition.


## PyMOL selections for all structures

### 3Q0P (PUF + RNA, Fig S1)

```pymol
select bb_contacts, resi 897+932+933+1041+1076+1120
select base_contacts, resi 860+864+867+897+899+900+903+935+936+939+972+975+1007+1008+1011+1040+1043+1044+1047+1077+1079+1080+1083+1119+1122+1123+1126+1156+1159
select oh_contacts, resi 933+968+1076
select mixed_contacts, resi 936+1123
```

SI = 0.991 (almost pure base contacts). ideal for demonstrating the
NAS-Bench contact decomposition concept.


### 1EC6 (KH + RNA, Fig 5b left)

```pymol
select bb_contacts, resi 23+24+25+45
select base_contacts, resi 14+15+17+18+19+21+22+28+39+40+41+42+43+44+47+54+83+88
select oh_contacts, resi 15+23+29+38+45+52
```

SI = 0.932 (strong base preference). KH domain bound to RNA.


### 1ZTG (KH + DNA, Fig 5b right)

```pymol
select bb_contacts, resi 23+29+32+33+36
select base_contacts, resi 22+26+27+29+30+40+48+49+51+57+78+82
select mixed_contacts, resi 29+30+31
```

SI = 0.866 (also base-dominated). KH domain bound to DNA.
IMPORTANT: align protein chains before rendering:
  align 1ZTG and chain A, 1EC6 and chain A


### 1MW8 (Topoisomerase + DNA, Fig S3)

```pymol
select bb_contacts, resi 32+33+115+168+180+191+192+193+194+195+196+197+321+495+496+507
select base_contacts, resi 70+169+173+176+177+181+184+189+190+499
select mixed_contacts, resi 40+172
```

SI_DNA = 0.647. for the "frozen DNA machinery" panel, overlay
the 22/23 identical ancestral residues in green. the single
substitution (col 43, H→R, PP=0.42) should be in yellow.


### 1URN (U1A RRM + RNA, Fig S4 left)

```pymol
select bb_contacts, resi 22+46+48+51
select base_contacts, resi 5+15+16+19+44+52+53+54+56+80+85+86+87+88+89+90+91+92
select oh_contacts, resi 50
select mixed_contacts, resi 13+49+52+92
```

SI = 0.912. classic RNA specialist — overwhelmingly base contacts.


### 1LMB (Lambda repressor + DNA, Fig S4 right)

```pymol
select bb_contacts, resi 5+19+22+26+32+33+42+43+50+52+56+58+61
select base_contacts, resi 3+4+34+44+45+46+55
select mixed_contacts, resi 2+55
```

SI = 0.490. DNA specialist — balanced backbone + base contacts, zero 2'OH.


## notes for the student

1. overlay SI values as text in the corner of each rendered image
2. use identical camera orientations for paired RNA/DNA panels (align first)
3. every residue selection comes from the TSV files — verify by spot-checking
4. the PDB files for 3Q0P, 1URN, and 1LMB are in this directory;
   others are in results/asr/pfXXXXX/*.pdb
5. for Fig 3b, the tRNA chain in 4YYE is chain C (verify in PyMOL)
6. save .pse session files for each figure for future editing
