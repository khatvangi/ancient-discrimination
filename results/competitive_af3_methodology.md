# competitive AF3 binding: methodology and metrics guide

## design

- 3 chains per prediction: protein + ssRNA (15 nt) + ssDNA (15 nt)
- RNA sequence: GACUGAUUCGAUCAG
- DNA sequence: GACTGATTCGATCAG (same sequence, T/U swapped)
- equal chains = fair competition — AF3 must choose which substrate to dock

## output files per job

### summary_confidences_{N}.json (N = 0-4, one per model)

| metric | type | what it tells you |
|--------|------|-------------------|
| `iptm` | scalar | global interface TM-score — MISLEADING for large proteins (diluted by length) |
| `ptm` | scalar | predicted TM-score of protein fold quality |
| `ranking_score` | scalar | AF3's composite ranking — use to pick best model |
| `chain_iptm` | array[3] | per-chain iptm: [protein, RNA, DNA] |
| `chain_ptm` | array[3] | per-chain fold quality |
| `chain_pair_iptm` | 3x3 matrix | pairwise interface TM-scores — [0][1]=protein-RNA, [0][2]=protein-DNA |
| `chain_pair_pae_min` | 3x3 matrix | min predicted aligned error (Angstroms) — LOWER = BETTER |
| `fraction_disordered` | scalar | fraction of disordered residues |
| `has_clash` | scalar | steric clash indicator |
| `num_recycles` | scalar | number of recycling iterations used |

### full_data_{N}.json (per-residue details — THE REAL METRICS)

| metric | type | what it tells you |
|--------|------|-------------------|
| `contact_probs` | NxN matrix | probability of contact between every token pair |
| `pae` | NxN matrix | predicted aligned error in Angstroms per residue pair |
| `atom_plddts` | array | per-atom pLDDT confidence (0-100) |
| `token_chain_ids` | array | chain assignment per token (A=protein, B=RNA, C=DNA) |
| `token_res_ids` | array | residue number per token |
| `atom_chain_ids` | array | chain assignment per atom |

## metric hierarchy (most to least informative)

### 1. contact count (from contact_probs, PRIMARY)
```
protein-RNA contacts = count(contact_probs[protein_tokens, rna_tokens] > 0.5)
protein-DNA contacts = count(contact_probs[protein_tokens, dna_tokens] > 0.5)
```
- directly shows which substrate the protein is interacting with
- ratio > 1.5 → RNA preference; ratio < 0.67 → DNA preference; else BOTH/TIE
- if both = 0 → AF3 failed to dock either substrate (NO_DOCK)

### 2. per-chain pLDDT (from atom_plddts, QUALITY FILTER)
```
rna_plddt = mean(atom_plddts[atom_chain_ids == 'B'])
dna_plddt = mean(atom_plddts[atom_chain_ids == 'C'])
```
- pLDDT > 70: well-placed chain, confident docking
- pLDDT 50-70: moderate confidence
- pLDDT < 50: chain is essentially floating, not docked → UNINFORMATIVE
- use as quality filter: if both chains < 50, the prediction tells you nothing

### 3. chain_pair_pae_min (from summary_confidences, INTERFACE CONFIDENCE)
```
pae_protein_rna = chain_pair_pae_min[0][1]
pae_protein_dna = chain_pair_pae_min[0][2]
```
- lower PAE = higher confidence in the interface
- PAE < 3 Angstroms: high-confidence interface
- PAE 3-8: moderate
- PAE > 10: no confident interface
- more discriminating than iptm for large proteins

### 4. chain_pair_iptm (from summary_confidences, SUPPLEMENTARY)
```
iptm_protein_rna = chain_pair_iptm[0][1]
iptm_protein_dna = chain_pair_iptm[0][2]
discrimination_index = iptm_protein_rna - iptm_protein_dna
```
- threshold: |DI| > 0.02 to call preference
- CAUTION: ipTM is normalized by full protein length
- for large proteins (>700 aa) with a 15-mer probe, ipTM gets diluted
- can show "TIE" even when contacts clearly favor one substrate
- see: "Res ipSAE loquunt" (Dunbrack 2025) for ipTM limitations

## classification scheme

| call | criteria |
|------|----------|
| **RNA** | rna_contacts > dna_contacts × 1.5, rna_plddt > 50 |
| **DNA** | dna_contacts > rna_contacts × 1.5, dna_plddt > 50 |
| **BOTH** | both have contacts, ratio within 0.67-1.5 |
| **WEAK** | zero contacts but pLDDT > 50 for at least one chain |
| **NO_DOCK** | zero contacts AND pLDDT < 50 for both chains |

## known pitfalls

### 1. ipTM dilution for large proteins
- ipTM normalizes by protein length. a 1100 aa protein with 15 nt probe has only ~2% of tokens at the interface
- the global iptm can be 0.4-0.5 (looks decent) while zero contacts exist (actual: no docking)
- ALWAYS check contacts + pLDDT, never rely on ipTM alone

### 2. NO_DOCK ≠ "binds both equally"
- 62% of server predictions (>700 aa proteins) had zero contacts with either substrate
- these were initially classified as "TIE" using ipTM threshold — WRONG
- they mean "AF3 couldn't place the 15-mer on this large protein"
- expected: local batch (≤700 aa) should have better docking rates

### 3. probe sequence matters
- we use a generic 15-mer (GACTGATTCGATCAG / GACUGAUUCGAUCAG)
- proteins that need specific substrates (tRNA, mismatched DNA, specific motifs) may fail to dock
- tRNA synthetases: 6/7 NO_DOCK (need tRNA structure, not a 15-mer)
- MutS: all NO_DOCK or WEAK (needs mismatched dsDNA, not ssDNA 15-mer)

### 4. model selection
- AF3 server gives 5 models per job (model_0 through model_4)
- use the TOP-RANKED model (highest ranking_score) as primary reference
- server models are sorted by ranking_score (model_0 = best)
- for robust conclusions, check all 5 and report consistency

### 5. from aaRS project (competitive binding lessons)
- threshold for "good binding": iptm > 0.5
- "promiscuous" = both cognate AND competitor iptm > 0.5
- use Δ ipTM as primary, pocket IoU as structural confirmation
- always run negative controls to validate thresholds
- track fraction_disordered — high disorder invalidates discrimination claims

## server results — complete (45/45 predictions, 2026-02-23)

### overall: 5 RNA / 10 DNA / 1 BOTH / 3 WEAK / 26 NO_DOCK
- meaningful (RNA+DNA+BOTH): 16/45 = 36%
- uninformative (WEAK+NO_DOCK): 29/45 = 64%

### per-family breakdown:

| family | known | n | RNA | DNA | BOTH | WEAK | NO_DOCK | notes |
|--------|-------|---|-----|-----|------|------|---------|-------|
| DNA_pol_B (PF00136) | DNA | 6 | 3 | 2 | 1 | 0 | 0 | best family — all 6 dock. avg contacts: 9.7 RNA, 6.3 DNA |
| Piwi (PF02171) | RNA | 3 | 0 | 3 | 0 | 0 | 0 | all dock, all DNA. correct biology (piRNA→DNA) |
| KH (PF00013) | RNA | 3 | 0 | 2 | 0 | 0 | 1 | 2/3 dock, both DNA |
| DEAD_helic (PF00270) | RNA/DNA | 5 | 1 | 1 | 0 | 0 | 3 | mixed, mostly no dock |
| S4_RNA-bd (PF00575) | RNA | 4 | 0 | 1 | 0 | 0 | 3 | only bsu docks → DNA |
| DNA_topo (PF01131) | DNA | 3 | 1 | 0 | 0 | 0 | 2 | only hsa docks → RNA (surprise) |
| tRNA-synt (PF00133) | RNA | 8 | 0 | 1 | 0 | 0 | 7 | only mja docks → DNA. needs tRNA structure |
| MutS (PF00488) | DNA | 6 | 0 | 0 | 0 | 3 | 3 | zero contacts anywhere. needs mismatched dsDNA |
| UvrD_helic (PF00580) | DNA | 4 | 0 | 0 | 0 | 0 | 4 | all NO_DOCK |
| SpoU_MeTase (PF00588) | RNA | 2 | 0 | 0 | 0 | 0 | 2 | all NO_DOCK |
| OB_fold (PF01336) | DNA/RNA | 1 | 0 | 0 | 0 | 0 | 1 | NO_DOCK |

### key observations:
1. DNA Pol B is the star — only family where all species dock confidently
2. 64% of large-protein predictions are uninformative (probe too small)
3. known RNA binders that DO dock often show DNA preference (KH, Piwi, S4)
4. Piwi DNA preference is correct biology, not artifact
5. proteins needing specific substrates (tRNA-synt, MutS) almost never dock
6. local batch (≤700 aa) expected to have much higher docking rate

## files

| file | what |
|------|------|
| `structures/af3_competitive_server/all_server_competitive.json` | all 45 server jobs |
| `structures/af3_competitive_server/remaining_server_competitive.json` | 19 remaining jobs to submit |
| `structures/af3_competitive_server/unzipped/` | completed server results |
| `structures/af3_competitive_local_batch/` | 93 local job inputs |
| `structures/af3_competitive_output/` | local batch results (running) |
| `scripts/14_generate_competitive_inputs.py` | server input generator |
| `scripts/14b_generate_local_competitive_inputs.py` | local input generator |
| `scripts/16_run_competitive_batch.sh` | local batch runner |
