# Project Journal — Ancient Discrimination

**Dates**: 2026-02-21 to present
**Status**: Phase 3 in progress (AF3 structural predictions)

---

## Headline Numbers

| Metric | GO (lower bound) | RBPWorld (primary) | Union (upper) |
|--------|------------------|--------------------|---------------|
| LUCA × RNA-binding | 30 (3.4%) | **54 (6.2%)** | 54 (6.2%) |
| pre-LUCA × RNA-binding | 6 (5.9%) | **8 (7.8%)** | 8 (7.8%) |
| ancient × RNA-binding | 36 (3.7%) | **62 (6.4%)** | 62 (6.4%) |
| LUCA × DNA-binding | 21 (2.4%) | 21 (2.4%) | 21 (2.4%) |
| LUCA × dual (RNA+DNA) | 0 | 0 | 0 |

**Primary estimate: 54 of 871 LUCA Pfam domains (6.2%) are RNA-binding.**
**8 of 102 pre-LUCA Pfam domains (7.8%) are RNA-binding — higher fraction.**

---

## Data Sources Used

### Wehbi et al. 2024 LUCA Pfam classification
- Source: github.com/sawsanwehbi/Pfam-age-classification
- File: ClassifiedPFAMs.csv (robust classifications)
- Result: 871 LUCA + 102 pre-LUCA + 8,288 total classified Pfams
- **Note**: Paper reports 969 LUCA; we use 871 (robust classification excludes
  98 domains reclassified as "unclassifiable" by the authors post-publication)

### RNA-binding Pfam domains (three tiers)
- **GO-derived (487)**: Pfam2GO mapping + QuickGO hierarchy (164 descendant
  terms of GO:0003723). Conservative lower bound.
- **RBPWorld/EuRBPDB (916)**: Combined from RBPWorld API (849 families) +
  EuRBPDB HMM archive (791 domains with ACC lines). Experimentally curated.
- **Union (968)**: All unique Pfam IDs from both sources.
- Note: RBPWorld claims 998 RBDs but this counts HMM profiles; only 849 had
  protein hits. Combined with EuRBPDB, we get 916 unique Pfam IDs. Union
  with GO-based gives 968.

### DNA-binding Pfam domains
- Source: Pfam2GO + QuickGO (131 descendants of GO:0003677)
- Result: 385 DNA-binding Pfam domains

### GREEN-51 families
- Source: InterPro API mapping of 51 IPR IDs to Pfam
- Result: only 25/51 had Pfam member databases; the rest use CDD/SMART/other
- **KILL SWITCH TRIGGERED**: 0 of 25 Pfam-mapped GREEN families are LUCA-age
  (most are LBCA = bacterial origin)

---

## Key Findings

### 1. LUCA's RNA-binding proteome is translation-dominated
The 54 LUCA RNA-binding domains are overwhelmingly:
- tRNA synthetases (13+ domains across both classes)
- Ribosomal proteins (L4, L6, L10, S4)
- RNA modification enzymes (methyltransferases, pseudouridine synthases)
- RNA helicases (DEAD-box, UvrD)
- Signal recognition particle (SRP54)

### 2. Pre-LUCA is even more RNA-centric
The 8 pre-LUCA RNA-binding domains include EF-Tu (GTP_EFTU) and core tRNA
synthetase domains — the absolute oldest translation machinery. The pre-LUCA
RNA-binding fraction (7.8%) exceeds LUCA (6.2%), consistent with "older = more
RNA-associated."

### 3. LUCA's DNA-binding proteome is repair-dominated
The 21 LUCA DNA-binding domains are:
- DNA repair: MutS, HHH, ERCC4
- DNA modification: methyltransferases (N6_N4_Mtase, Methylase_S)
- DNA topology: topoisomerase, gyrase
- DNA replication: DNA_pol_B
- Transcription factors: HTH variants

### 4. Zero annotated dual-binders
No LUCA domain is annotated as binding BOTH RNA and DNA. This likely reflects
modern functional annotation bias — domains are classified based on their
current function, not ancestral capability. This is exactly the gap our
structural prediction phases (3-5) are designed to address.

### 5. GREEN-51 kill switch triggered
Our 51 GREEN families are NOT LUCA-age. They are predominantly LBCA (bacterial
origin) that spread to archaea and eukaryotes via horizontal gene transfer.
**Universality ≠ LUCA-age.** This means Phase 2 must select from the Wehbi
54 LUCA RNA-binding domains, not our GREEN families.

### 6. Functional classification (Task 6-7)

**54 LUCA RNA-binding domains:**

| Category | Count | % |
|----------|-------|---|
| translation | 20 | 37.0% |
| other | 11 | 20.4% |
| rna_processing | 9 | 16.7% |
| rna_modification | 7 | 13.0% |
| rna_binding_structural | 7 | 13.0% |

**21 LUCA DNA-binding domains:**

| Category | Count | % |
|----------|-------|---|
| dna_repair | 9 | 42.9% |
| dna_modification | 4 | 19.0% |
| transcription | 3 | 14.3% |
| dna_topology | 2 | 9.5% |
| dna_replication | 1 | 4.8% |
| other | 1 | 4.8% |
| rna_processing | 1 | 4.8% |

**Key contrast:** RNA = translation-dominated (37%), DNA = repair-dominated (43%).
LUCA's RNA proteome served the ribosome; its DNA proteome served genome maintenance.

### 7. Alva ancient vocabulary mapping (Task 9)

9 of 40 Alva fragments (22%) map to at least one LUCA NA-binding Pfam domain.
25 of 83 LUCA NA-binding Pfams (30%) are structurally related to Alva's ancient
peptide vocabulary.

**Notable cross-matches:**
- Fragment 2 (HhH motif, Weil-Ktorza "ambidextrous"): maps to 6 LUCA domains
  (5 RNA methyltransferases + 1 DNA methyltransferase). Both RNA and DNA
  modification enzymes share this ancient structural element.
- Fragment 3 (P-loop NTPase): maps to 9 LUCA domains (7 RNA + 2 DNA). The most
  prolific ancient motif in LUCA's NA-binding proteome.
- Fragment 7 (Rossmann fold): maps to 5 class I aaRS domains (all RNA-binding).
  Alva annotated this as "DNA-binding" based on fold-level capability, but in LUCA
  the actual function was tRNA synthetase — an RNA interaction.

**Coverage by function:**
- rna_modification: 50% (4/8 domains contain ancient vocabulary)
- translation: 33% (8/24)
- transcription: 100% (3/3 — all HTH variants)
- dna_repair: 11% (1/9)

### 8. Annotation bias check

DNA-binding Pfam annotations (385) come only from GO:0003677 descendants.
RNA-binding annotations (916) incorporate experimental evidence (CLIP, RIC, OOPS).
The 2.4:1 ratio of annotation sources means the observed 2.6:1 LUCA RNA:DNA ratio
may partly reflect annotation asymmetry. However, the functional profiles are
structurally different (translation vs repair) and robust to coverage gaps.

---

## Cross-validation

| Reference | Their estimate | Our estimate | Status |
|-----------|---------------|--------------|--------|
| Anantharaman et al. 2002 | ~40-45 RNA metabolism domains in LUCA | 54 RNA-binding | Consistent (ours slightly higher, includes structural domains) |
| Crapitto et al. 2022 | 111 GO terms, translation-dominated | Translation-dominated | Consistent |
| Wehbi et al. 2024 | 969 (paper) / 871 (robust) LUCA Pfams | 871 used | Using robust classification |

---

## Kill Switch Evaluation

| Gate | Threshold | Result | Decision |
|------|-----------|--------|----------|
| LUCA × RNA count | > 30 | 54 | **PASS** — proceed |
| GREEN-51 × LUCA | > 10 | 0 | **FAIL** — PIVOT to Wehbi LUCA domains |

---

## Methodological Notes

### Pfam2GO coverage gap
Only ~5,000 of ~20,000+ Pfam entries have GO annotations in the pfam2go
mapping. This systematically undercounts RNA-binding for the GO-only tier.
The RBPWorld/EuRBPDB tier compensates by using experimentally curated domain
lists (CLIP, RIC, OOPS evidence).

### InterPro → Pfam mapping loss
26 of 51 GREEN InterPro entries (51%) have no Pfam member databases. This is
because InterPro integrates multiple source databases; some entries are only
in CDD, SMART, PROSITE, or are directly curated at InterPro level. This is
a systematic limitation when bridging between InterPro and Pfam-based analyses.

### Pre-LUCA deserves special attention
The 102 pre-LUCA domains (8 RNA-binding, 7.8%) are the strongest candidates
for testing the "ancient generalism" hypothesis. These predate LUCA and are
likely the most ancient protein domains still recognizable. They should be
prioritized in Phase 2 family selection.

---

## Files Produced

| File | Description | Rows |
|------|-------------|------|
| `data/wehbi_luca/luca_pfams.tsv` | All Pfam age classifications | 8,288 |
| `data/wehbi_luca/raw/ClassifiedPFAMs.csv` | Original download | — |
| `data/rbpworld/rna_binding_combined.tsv` | Three-tier RNA-binding reference | 968 |
| `data/rbpworld/rbd_pfam_ids.txt` | RBPWorld/EuRBPDB Pfam IDs | 916 |
| `data/rbpworld/dna_binding_pfams.tsv` | DNA-binding Pfam domains | 385 |
| `data/rbpworld/na_binding_pfams.tsv` | All NA-binding (union) | 855 |
| `data/alva_fragments/fragments.tsv` | Alva 2015 40 ancient peptides | 40 |
| `data/ipr_to_pfam_map.tsv` | IPR → Pfam mappings | 253 |
| `results/phase1_intersections.json` | All intersection counts and IDs | — |
| `results/phase1_census.tsv` | Per-Pfam annotation table | 9,065 |
| `results/phase1_functional_composition.tsv` | Functional classification of 83 domains | 83 |
| `results/phase1_alva_luca_mapping.tsv` | Alva fragment → LUCA mapping | 40 |
| `results/phase1_alva_summary.json` | Alva mapping summary statistics | — |

---

---
---

# Phase 2 — Family Selection and Ortholog Identification

**Date**: 2026-02-21
**Status**: Complete (Tasks 8-10)

## Decisions

- selected 21 families across 3 tiers (8 dual-binders, 7 controls, 6 thesis-critical)
- used PDB structure survey to identify families with BOTH RNA and DNA complex structures
- species panel: 9 organisms (3 Bac, 3 Arc, 3 Euk), 148 orthologs total
- prediction manifest: 296 jobs (148 proteins × RNA + DNA)

## Key outcomes

- 19/21 families pass ≥5 orthologs / ≥2 domains
- Piwi and Phage_integrase below threshold but retained (thesis-critical / control)
- H. volcanii is sparse in Swiss-Prot — missing from 15/21 families
- 91 orthologs (61%) have existing PDB structures for validation

## Files produced

see `results/phase2_summary.md` for full detail.

---
---

# Phase 3 — AF3 Structural Predictions

**Date**: 2026-02-22
**Status**: In progress

## Session log

### 1. AF3 input generation (Task 11)

generated 296 AF3 JSON input files in `structures/af3_inputs/`.

**AF3 native format** (for local runs):
- `dialect: "alphafold3"`, `version: 1`
- chain keys: `protein`, `rna`, `dna` with `id` fields ("A", "B")
- substrates: poly-U 10mer (RNA), poly-dT 10mer (DNA)

**AF3 Server format** (for large proteins):
- `dialect: "alphafoldserver"`, `version: 1`
- chain keys: `proteinChain`, `rnaSequence`, `dnaSequence` with `count` field
- top-level is a list of jobs

### 2. pilot v1 — bugs encountered and fixed

launched 6 pilot jobs on 2 GPUs (2× Titan RTX 24GB). three bugs surfaced:

**bug 1: AF3 JSON format mismatch**
- used AlphaFold Server format (`proteinChain`/`rnaSequence`/`dnaSequence`) instead of native format
- fix: regenerated inputs with `dialect: "alphafold3"`, `protein`/`rna`/`dna` keys + `id` fields

**bug 2: mmseqs not on PATH**
- `mmseqs` binary at `/home/kiran/miniforge3/bin/mmseqs` not in af3_mmseqs2 conda env
- fix: `export PATH="/home/kiran/miniforge3/bin:$PATH"` in run scripts

**bug 3: absl module not found**
- `conda activate af3_mmseqs2` doesn't work in bash subshells
- fix: use direct python path `/home/kiran/miniforge3/envs/af3_mmseqs2/bin/python`

### 3. pilot v1 — mmseqs createdb race condition

both GPU processes started `mmseqs createdb` on UniRef90 simultaneously (first-time
database build). classic TOCTOU race: both checked if DB existed → neither saw the
other → both wrote to same output directory. GPU 0's files got overwritten by GPU 1's
writes (visible as `(deleted)` file handles in /proc). GPU 0 crashed, GPU 1 survived.

**one-time cost**: mmseqs databases created for all 4 search databases:
- UniRef90 (67 GB FASTA → 73 GB index, ~25 min)
- MGnify (120 GB → built by GPU 1 after UniRef90)
- small_BFD (17 GB)
- UniProt (102 GB)

total first-time database creation: ~4 hours. subsequent runs skip this entirely.

### 4. pilot v1 — hmmbuild not found

GPU 1 completed all 4 MSA searches (UniRef90, MGnify, small_BFD, UniProt) successfully,
then crashed at the template search step:
```
RuntimeError: hmmbuild binary not found at hmmbuild
```

root cause: PATH included `/home/kiran/miniforge3/bin` (has `mmseqs`) but NOT
`/home/kiran/miniforge3/envs/af3_mmseqs2/bin` (has `hmmbuild`, `hmmsearch`, `hmmalign`,
`nhmmer`, `jackhmmer`).

GPU 0 rerun (after race condition fix) would have hit the same error.

### 5. pilot v2 — PATH fixed, relaunched

created `scripts/12_run_af3_pilot_v2.sh` with corrected PATH:
```bash
export PATH="/home/kiran/miniforge3/envs/af3_mmseqs2/bin:/home/kiran/miniforge3/bin:$PATH"
```

also removed `set -e` (replaced with `set -uo pipefail`) so individual job failures
don't kill the entire GPU batch.

verified all binaries found (mmseqs, hmmbuild, hmmsearch, hmmalign) before launching.
pilot v2 started 2026-02-22 09:11 CST. both GPUs allocated (23 GB each).

as of last check (~10:00 CST):
- GPU 0 (S4-RNA): searching UniProt (4th/last database)
- GPU 1 (HHH-DNA): searching small_BFD (3rd/4th database)

MSA searches running ~8-10 min per database (vs hours in v1 when building DBs).
estimated first results: ~2-3 hours after launch.

### 6. GPU memory analysis — large protein strategy

analyzed 148 proteins against 24GB GPU memory limit:

| category | count | tokens | GPU needed |
|----------|-------|--------|------------|
| fits easily | 70 | ≤512 | 24 GB |
| fits tight | 49 | 513-1024 | 24 GB (test) |
| risky | 16 | 1025-1536 | may OOM |
| too large | 3 | >1536 | needs A100 |

20 proteins >1000 aa are problematic on 24GB. largest: DNA_pol_B from A. thaliana (1890 aa = 1900 tokens).

### 7. AF3 Server for large proteins

generated 40 AF3 Server JSON inputs (`structures/af3server_inputs/`) for the 20 proteins
>1000 aa. all under the 5000-token server limit (largest: 1900 tokens).

AF3 Server: alphafoldserver.com, 30 jobs/day, non-commercial, Google account.
40 jobs = 2 days (30 today + 10 tomorrow). user submitted first 30 batch.

**parallel execution strategy:**
- local 24GB GPUs: pilot (6 jobs) → easy batch (140 jobs) → tight batch (98 jobs, test for OOM)
- AF3 Server: 40 jobs for proteins >1000 aa (submitted 2026-02-22)
- RunPod A100 (backup): if tight batch has many OOM failures

### 8. split pipeline capability confirmed

AF3 supports `--run_data_pipeline=true --run_inference=false` and vice versa.
this means MSA can run locally (CPU-only, no GPU memory limit) and inference can
be shipped to RunPod A100. however, AF3 Server doesn't accept pre-computed MSA —
it only takes raw sequences. so the split is only useful for RunPod, not for the Server.

## pilot jobs (awaiting results)

| GPU | job | protein | size | annotated | thesis prediction |
|-----|-----|---------|------|-----------|-------------------|
| 0 | PF01479_bsu_P37557_RNA | S4 ribosomal | 86 aa | RNA | specialist (RNA) |
| 0 | PF01479_bsu_P37557_DNA | S4 ribosomal | 86 aa | RNA | low DNA ipTM |
| 0 | PF00633_eco_P0AB83_RNA | HHH motif | 211 aa | DNA | **key test**: generalist? |
| 1 | PF00633_eco_P0AB83_DNA | HHH motif | 211 aa | DNA | should bind DNA well |
| 1 | PF00575_mja_Q57840_RNA | S1 OB-fold | 187 aa | RNA | known dual-binder |
| 1 | PF00575_mja_Q57840_DNA | S1 OB-fold | 187 aa | RNA | may also bind DNA |

**what to look for:**
- S4: expect DI > 0.1 (RNA specialist). positive control.
- HHH: if DI ≈ 0 → supports thesis (ancient generalist). **key test.**
- S1: known dual-binder, expect DI near 0 or mildly positive. positive control.
- all ipTM > 0.3 for sanity (if lower, substrate too short or predictions unreliable)

## scripts produced in phase 3

| script | purpose |
|--------|---------|
| `scripts/11_generate_af3_inputs.py` | generate local AF3 JSON inputs (296 files) |
| `scripts/11b_generate_af3server_inputs.py` | generate AF3 Server JSON inputs (40 files) |
| `scripts/12_run_af3_pilot.sh` | pilot v1 (had PATH bugs) |
| `scripts/12_run_af3_pilot_v2.sh` | pilot v2 (fixed, currently running) |
| `scripts/12_rerun_gpu0.sh` | GPU 0 rerun (superseded by v2) |
| `scripts/13_analyze_af3_pilot.py` | pilot analysis: DI + go/no-go decision |

## next steps (after pilot completes)

1. run `python scripts/13_analyze_af3_pilot.py` → go/no-go on full batch
2. if GO: launch local batch for ≤1000 aa proteins (~256 jobs on 2 GPUs)
3. try "tight" batch (513-1024 tokens) on local GPUs → check for OOM
4. collect AF3 Server results for large proteins (2 days)
5. if tight batch has many OOM: use RunPod A100 with split pipeline
6. Task 14: quality filter + interface extraction once all predictions complete
