# ──────────────────────────────────────────────────────────────────────────────
# Makefile for "The Evolution of Nucleic Acid Discrimination"
#
# Pipeline overview:
#   1. census   — Phase 1: download datasets, compute LUCA x RNA/DNA intersection
#   2. nasbench — NAS-Bench contact analysis for LUCA and modern control families
#   3. asr      — ancestral state reconstruction + substitution analysis (19 families)
#   4. figures  — generate 5 publication figures (PNG + PDF)
#   5. tables   — generate 5 publication tables (TSV)
#
# Prerequisites:
#   conda activate phylo_asr   (for ASR target — needs MAFFT, IQ-TREE, HMMER, CD-HIT)
#   pip install -r requirements.txt (for all other targets)
#
# Usage:
#   make census       # run Phase 1 intersection
#   make nasbench     # run NAS-Bench analysis
#   make asr          # run ASR for all 18 families (hours; needs phylo_asr env)
#   make figures      # generate all 5 publication figures
#   make tables       # generate all 5 publication tables
#   make all          # run everything in order
# ──────────────────────────────────────────────────────────────────────────────

PYTHON := python3
SCRIPTS := scripts
RESULTS := results
FIGURES := figures

.PHONY: all census nasbench asr figures tables clean

all: census nasbench asr figures tables

# ── Phase 1: Census ──────────────────────────────────────────────────────────
# downloads Wehbi LUCA domains, RBPWorld RNA-binding families, Alva fragments;
# maps InterPro to Pfam; computes intersections and functional classification.

census: $(RESULTS)/phase1_intersections.json

$(RESULTS)/phase1_intersections.json:
	$(PYTHON) $(SCRIPTS)/01_download_wehbi.py
	$(PYTHON) $(SCRIPTS)/02_download_rbpworld.py
	$(PYTHON) $(SCRIPTS)/03_download_alva.py
	$(PYTHON) $(SCRIPTS)/04_map_ipr_to_pfam.py
	$(PYTHON) $(SCRIPTS)/05_compute_intersections.py
	$(PYTHON) $(SCRIPTS)/06_functional_classification.py
	$(PYTHON) $(SCRIPTS)/07_alva_luca_mapping.py

# ── NAS-Bench analysis ───────────────────────────────────────────────────────
# computes Specificity Index (SI) and Discrimination Index (DI) from PDB
# co-crystal structures for modern controls and LUCA-age families.

nasbench: $(RESULTS)/nasbench_full_luca.tsv $(RESULTS)/nasbench_modern_controls.tsv

$(RESULTS)/nasbench_full_luca.tsv:
	$(PYTHON) $(SCRIPTS)/nasbench_full_luca.py

$(RESULTS)/nasbench_modern_controls.tsv:
	$(PYTHON) $(SCRIPTS)/nasbench_modern_controls.py

# ── ASR (ancestral state reconstruction) ─────────────────────────────────────
# runs MAFFT + IQ-TREE + ASR + convergence analysis for all 18 families.
# requires conda activate phylo_asr (MAFFT 7.525, IQ-TREE 3.0.1, HMMER 3.4, CD-HIT 4.8.1).
# this target takes several hours on 16 cores.

asr: $(RESULTS)/asr/convergence_master.tsv

$(RESULTS)/asr/convergence_master.tsv:
	bash $(SCRIPTS)/batch_asr_all.sh
	$(PYTHON) $(SCRIPTS)/substitution_analysis.py

# ── Figures ──────────────────────────────────────────────────────────────────
# generates 5 publication figures (PNG + PDF) from precomputed results.
#   fig1_validation       — NAS-Bench modern vs LUCA SI/DI
#   fig2_ancestral_chemistry — ancestral contact substitution patterns
#   fig3_lysine_claw      — PF03129 case study
#   fig4_cross_tool       — cross-tool disagreement heatmap
#   fig5_kh_phylogeny     — KH domain ProNA2020 P_RNA vs P_DNA

figures: $(FIGURES)/fig1_validation.pdf

$(FIGURES)/fig1_validation.pdf: $(RESULTS)/nasbench_full_luca.tsv $(RESULTS)/nasbench_modern_controls.tsv $(RESULTS)/asr/convergence_master.tsv
	$(PYTHON) $(SCRIPTS)/make_publication_figures.py

# ── Tables ───────────────────────────────────────────────────────────────────
# generates 5 publication tables (TSV) from precomputed results.
#   table1_census         — LUCA NA-binding domain census summary
#   table2_nasbench_asr   — SI/DI + ASR root states for 18 families
#   table3_substitutions  — substitution analysis by contact type
#   table4_cross_tool     — cross-tool comparison matrix
#   table5_modern_controls — modern specialist validation

tables: $(RESULTS)/tables/table1_census.tsv

$(RESULTS)/tables/table1_census.tsv: $(RESULTS)/nasbench_full_luca.tsv $(RESULTS)/asr/convergence_master.tsv
	$(PYTHON) $(SCRIPTS)/make_publication_tables.py

# ── Clean ────────────────────────────────────────────────────────────────────

clean:
	rm -f $(FIGURES)/fig*.png $(FIGURES)/fig*.pdf
	rm -f $(RESULTS)/tables/table*.tsv
