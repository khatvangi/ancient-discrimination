# The Evolution of Nucleic Acid Discrimination

Computational analysis showing that ancient protein domains in LUCA's
proteome exhibit widespread weak RNA/DNA discrimination in structural-contact
analyses, consistent with ancestral generalism. Binding specificity appears
to be a derived trait.

## Key findings

- **871** LUCA-age Pfam domains identified (Wehbi et al. 2024 classification)
- **54** RNA-binding, **21** DNA-binding families at LUCA age (2.6:1 ratio)
- **19** rows (18 independent families) with paired RNA/DNA structural data analyzed via NAS-Bench + ASR
  (PF07521 and PF10996 share the same source protein 9BCU)
- **72%** of ancestral NA-contacting residues differ from modern specialists
- **59%** radical substitutions — the ancestral binding mode is largely extinct

## Repository structure

```
├── THESIS.md                    thesis and claims
├── CITATION.cff                 citation metadata
├── LICENSE                      MIT license
├── .zenodo.json                 Zenodo upload metadata
├── data/
│   ├── wehbi_luca/              Wehbi et al. 2024 LUCA Pfam classification
│   ├── rbpworld/                RBPWorld RNA-binding Pfam domains
│   ├── alva_fragments/          Alva et al. 2015 ancient peptide data
│   └── sequences/               representative protein sequences
├── scripts/                     analysis scripts (49 files)
│   ├── nasbench.py              core NAS-Bench contact analysis
│   ├── asr_run_family.py        ASR + convergence pipeline per family
│   ├── substitution_analysis.py conservative/radical substitution classifier
│   ├── make_publication_figures.py   5 publication figures
│   └── make_publication_tables.py    5 publication tables
├── results/
│   ├── nasbench_full_luca.tsv   NAS-Bench SI/DI for 19 rows (18 independent families)
│   ├── phase1_intersections.json census intersection data
│   ├── tables/                  5 publication tables (TSV)
│   └── asr/                     per-family ASR results (19 rows, 18 independent families)
│       ├── convergence_master.tsv   master convergence table
│       └── {pfam_id}/              convergence, substitution, tree per family
└── figures/
    ├── fig1-fig5 (PNG + PDF)    5 publication figures
    └── structure_data/          per-residue contact TSVs for PyMOL work
```

**Note:** AF3/Boltz-2 structural predictions (Phase 3) are planned but not
yet completed. Pilot feasibility data and input-generation scripts are
included. Predicted structures, when available, will be deposited separately.

## Reproducing the analysis

### Phase 1: Census (no structural tools needed)

```bash
# compute LUCA x RNA-binding Pfam intersection
python scripts/nasbench_full_luca.py

# generate publication tables and figures
python scripts/make_publication_tables.py
python scripts/make_publication_figures.py
```

### Phase 2: ASR convergence (requires IQ-TREE, MAFFT, HMMER, CD-HIT)

```bash
# run ASR pipeline for a single family
conda activate phylo_asr
python scripts/asr_run_family.py PF00013 1EC6 1ZTG

# run substitution analysis across all families
python scripts/substitution_analysis.py
```

### Dependencies

- Python 3.9+ with: biopython, pandas, numpy, scipy, matplotlib, seaborn, requests
- conda environment `phylo_asr` with: MAFFT, IQ-TREE 2, HMMER 3, CD-HIT

## Figure / Table provenance

Every output figure and table can be regenerated from the data in `results/`.

| Output | Script | Key inputs |
|--------|--------|------------|
| Fig 1 (`fig1_validation`) | `make_publication_figures.py` | `nasbench_full_luca.tsv`, `nasbench_modern_controls.tsv` |
| Fig 2 (`fig2_ancestral_chemistry`) | `make_publication_figures.py` | `asr/pf*/substitution_analysis.tsv`, `asr/convergence_master.tsv` |
| Fig 3 (`fig3_lysine_claw`) | `make_publication_figures.py` | `asr/pf03129/PF03129_convergence_summary.tsv` |
| Fig 4 (`fig4_cross_tool`) | `make_publication_figures.py` | `nasbench_full_luca.tsv`, `prona2020_luca_domains.tsv`, `drbp_edp_luca_domains.tsv` |
| Fig 5 (`fig5_kh_phylogeny`) | `make_publication_figures.py` | `asr/kh/kh_prona_all_nodes.tsv` |
| Table 1 (`table1_census`) | `make_publication_tables.py` | `phase1_intersections.json`, `phase1_census.tsv`, `luca_pdb_census.tsv`, `nasbench_full_luca.tsv`, `asr/convergence_master.tsv` |
| Table 2 (`table2_nasbench_asr`) | `make_publication_tables.py` | `nasbench_full_luca.tsv`, `prona2020_luca_domains.tsv`, `drbp_edp_luca_domains.tsv`, `asr/convergence_master.tsv` |
| Table 3 (`table3_substitutions`) | `make_publication_tables.py` | `asr/pf*/substitution_analysis.tsv`, `asr/convergence_master.tsv` |
| Table 4 (`table4_cross_tool`) | `make_publication_tables.py` | `nasbench_full_luca.tsv`, `prona2020_luca_domains.tsv`, `drbp_edp_luca_domains.tsv` |
| Table 5 (`table5_modern_controls`) | `make_publication_tables.py` | `nasbench_modern_controls.tsv` |

All input paths are relative to `results/`. Figures are written to `figures/` as PNG + PDF.

## Key references

- Wehbi et al. 2024 PNAS — 871 LUCA-age Pfam domain classification
- Liao et al. 2025 NAR — RBPWorld database (998 RNA-binding Pfam domains)
- Alva et al. 2015 eLife — 40 ancient peptide vocabulary fragments
- Moody et al. 2024 Nature Ecol Evol — 2,600 LUCA proteins

## Citation

See `CITATION.cff` for citation information.

## License

MIT — see `LICENSE`.
