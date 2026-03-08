# The Evolution of Nucleic Acid Discrimination

Computational analysis demonstrating that ancient protein domains in LUCA's
proteome were generalist nucleic acid binders that did not distinguish between
RNA and DNA. Binding specificity is a derived trait.

## Key findings

- **871** LUCA-age Pfam domains identified (Wehbi et al. 2024 classification)
- **54** RNA-binding, **21** DNA-binding families at LUCA age (2.6:1 ratio)
- **19** families with paired RNA/DNA structural data analyzed via NAS-Bench + ASR
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
│   ├── nasbench_full_luca.tsv   NAS-Bench SI/DI for 19 paired families
│   ├── phase1_intersections.json census intersection data
│   ├── tables/                  5 publication tables (TSV)
│   └── asr/                     per-family ASR results (19 families)
│       ├── convergence_master.tsv   master convergence table
│       └── {pfam_id}/              convergence, substitution, tree per family
└── figures/
    ├── fig1-fig5 (PNG + PDF)    5 publication figures
    └── structure_data/          per-residue contact TSVs for PyMOL work
```

**Note:** AF3/Boltz-2 predicted structures (14 GB) are not included in this
repository. They are available upon request from the corresponding author.

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

## Key references

- Wehbi et al. 2024 PNAS — 871 LUCA-age Pfam domain classification
- Liao et al. 2025 NAR — RBPWorld database (998 RNA-binding Pfam domains)
- Alva et al. 2015 eLife — 40 ancient peptide vocabulary fragments
- Moody et al. 2024 Nature Ecol Evol — 2,600 LUCA proteins

## Citation

See `CITATION.cff` for citation information.

## License

MIT — see `LICENSE`.
