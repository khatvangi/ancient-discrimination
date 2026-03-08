# ancient-discrimination/

Dedicated workspace for the project: **"The Evolution of Nucleic Acid Discrimination"**

## Thesis
Ancient protein domains did not distinguish between RNA and DNA. Nucleic acid
binding specificity is a derived trait that emerged through domain diversification
after LUCA.

## Structure
```
ancient-discrimination/
├── THESIS.md              this project's thesis and claims
├── README.md              this file
├── data/                  downloaded reference datasets
│   ├── wehbi_luca/        Wehbi et al. 2024 LUCA Pfam classification
│   ├── rbpworld/          RBPWorld RNA-binding Pfam domains
│   └── alva_fragments/    Alva et al. 2015 ancient peptide data
├── scripts/               analysis scripts (self-contained)
├── results/               computed outputs (TSV, JSON)
├── structures/            AF3/Boltz-2 predicted complexes
└── figures/               publication figures
```

## Key documents
- `../docs/plans/2026-02-21-ancient-discrimination-design.md` — full design
- `../GAP_ANALYSIS_2026-02-21.md` — literature gap analysis
- `../PROJECT_THESIS_NOTES_2026-02-21.md` — original thesis notes

## Quick start
Phase 1 (census) can be run without any structural prediction tools.
See design doc for full phase breakdown.
