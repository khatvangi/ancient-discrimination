# Phase 1 Summary — Census of LUCA's NA-binding Proteome

## Headline Number

**54 of 871 LUCA Pfam domains (6.2%) are RNA-binding.**

This number did not exist in the literature prior to this computation.

## Key Results

| Metric | Value |
|--------|-------|
| LUCA RNA-binding domains | 54 (6.2% of 871) |
| pre-LUCA RNA-binding domains | 8 (7.8% of 102) |
| LUCA DNA-binding domains | 21 (2.4% of 871) |
| LUCA dual-binding (annotated) | 0 |
| RNA:DNA ratio in LUCA | 2.6:1 |
| Annotation coverage ratio (RNA:DNA sources) | 2.4:1 |

## Functional Profiles

**RNA-binding (54 domains):** translation-dominated
- translation: 20 (37%) — tRNA synthetases, ribosomal proteins, SRP, EF-Tu
- rna_processing: 9 (17%) — DEAD-box helicases, PIN ribonucleases, RNA ligase
- rna_modification: 7 (13%) — methyltransferases, pseudouridine synthases
- rna_binding_structural: 7 (13%) — KH, S1, S4, PUA domains
- other: 11 (20%) — moonlighters, uncharacterized

**DNA-binding (21 domains):** repair-dominated
- dna_repair: 9 (43%) — MutS, integrases, ERCC4, HHH
- dna_modification: 4 (19%) — methyltransferases, restriction enzymes
- transcription: 3 (14%) — HTH variants
- dna_topology: 2 (10%) — topoisomerase, gyrase
- dna_replication: 1 (5%) — DNA polymerase B
- other: 2 (10%)

## Alva Ancient Vocabulary

25 of 83 LUCA NA-binding Pfam domains (30%) contain structural elements from
Alva et al. 2015's 40 oldest known protein motifs.

Most thesis-relevant fragments:
- Fragment 2 (HhH): appears in both RNA and DNA methyltransferases at LUCA age
- Fragment 3 (P-loop): maps to 9 LUCA domains (7 RNA + 2 DNA)
- Fragment 7 (Rossmann): annotated DNA-binding at fold level, but LUCA function is tRNA synthetase (RNA)

## Kill Switches

| Gate | Result | Decision |
|------|--------|----------|
| LUCA × RNA > 30 | 54 | PASS — proceed |
| GREEN-51 × LUCA > 10 | 0 | FAIL — pivot to Wehbi LUCA domains |

## Caveats

1. The 2.6:1 RNA:DNA ratio partly reflects annotation asymmetry (916 RNA-binding
   Pfams from RBPWorld vs 385 DNA-binding from GO only). No curated DNA-binding
   equivalent of RBPWorld exists.
2. Zero dual-binders is annotation artifact — domains are classified by modern
   function, not ancestral capability.
3. Pre-LUCA 7.8% vs LUCA 6.2% is not statistically significant at these sample sizes.

## Cross-validation

| Reference | Their estimate | Ours | Status |
|-----------|---------------|------|--------|
| Anantharaman et al. 2002 | ~40-45 RNA metabolism in LUCA | 54 RNA-binding | consistent |
| Crapitto et al. 2022 | translation-dominated | translation 37% | consistent |
| Wehbi et al. 2024 | 969/871 LUCA Pfams | 871 (robust) | using robust set |

## Files Produced

| File | Description |
|------|-------------|
| `data/wehbi_luca/luca_pfams.tsv` | 8,288 Pfam age classifications |
| `data/rbpworld/rna_binding_combined.tsv` | 968 three-tier RNA-binding reference |
| `data/rbpworld/dna_binding_pfams.tsv` | 385 DNA-binding Pfam domains |
| `data/alva_fragments/fragments.tsv` | Alva 2015 40 ancient peptides |
| `data/ipr_to_pfam_map.tsv` | 253 IPR → Pfam mappings |
| `results/phase1_intersections.json` | all intersection counts and IDs |
| `results/phase1_census.tsv` | 9,065-row per-Pfam annotation table |
| `results/phase1_functional_composition.tsv` | 83 classified LUCA NA-binding domains |
| `results/phase1_alva_luca_mapping.tsv` | fragment → LUCA domain mapping |
| `results/phase1_journal.md` | detailed methodology journal |
