# Phase 2 Summary — Family Selection and Ortholog Identification
# NOTE: this phase selected families and orthologs for Phase 3 AF3/Boltz-2
# structural predictions. Phase 3 is planned but not yet completed.
# the prediction manifest below describes intended work, not completed runs.

## Selected Families

21 domain families selected for AF3/Boltz-2 structural prediction, organized in 3 tiers.

### Tier 1: Experimentally dual (RNA + DNA complex structures in PDB)

These directly test the thesis — domains annotated as single-substrate specialists
that already have experimental evidence of binding both RNA and DNA.

| Pfam | Name | RNA structs | DNA structs | Function | Age |
|------|------|-------------|-------------|----------|-----|
| PF00575 | S1 (OB-fold) | 308 | 257 | rna_binding_structural | LUCA |
| PF00270 | DEAD | 284 | 56 | rna_processing | LUCA |
| PF02171 | Piwi | 75 | 11 | other | LUCA |
| PF00013 | KH_1 | 23 | 10 | rna_binding_structural | LUCA |
| PF01336 | tRNA_bind | 14 | 11 | translation | preLUCA |
| PF00398 | RrnaAD | 27 | 15 | rna_modification | LUCA |
| PF01479 | S4 | 1689 | 68 | rna_binding_structural | LUCA |
| PF01131 | Topoisom_bac | 8 | 27 | dna_topology | LUCA |

### Tier 2: Controls (RNA-only and DNA-only specialists)

| Pfam | Name | RNA structs | DNA structs | Role |
|------|------|-------------|-------------|------|
| PF01472 | PUA | 29 | 0 | RNA-only control |
| PF01509 | TruB_N | 30 | 0 | RNA-only control |
| PF00448 | SRP54 | 17 | 0 | RNA-only control |
| PF00588 | SpoU_methylase | 7 | 0 | RNA-only control |
| PF00580 | UvrD-helicase | 0 | 32 | DNA-only control |
| PF00589 | Phage_integrase | 0 | 33 | DNA-only control |
| PF00136 | DNA_pol_B | 9 | 196 | DNA-specialist control |

### Tier 3: Thesis-critical (sparse PDB but important for hypothesis)

| Pfam | Name | Why critical |
|------|------|-------------|
| PF00633 | HHH | Weil-Ktorza 2025 "ambidextrous" motif; Alva fragment 2 |
| PF00009 | GTP_EFTU | Pre-LUCA; EF-Tu core translation; Alva fragment 3 |
| PF00488 | MutS_V | DNA repair with 1 anomalous RNA structure |
| PF00133 | tRNA-synt_1 | Pre-LUCA class I aaRS; Rossmann fold (Alva fragment 7) |
| PF00152 | tRNA-synt_2 | Class II aaRS core; translation |
| PF00347 | Ribosomal_L6 | Ribosomal protein; high structure coverage |

## Ortholog Panel

- **148 orthologs** across 9 species (3 Bacteria, 3 Archaea, 3 Eukarya)
- **19 of 21 families** pass the ≥5 orthologs / ≥2 domains threshold
- **91 orthologs** (61%) have existing PDB structures for validation
- **0 download failures** from UniProt

### Species panel

| Domain | Species | Short | Taxon ID |
|--------|---------|-------|----------|
| Bacteria | E. coli K-12 | eco | 83333 |
| Bacteria | B. subtilis 168 | bsu | 224308 |
| Bacteria | T. thermophilus | tth | 274 |
| Archaea | M. jannaschii | mja | 243232 |
| Archaea | H. volcanii | hvo | 2246 |
| Archaea | S. solfataricus P2 | sso | 273057 |
| Eukarya | S. cerevisiae S288C | sce | 559292 |
| Eukarya | H. sapiens | hsa | 9606 |
| Eukarya | A. thaliana | ath | 3702 |

### Families that failed threshold

- **PF02171 (Piwi)**: 3 orthologs only (1 Arc + 2 Euk). No bacterial Argonaute in Swiss-Prot. Retained because thesis-critical (RNA-guided DNA cleavage).
- **PF00589 (Phage_integrase)**: 4 orthologs (2 Bac + 2 Arc). No eukaryotic integrases. Retained as DNA-only control.

### Coverage gaps

H. volcanii (hvo) is missing from 15 of 21 families — sparse Swiss-Prot annotation for this species. All other species have good coverage.

## Prediction Manifest

- **296 total jobs**: 148 protein × RNA + 148 protein × DNA
- Protein sizes: 86–1890 aa (median 497 aa)
- 20 proteins >1000 aa — may need domain extraction for Boltz-2 efficiency
- Substrates: poly-U 10-mer (RNA), poly-dT 10-mer (DNA)

## Experimental Design

For each ortholog, AF3/Boltz-2 will predict:
1. Protein + poly-U RNA → ipTM score, interface contacts
2. Protein + poly-dT DNA → ipTM score, interface contacts

The **Discrimination Index** = ipTM(annotated substrate) - ipTM(other substrate).

If ancient domains were generalist NA binders:
- Tier 1 domains should show DI ≈ 0 (equal affinity for RNA and DNA)
- Tier 2 RNA controls should show DI > 0 (prefer RNA)
- Tier 2 DNA controls should show DI < 0 (prefer DNA)
- Cross-domain comparison: archaeal/bacterial DI should be closer to 0 than eukaryotic DI

## Files Produced

| File | Description | Rows |
|------|-------------|------|
| `results/phase2_selected_families.tsv` | 21 families with rationale | 21 |
| `results/phase2_orthologs.tsv` | ortholog panel | 148 |
| `results/phase2_prediction_manifest.tsv` | prediction job manifest | 296 |
| `results/pdb_na_complex_survey.tsv` | PDB structure survey | 36 |
| `data/sequences/*.fasta` | protein sequences | 148 |
| `data/substrates.fasta` | RNA and DNA substrates | 2 |

## Next Steps

1. **Phase 3, Task 11**: Install and benchmark Boltz-2 on known PDB complexes
2. **Phase 3, Task 12-13**: Run the 296 predictions (needs GPU)
3. Domain extraction for 20 proteins >1000 aa before prediction
