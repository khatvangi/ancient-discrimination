# Literature Summary for Ancient Discrimination Project

Condensed from `../GAP_ANALYSIS_2026-02-21.md` — read that file for full citations.

---

## THE HEADLINE GAP

Nobody has computed: **what fraction of LUCA's ~969 protein domains were RNA-binding?**

Two datasets exist with the same ID system (Pfam). Their intersection is trivial
but unpublished:
- Wehbi et al. 2024 (PNAS): 969 LUCA Pfam domains, 101 pre-LUCA
- RBPWorld (Liao et al. 2025, NAR): 998 RNA-binding Pfam domains

## LUCA PROTEOME RECONSTRUCTIONS

| Study | Year | LUCA gene count | Method | ID system |
|-------|------|-----------------|--------|-----------|
| Koonin 2003 | 2003 | 500-600 | Comparative genomics | COG |
| Delaye et al. 2005 | 2005 | 115 Pfam motifs | Phylogenomics | Pfam |
| Ranea et al. 2006 | 2006 | 140 domains | SCOP structural | SCOP |
| Weiss et al. 2016 | 2016 | 355 families | Phylogenetic (strict) | COG |
| Crapitto et al. 2022 | 2022 | 366 clusters (consensus of 8 studies) | Meta-analysis | eggNOG |
| Wehbi et al. 2024 | 2024 | 969 Pfam / 445 clans | Gene-tree reconciliation | **Pfam** |
| Moody et al. 2024 | 2024 | 399 (strict) to 2,855 (weighted) | Tree reconciliation | Gene families |

**Best dataset for our purposes: Wehbi 2024** — uses Pfam IDs, downloadable, largest count.

## ANCIENT NUCLEIC ACID BINDING

### What fraction of ancient folds are NA-binding?
- Alva et al. 2015: **13/40 (33%)** of ancient peptide fragments bind nucleic acids
- Of the **9 most ancient folds, 6 (67%)** contain NA-binding fragments
- 8 of the 13 NA-binding fragments are in ribosomal proteins

### The 9 most ancient folds (Caetano-Anollés):
1. P-loop NTPase — nucleotide hydrolysis
2. DNA/RNA-binding 3-helical bundle — NA binding
3. TIM barrel — metabolic enzymes
4. Rossmann fold — oxidoreductase
5. Ferredoxin-like — iron-sulfur electron transfer
6. Flavodoxin-like — electron transfer
7. RNase H-like — nucleic acid processing
8. OB-fold — NA binding (ssDNA, RNA, tRNA)
9. SAM-dependent methyltransferases — methylation

**Split: ~33% NA-related, ~44% metabolic, ~22% bridging/dual**

## EVIDENCE FOR ANCIENT GENERALISM

| Study | Year | Finding |
|-------|------|---------|
| Tran et al. (Nature Biotech) | 2024 | Ancestral Cas12a (~3 Gyr) bound RNA, DNA, ssDNA equally |
| Yagi & Tagami (Nature Comms) | 2024 | Ancestral beta-barrels retained both RNA and DNA affinity |
| Weil-Ktorza et al. (Angew Chem) | 2025 | HhH motif is "ambidextrous" — binds DNA of either chirality |
| Aravind et al. (Genome Biol) | 2003 | Alba: RNA-binding ancestral, DNA-binding derived (crenarchaea) |
| Jiang et al. (JBC) | 1997 | CSD: RNA chaperone (bacteria) → DNA TF (eukaryotes) |

**Emerging consensus: ancient proteins were generalist nucleic acid binders.**

## KNOWN DUAL RNA/DNA-BINDING DOMAIN FAMILIES

From Hudson & Ortlund 2014 (Nature Rev Mol Cell Biol):

| Domain | RNA function | DNA function |
|--------|-------------|--------------|
| OB-fold | tRNA binding, rRNA binding, mRNA binding | ssDNA binding (RPA, SSB, POT1) |
| CSD | RNA chaperone (CspA) | Y-box TF (YB-1) |
| KH | mRNA binding (hnRNP K, FMRP) | ssDNA binding (PCBP) |
| Alba | RNase P/MRP subunits (eukarya) | Chromatin protein (crenarchaea) |
| HTH | RNA binding (Aca2, 6S RNA) | DNA TF (primary function) |
| S1 | mRNA binding on ribosome | Minor ssDNA binding |
| Zinc finger (C2H2) | 5S rRNA binding (TFIIIA F4-6) | 5S rDNA binding (TFIIIA F1-3) |
| RRM | Pre-mRNA splicing, mRNA stability | ssDNA in some contexts |
| Sm/Lsm | snRNA binding | RNA-exclusive (control) |

## WHAT AF3/BOLTZ-2 ENABLES (NEW)

- Can predict protein-RNA and protein-DNA complexes (AF2 could not)
- PNAbind discriminates RNA vs DNA binding with AUROC 0.92
- NA-MPNN (Baker lab, 2025) unifies protein-DNA and protein-RNA analysis
- Boltz-2 (MIT license) matches AF3 accuracy, freely installable
- **Nobody has used these tools to compare RNA vs DNA binding across species**

## KEY LIMITATION

AF3 was trained on modern structures. Truly ancient structural states (extinct
folds) may be invisible to it. The DZBB fold (Yagi & Tagami 2024) was not
predicted by any ML method. Results for highly divergent ancient sequences should
be interpreted with caution.

## FUNCTIONAL COMPOSITION OF LUCA (partial knowledge)

From Crapitto et al. 2022 (consensus of 8 studies):
- Top enriched GO terms: mRNA binding, rRNA binding, aminoacyl-tRNA ligase activity
- Translation-related terms dominate
- 169 nonredundant EC codes in LUCA
- **But nobody computed a clean percentage breakdown**

From Anantharaman et al. 2002:
- ~100 protein domains involved in RNA metabolism across all life
- ~40-45 of those trace to LUCA
- RNA metabolism = 3-11% of modern proteomes (highest in small-genome bacteria)

## THE UNRESOLVED DEBATE

**Metabolism-first** (Caetano-Anollés): oldest folds are TIM barrel, Rossmann → metabolic.
Translation machinery evolved later.

**RNA-cofactor-first** (Lupas/Alva): oldest peptide fragments are ribosomal RNA-binding.
Proteins began as RNA cofactors.

Our project could contribute: within the LUCA proteome, are RNA-binding domains
systematically older or younger than metabolic domains? (Testable with Wehbi age strata.)
