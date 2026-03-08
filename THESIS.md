# The Evolution of Nucleic Acid Discrimination

## One-sentence claim

Ancient nucleic-acid-binding domains show widespread weak RNA/DNA discrimination
in structural-contact analyses, consistent with ancestral generalism — binding
specificity appears to be a derived trait that emerged through domain
diversification after LUCA.

## Scope and caveats

This claim rests on three evidence types, each with inherent limitations:

1. **Annotation-based census**: the LUCA RNA/DNA-binding fractions (54 RNA, 21 DNA
   out of 871 LUCA-age Pfam domains) use modern functional annotation (Pfam, GO,
   RBPWorld). RNA annotations draw on experimental CLIP/RIC/OOPS data while DNA
   annotations rely solely on GO-derived terms, creating an asymmetry in coverage.

2. **Structural contact analysis**: the NAS-Bench paired SI comparison covers 19
   Pfam rows (18 independent domain families) that have both RNA and DNA co-crystal
   structures in PDB. This sample is biased toward experimentally tractable,
   well-studied families — it is not a random draw from LUCA's proteome.

3. **Ancestral state reconstruction**: ASR infers root-state amino acid identities
   from extant sequences using probabilistic models (IQ-TREE, LG+G4). These are
   statistical reconstructions of the most likely ancestral states, not direct
   observations of extinct proteins. Conclusions concern the surviving
   representatives of ancient lineages.

## The philosophical thread

This connects to the deeper question in origin-of-life biology: how did molecular
specificity emerge from molecular promiscuity? Just as aminoacyl-tRNA synthetases
(aaRS) evolved to discriminate between similar amino acids, nucleic-acid-binding
domains evolved to discriminate between RNA and DNA. The evolution of discrimination
is a universal theme in molecular evolution — it is how biology creates meaning
from chemistry.

## What we claim (novel)

1. The RNA-binding fraction of LUCA's proteome has never been quantified.
   We provide that number for the first time, using the intersection of
   Wehbi et al. 2024 LUCA Pfam classifications with RBPWorld/EuRBPDB
   RNA-binding domain lists. This is an annotation-based census, not a
   direct measurement of binding.

2. Among LUCA-age domains with paired RNA and DNA co-crystal structures,
   structural contact analysis (NAS-Bench) reveals that 7/18 (39%) show
   near-identical contact profiles for RNA and DNA (DI < 0.10), and 13/18
   (72%) show weak discrimination (DI < 0.25). This is consistent with
   ancestral generalism, though the PDB-available sample is not exhaustive.

3. The timing and mechanism of the RNA/DNA discrimination transition may be
   traceable computationally using structural prediction across orthologs
   from all three domains of life. This remains a hypothesis for future
   structural prediction work (Phase 3, planned but not yet completed).

## What supports this (existing evidence)

- Tran et al. 2024 (Nature Biotech): ancestral Cas12a bound RNA, DNA, ssDNA equally
- Yagi & Tagami 2024 (Nature Comms): ancestral beta-barrels retained dual binding
- Weil-Ktorza et al. 2025 (Angew Chem): HhH motif is "ambidextrous"
- Alba domain: RNA-binding ancestral, DNA-binding derived (Aravind 2003)
- CSD: RNA chaperone in bacteria, DNA TF in eukaryotes
- Alva et al. 2015: 33% of ancient peptide fragments bind nucleic acids

## What doesn't exist yet (our gaps to fill)

- GAP-1: LUCA RNA-binding fraction (Wehbi 969 Pfam x RBPWorld 998 Pfam intersection)
- GAP-2: Systematic AF3 differential binding (same domain + RNA vs + DNA) across species
- GAP-3: Ancestral generalism quantification via PNAbind scoring
- GAP-4: Alva 40 peptides mapped onto LUCA domains
- GAP-5: Functional composition of LUCA proteome (pie chart)

## Target venue

eLife, MBE, or Nature Communications — depending on strength of results.
Current evidence (census + NAS-Bench + ASR) supports a census-focused paper.
AF3/Boltz-2 structural predictions (Phase 3) are planned but not yet completed;
they would strengthen the paper if successful but are not required for the
core census + structural-contact claims.

## Key datasets (all downloadable)

- Wehbi et al. 2024: 969 LUCA Pfam domains (GitHub/Figshare)
- RBPWorld / Liao et al. 2025: 998 RNA-binding Pfam domains
- Alva et al. 2015: 40 ancient peptide fragments (eLife supplement)
- Kolodny et al. 2021: 525 bridging themes (MBE supplement)
- AlphaFold Database: predicted structures
- Our 51 GREEN families: already collected and QC'd
