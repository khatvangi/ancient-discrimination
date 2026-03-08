# The Evolution of Nucleic Acid Discrimination

## Central claim

Many ancient nucleic-acid-binding domain families lack the specialist-like
atom-type contact chemistry associated with strong RNA/DNA discrimination,
consistent with ancestral generalism.

## Scope and caveats

This claim rests on three evidence types, each with inherent limitations:

1. **Annotation-based census**: the LUCA RNA/DNA-binding fractions (54 RNA, 21 DNA
   out of 871 LUCA-age Pfam domains) use modern functional annotation (Pfam, GO,
   RBPWorld). RNA annotations draw on experimental CLIP/RIC/OOPS data while DNA
   annotations rely solely on GO-derived terms, creating an asymmetry in coverage.

2. **Structural contact analysis**: the NAS-Bench paired SI comparison covers 19
   Pfam rows (18 independent domain families, after collapsing PF07521/PF10996
   which share source protein 9BCU) that have both RNA and DNA co-crystal
   structures in PDB. This is 28% of the 87 LUCA+preLUCA NA-binding domains —
   a biased sample of experimentally tractable families, not a random draw.

3. **Ancestral state reconstruction**: ASR infers root-state amino acid identities
   from extant sequences using probabilistic models (IQ-TREE, LG+G4). These are
   statistical reconstructions of the most likely ancestral states, not direct
   observations of extinct proteins. Conclusions concern the surviving
   representatives of ancient lineages.

See `SCOPE_AND_LIMITATIONS.md` for the full treatment of biases and boundaries.

## Novel contributions

1. **First census of LUCA's RNA-binding proteome.** The intersection of
   Wehbi et al. 2024 LUCA Pfam classifications with RBPWorld/EuRBPDB
   RNA-binding domain lists. This is an annotation-based census, not a
   direct measurement of binding.

2. **Systematic structural contact benchmark.** Among the 18 independent
   LUCA-age domain families with paired RNA and DNA co-crystal structures,
   7/18 (39%) show near-identical atom-type contact profiles (DI < 0.10),
   and 13/18 (72%) show weak discrimination (DI < 0.25). DI thresholds are
   heuristic operational cutoffs supported by sensitivity analysis (see
   `results/nasbench_validation.tsv`), not statistically derived boundaries.

3. **ASR convergence analysis.** Root-state residues at NA-contacting positions
   differ from modern specialist equivalents in 72% of cases, with 59% radical
   substitutions — the ancestral binding mode is largely extinct.

## Supporting literature

- Tran et al. 2024 (Nature Biotech): ancestral Cas12a bound RNA, DNA, ssDNA equally
- Yagi & Tagami 2024 (Nature Comms): ancestral beta-barrels retained dual binding
- Weil-Ktorza et al. 2025 (Angew Chem): HhH motif is "ambidextrous"
- Alba domain: RNA-binding ancestral, DNA-binding derived (Aravind 2003)
- CSD: RNA chaperone in bacteria, DNA TF in eukaryotes
- Alva et al. 2015: 33% of ancient peptide fragments bind nucleic acids

## Key datasets

- Wehbi et al. 2024: 871 LUCA-age Pfam domains (GitHub/Figshare)
- RBPWorld / Liao et al. 2025: 998 RNA-binding Pfam domains
- Alva et al. 2015: 40 ancient peptide fragments (eLife supplement)
- Kolodny et al. 2021: 525 bridging themes (MBE supplement)
