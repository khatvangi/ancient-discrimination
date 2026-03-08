# Scope and Limitations

This document explicitly states the boundaries and biases of the analysis,
as required for transparent computational biology.

---

## 1. Census limitations (Phase 1) — annotation asymmetry

The census of LUCA's RNA-binding and DNA-binding proteome uses modern
functional annotation systems: Pfam domain classifications, Gene Ontology
(GO) terms, and experimentally curated databases (RBPWorld, EuRBPDB).

**RNA-annotation asymmetry**: RNA-binding annotations draw on three tiers of
evidence — GO-derived terms (487 Pfam domains), experimentally curated
CLIP/RIC/OOPS data from RBPWorld and EuRBPDB (916 domains), and their union
(968 domains). DNA-binding annotations rely solely on GO:0003677 descendants
(385 domains). This 2.4:1 ratio of annotation source coverage means the
observed 2.6:1 LUCA RNA:DNA ratio may partly reflect annotation asymmetry
rather than biological reality. The functional profiles (translation vs.
repair) are structurally different and more robust to this bias than the
raw counts.

**Modern annotation applied to ancient domains**: Pfam domains are classified
by their function in extant organisms. A domain called "RNA-binding" today
may have bound both RNA and DNA in LUCA. The census counts modern functional
labels, not ancestral biochemical activities.

**Zero observed dual-binders is an artifact**: No LUCA domain is annotated as
both RNA-binding and DNA-binding. This reflects classification conventions
(each domain gets one primary annotation), not the absence of dual function.

---

## 2. Structural sample bias (NAS-Bench)

The NAS-Bench paired structural contact analysis covers 19 convergence
table rows representing 18 independent domain families (see section 5
below for the PF07521/PF10996 deduplication). These are families that have
BOTH RNA and DNA co-crystal structures deposited in the PDB.

**PDB-availability bias**: only 24 of 87 LUCA+preLUCA NA-binding domains
(28%) have both RNA and DNA co-crystal structures. This sample is biased
toward experimentally tractable, well-studied families — typically those
with clear biological interest, good expression systems, and amenable
crystallization properties. The 63 families without paired structures are
not characterized by this analysis. It is an informative discovery set,
not a neutral sample of ancient NA-binding domains.

**Within-family PDB variation**: different PDB entries for the same Pfam
domain can yield different SI/DI values because they represent different
proteins, different organisms, and different experimental conditions.
Cross-validation between hand-curated (v2) and automated (full_luca)
PDB selections shows sign agreement in 6/9 overlapping families, with
the 3 flips occurring near the generalist boundary (|delta-SI| near 0).

**Bootstrap CI width**: Many families have wide 95% CIs on DI due to small
numbers of contacting residues. Category assignments (generalist/moderate/
specialist) are robust for extreme DI values but uncertain near thresholds
(0.10, 0.25).

---

## 3. ASR limitations — probabilistic reconstruction, not direct observation

ASR infers root-state amino acid identities from extant sequences using
maximum-likelihood phylogenetic models (IQ-TREE 2, LG+G4 substitution
model). These are probabilistic reconstructions, not direct observations
of extinct ancestors.

**Model assumptions**: the LG+G4 model assumes time-reversible amino acid
substitution with gamma-distributed rate variation across sites. This may
not capture lineage-specific compositional biases or functional constraints
on ancient NA-binding residues.

**Surviving representatives**: all conclusions concern the surviving
representatives of ancient protein lineages. Extinct domain families that
left no modern descendants are invisible to this analysis. The domains we
reconstruct are the "winners" of evolutionary history.

**Root node uncertainty**: Root posterior probabilities vary by family
(mean PP range: 0.545-0.829). Sites with PP < 0.7 are included but
flagged. Convergence conclusions are strongest for high-PP positions.

---

## 4. AF3/Boltz-2 structural predictions (Phase 3) — planned, not yet completed

The AF3/Boltz-2 structural prediction component is planned but not yet
completed. What exists in this repository:

- **Phase 2 outputs**: family selection, ortholog identification, and
  prediction manifests (296 planned jobs across 21 domain families)
- **Phase 3 scripts**: input generation and batch run scripts for AF3
  (both local Docker and AF3 Server formats)
- **Pilot feasibility data**: 45 AF3 Server competitive binding predictions
  (documented in results/competitive_af3_methodology.md and results/journal.md)
- **Pilot finding**: 36% docking rate with systematic DNA bias; 64% of
  large-protein predictions were uninformative (probe too small)

The full local AF3 batch (93 jobs for proteins <=700 aa) was partially run
but not completed. No Phase 3 results are used in the current census,
NAS-Bench, or ASR analyses — those analyses rely entirely on existing PDB
experimental structures and extant/reconstructed sequences.

---

## 5. The 18 vs 19 families clarification

The convergence_master.tsv file contains 19 rows. However, PF07521 (RMMBL,
metallo-beta-lactamase fold) and PF10996 (Beta-Casp, beta-CASP domain) are
different Pfam domains from different clans that co-occur in the same
protein and were analyzed using the same PDB structure (9BCU). They have
identical DI values (both 0.4559) because they share the same source
protein and the same RNA/DNA co-crystal contacts.

**Correct count**: 19 rows in the convergence table correspond to 18
independent domain families. All population-level statistics (e.g., "7/18
generalist", "13/18 weak discrimination") use the deduplicated count of 18.

---

## 6. Discrimination Index (DI) metric limitations

The DI metric (|SI_RNA - SI_DNA|) classifies the contact chemistry at the
protein-nucleic acid interface. It measures what fraction of contacts
involve sugar/base/phosphate atoms and how similar that profile is between
RNA and DNA co-crystal structures.

**What DI measures**: structural contact geometry — whether a domain grips
the backbone (low SI, chemically identical for RNA and DNA) or reads bases
(high SI, potentially substrate-specific).

**What DI does not measure**:
- Binding affinity (Kd) or thermodynamic preference
- Biological specificity in vivo (which depends on cellular context,
  competition, localization, and post-translational modifications)
- Kinetic discrimination (kon/koff differences)
- The effect of flanking domains, cofactors, or multimerization

A domain with DI = 0 (structural generalist) could still show strong
RNA/DNA preference in vivo through mechanisms not captured by static
crystal-structure contacts.

---

## What the data support

The evidence supports this claim: "A substantial fraction of ancient
NA-binding domain families with available paired structural evidence show
weak RNA/DNA discrimination by structural-contact metrics, consistent with
ancestral generalism."

The evidence does NOT support: "All ancient proteins were generalist NA
binders" or "No ancient protein could distinguish RNA from DNA." Five of
eighteen families show specialist-level discrimination (DI > 0.25), and
the structural sample is not exhaustive.
