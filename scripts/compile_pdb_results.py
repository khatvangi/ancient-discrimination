#!/usr/bin/env python3
"""
compile_pdb_results.py
merge PDB NA-complex query results with InterPro mapping and cross-domain info.
produce final TSV table for family selection.
"""

# PDB query results (from query_pdb_na_complexes.py)
pdb_results = {
    "PF00013": {"name": "KH_1", "rna": 23, "dna": 10, "examples": "1EC6,2ANN,2ANR,2N8L,2N8M"},
    "PF00575": {"name": "S1", "rna": 308, "dna": 257, "examples": "1Y1W,1Y1Y,1Y77,2B63,2BX2"},
    "PF01479": {"name": "S4", "rna": 1689, "dna": 68, "examples": "1EG0,1FJG,1FKA,1HNW,1HNX"},
    "PF01472": {"name": "PUA", "rna": 29, "dna": 0, "examples": "1R3E,1ZE2,2AB4,2HVY,2RFK"},
    "PF00633": {"name": "HHH", "rna": 0, "dna": 8, "examples": "2BGW,2OWO,4GLX,6SXB,7EF8"},
    "PF02171": {"name": "Piwi", "rna": 75, "dna": 11, "examples": "1YTU,2BGG,2F8S,2F8T,4F1N"},
    "PF00009": {"name": "GTP_EFTU", "rna": 390, "dna": 14, "examples": "1LS2,1MJ1,1OB2,1ZC8,1ZN0"},
    "PF00133": {"name": "tRNA-synt_1", "rna": 16, "dna": 0, "examples": "3ZGZ,3ZJT,3ZJU,3ZJV,4AQ7"},
    "PF00152": {"name": "tRNA-synt_2", "rna": 12, "dna": 0, "examples": "1ASY,1ASZ,1C0A,1IL2,3KFU"},
    "PF00347": {"name": "Ribosomal_L6", "rna": 1663, "dna": 53, "examples": "1EG0,1FFK,1JJ2,1K73,1K8A"},
    "PF00573": {"name": "Ribosomal_L4", "rna": 1868, "dna": 60, "examples": "1FFK,1J5A,1JJ2,1JZX,1JZY"},
    "PF00466": {"name": "Ribosomal_L10", "rna": 538, "dna": 34, "examples": "2QA4,3J16,3J77,3J78,3J7P"},
    "PF00565": {"name": "Ribosomal_S2", "rna": 0, "dna": 0, "examples": "none"},
    "PF00749": {"name": "tRNA-synt_1c", "rna": 27, "dna": 0, "examples": "1EUQ,1EUY,1EXD,1G59,1GSG"},
    "PF00750": {"name": "tRNA-synt_1d", "rna": 4, "dna": 0, "examples": "1F7U,1F7V,5B63,5YYN"},
    "PF00270": {"name": "DEAD", "rna": 284, "dna": 56, "examples": "2DB3,2HYI,2J0Q,2J0S,2XB2"},
    "PF00448": {"name": "SRP54", "rna": 17, "dna": 0, "examples": "1QZW,2J28,2V3C,2XXA,3NDB"},
    "PF00580": {"name": "UvrD-helicase", "rna": 0, "dna": 32, "examples": "1UAA,1W36,2IS1,2IS2,2IS4"},
    "PF00587": {"name": "tRNA-synt_2b", "rna": 18, "dna": 0, "examples": "1H4Q,1H4S,1KOG,1QF6,3W3S"},
    "PF00488": {"name": "MutS_V", "rna": 1, "dna": 43, "examples": "7QV3,1E3M,1NG9,1OH5,1OH6"},
    "PF00136": {"name": "DNA_pol_B", "rna": 9, "dna": 196, "examples": "4FXD,4Q5V,7N2M,8FOH,8FOJ"},
    "PF09339": {"name": "HTH_IclR", "rna": 0, "dna": 0, "examples": "none"},
    "PF01555": {"name": "N6_N4_Mtase", "rna": 0, "dna": 1, "examples": "6PBD"},
    "PF00398": {"name": "RrnaAD", "rna": 27, "dna": 15, "examples": "3FTF,4ADV,6YMW,7V2M,7V2P"},
    "PF00579": {"name": "tRNA-synt_1b", "rna": 11, "dna": 0, "examples": "1J1U,2AKE,2AZX,2DLC,2DR2"},
    "PF00588": {"name": "SpoU_methylase", "rna": 7, "dna": 0, "examples": "6YXY,7AM2,7OI6,9HCC,9HCD"},
    "PF00589": {"name": "Phage_integrase", "rna": 0, "dna": 33, "examples": "1CRX,1DRG,1F44,1KBU,1MA7"},
    "PF01131": {"name": "Topoisom_bac", "rna": 8, "dna": 27, "examples": "9CA4,9CAG,9GDA,9GDB,9GDC"},
    "PF01336": {"name": "tRNA_bind", "rna": 14, "dna": 11, "examples": "1ASY,1ASZ,1C0A,1IL2,3AMT"},
    "PF01509": {"name": "TruB_N", "rna": 30, "dna": 0, "examples": "1K8W,1R3E,1ZE2,1ZL3,2AB4"},
    "PF01926": {"name": "MMR_HSR1", "rna": 86, "dna": 0, "examples": "1X1L,3IEV,3J8G,3R9W,3R9X"},
    "PF02272": {"name": "DHHA1", "rna": 5, "dna": 1, "examples": "3WQY,3WQZ,5JJU,5O4Z,5O58"},
    "PF01042": {"name": "Ribonuclease_T2", "rna": 0, "dna": 0, "examples": "none"},
    "PF03372": {"name": "Exo_endo_phos", "rna": 4, "dna": 89, "examples": "8BVH,8UW3,9HDR,9K6I,1DE8"},
    "PF01548": {"name": "Transglut_core", "rna": 0, "dna": 0, "examples": "none"},
}

# cross-domain info from InterPro mapping
# all LUCA-age domains are by definition in bacteria + archaea
# for those in audit, we have actual counts
audit_data = {
    "PF00013": {"ipr": "IPR004088", "bac": 61156, "arc": 5752, "euk": 90934},
    "PF00575": {"ipr": "IPR003029", "bac": 261478, "arc": 11949, "euk": 45160},
    "PF01479": {"ipr": "IPR002942", "bac": 225599, "arc": 5412, "euk": 45079},
    "PF01472": {"ipr": "IPR002478", "bac": 37636, "arc": 10608, "euk": 17305},
    "PF02171": {"ipr": "IPR003165", "bac": 1395, "arc": 241, "euk": 34147},
    "PF00270": {"ipr": "IPR011545", "bac": 396324, "arc": 20246, "euk": 281909},
    "PF00136": {"ipr": "IPR006134", "bac": 6669, "arc": 5728, "euk": 16654},
    "PF01555": {"ipr": "IPR002941", "bac": 68276, "arc": 5742, "euk": 450},
    "PF00588": {"ipr": "IPR001537", "bac": 150024, "arc": 2712, "euk": 11599},
    "PF00589": {"ipr": "IPR002104", "bac": 358823, "arc": 10219, "euk": 4774},
    "PF01131": {"ipr": "IPR013497", "bac": 69856, "arc": 4433, "euk": 9488},
    "PF01336": {"ipr": "IPR004365", "bac": 159228, "arc": 11283, "euk": 29013},
    "PF02272": {"ipr": "IPR003156", "bac": 107631, "arc": 9829, "euk": 5061},
}

# functional category annotation
categories = {
    "PF00013": "structural_RNA_binding",
    "PF00575": "structural_RNA_binding",
    "PF01479": "structural_RNA_binding",
    "PF01472": "structural_RNA_binding",
    "PF00633": "dual_binder_candidate",
    "PF02171": "dual_binder_candidate",
    "PF00009": "translation",
    "PF00133": "translation",
    "PF00152": "translation",
    "PF00347": "ribosomal",
    "PF00573": "ribosomal",
    "PF00466": "ribosomal",
    "PF00565": "ribosomal",
    "PF00749": "translation",
    "PF00750": "translation",
    "PF00270": "RNA_processing",
    "PF00448": "RNA_targeting",
    "PF00580": "DNA_repair",
    "PF00587": "translation",
    "PF00488": "DNA_repair_control",
    "PF00136": "DNA_replication_control",
    "PF09339": "DNA_binding_control",
    "PF01555": "DNA_modification_control",
    "PF00398": "rRNA_modification",
    "PF00579": "translation",
    "PF00588": "RNA_modification",
    "PF00589": "DNA_recombination",
    "PF01131": "DNA_topology",
    "PF01336": "translation",
    "PF01509": "RNA_modification",
    "PF01926": "ribosome_GTPase",
    "PF02272": "nuclease",
    "PF01042": "nuclease",
    "PF03372": "nuclease",
    "PF01548": "negative_control",
}

# thesis relevance scoring
# higher = more important for the thesis
def computeRelevance(pfam_id, rna, dna):
    """score how relevant this domain is for structural prediction benchmarking"""
    score = 0

    # dual-binders are most valuable
    if rna > 0 and dna > 0:
        score += 10

    # more structures = better benchmarking
    total = rna + dna
    if total >= 100:
        score += 3
    elif total >= 20:
        score += 2
    elif total >= 5:
        score += 1

    # category bonuses
    cat = categories.get(pfam_id, "")
    if "dual" in cat:
        score += 5
    elif "translation" in cat or "ribosomal" in cat:
        score += 3
    elif "RNA" in cat:
        score += 2
    elif "DNA" in cat and "control" in cat:
        score += 4  # negative controls are important

    # cross-domain availability bonus
    if pfam_id in audit_data:
        d = audit_data[pfam_id]
        if d["bac"] > 1000 and d["arc"] > 100:
            score += 2

    return score


def main():
    # header
    cols = [
        "pfam_id", "pfam_name", "category",
        "has_rna", "has_dna", "is_dual",
        "rna_count", "dna_count",
        "cross_domain", "relevance_score",
        "example_pdbs"
    ]
    print("\t".join(cols))

    rows = []
    for pfam_id, info in pdb_results.items():
        rna = info["rna"]
        dna = info["dna"]
        has_rna = "yes" if rna > 0 else "no"
        has_dna = "yes" if dna > 0 else "no"
        is_dual = "DUAL" if rna > 0 and dna > 0 else ""
        cat = categories.get(pfam_id, "unknown")

        # cross-domain: LUCA-age means by definition Bac+Arc
        # but some (like Piwi) have very low archaeal counts
        if pfam_id in audit_data:
            d = audit_data[pfam_id]
            if d["bac"] > 1000 and d["arc"] > 100:
                xd = "Bac+Arc+Euk"
            elif d["bac"] > 100 and d["arc"] > 10:
                xd = "Bac+Arc"
            else:
                xd = "limited"
        else:
            xd = "LUCA(assumed)"

        score = computeRelevance(pfam_id, rna, dna)

        rows.append((score, pfam_id, info["name"], cat,
                     has_rna, has_dna, is_dual,
                     rna, dna, xd, score, info["examples"]))

    # sort by relevance score descending
    rows.sort(key=lambda x: x[0], reverse=True)

    for row in rows:
        vals = row[1:]  # skip the sort key
        print("\t".join(str(v) for v in vals))


if __name__ == "__main__":
    main()
