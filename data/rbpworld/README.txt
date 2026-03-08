RNA-Binding Domain (RBD) Pfam accessions
=========================================

Sources:
  RBPWorld (Liao et al. 2025, NAR) - Families API: 849 domains
  EuRBPDB (Liao et al. 2020, NAR) - 791 HMM profiles: 791 domains
  RBPWorld homepage claims 998 total HMM profiles

Combined unique domain names: 929
Combined unique Pfam IDs: 916

Source breakdown:
  In both: 711
  RBPWorld only: 138
  EuRBPDB only: 80

Files:
  rbd_pfam_mapping.tsv - domain name -> Pfam ID with source
  rbd_pfam_ids.txt - unique Pfam IDs only (one per line)

Note: The RBPWorld homepage says '998 Families' but the Families
browsing page (http://research.gzsys.org.cn/rbpworld/#/family)
only lists 849. The 998 figure refers to the total HMM profiles
used for searching, but only 849 found proteins in any of the
445 eukaryotic species in RBPWorld. The other ~149 profiles are
likely prokaryotic or viral RNA-binding domains with no eukaryotic hits.

Data obtained: 2026-02-21
Method: RBPWorld API (http://research.gzsys.org.cn/rbpworld/API/)
        + EuRBPDB HMM download (http://eurbpdb.gzsys.org.cn/)
