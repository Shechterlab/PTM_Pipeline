# Nuclear RNA neighbourhoods: pilot (not part of a release)

- `pilot_ri_barutcu.py`, `pilot_ri_controls.py`: are transcripts with introns whose retention changes after PRMT5
  knockdown (rMATS `RI.MATS.JCEC.txt`) associated with nuclear speckles or other nuclear bodies in APEX2-seq
  (Barutcu et al. 2022, GEO GSE176439, `GSE176439_NuclearBody_APEX_seq_gene_counts.txt.gz`)? Gene level; association is
  the bait pulldown over its own input, minus the nucleoplasmic NLS-APEX2 control, at least twofold in both replicates.
  Controls: expression and read-depth deciles (Mantel-Haenszel), and the direction of change among significant events.
- `pilot_maron_ms.py`: chromatin and chromatin-associated poly(A) RNA proteomics under GSK591 and MS023 (Maron et al.
  2022, eLife, Supplementary File 1): are methylated proteins over-represented among proteins that change?
  Uses `lib.py` beside it (the third argument is its folder).

The rMATS table is the lab's own and is not in the repository.
