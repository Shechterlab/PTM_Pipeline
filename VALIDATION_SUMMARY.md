# Validation summary

This package was re-audited after the first repair pass.

## What was validated

### Static validation
- `python -m py_compile scripts/*.py` completed successfully after the final patch set.

### Synthetic smoke validation
- `scripts/99_smoke_test_pipeline.py` completed successfully.
- It exercises:
  - `10_methyl_arg_clustering.py`
  - `12_cross_ptm_neighbor_enrichment.py`
  - `13_condensate_ptm_enrichment.py`
  - self-exclusion behavior in `clustering_utils.nearest_distances_to_other(..., exclude_self=True)`

### Actual bundled-data validation
Using the bundled base master, Maron S5, and ProMetheus mmc2 inputs:
- `01_parse_maron_s5.py` completed successfully
- `01_parse_prometheus_mmc2.py` completed successfully
- `02_integrate_arg_methyl_sources.py` direct integration path completed successfully
- `14_union_source_audit_and_review_figures.py` completed successfully on the integrated union

## Actual bundled-data integration summary
- union source rows: 20,079
- deduplicated methyl-Arg sites: 11,800
- deduplicated proteins: 4,371
- naive canonical accessions: 4,354
- exact multi-source sites: 1,462
- fuzzy multi-source sites: 1,604
- state-annotated sites: 1,303
- source rows needing sequence remap in direct mode: 40

## Actual exact-overlap highlights
- Maron exact sites: 11,511
- ProMetheus exact sites: 1,303
- UniProt exact sites: 741
- IEDB exact sites: 85
- Maron ∩ ProMetheus exact overlap: 1,094
- Maron ∩ UniProt exact overlap: 700
- ProMetheus ∩ UniProt exact overlap: 354
- IEDB ∩ Maron exact overlap: 45

## Newly added or strengthened components
- Global and per-family BH/FDR correction support
- Homotypic self-neighbor exclusion in directed nearest-neighbor analyses
- Stratified methyl-Arg neighbor analysis (disorder and local low-complexity strata)
- Robust condensate enrichment parsing, role aggregation, and q-value reporting
- Optional dbPTM ingestion path for expanding the methylarginine union
- Review-grade union audit figure module
