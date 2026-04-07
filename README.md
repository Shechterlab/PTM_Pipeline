# PTM pipeline v4

This package integrates the existing public human PTM master table with methylarginine-specific external resources.

Clean mental model:
- Raw sources: Maron Table S5, ProMetheusDB `mmc2.xlsx`, and the repo-local `human_ptm_master.tsv`
- Intermediate convenience layer: `human_ptm_master.tsv` is a synthetic public PTM background table, not a single ground-truth export
- Final integration layer: accession-level methylarginine union with per-site provenance, confidence tiers, canonical-accession audit fields, and state annotations when present

Included source parsers:
- Maron Table S5 Excel parser with optional Larsen exclusion
- ProMetheusDB supplementary `mmc2.xlsx` parser (`RK-methylsites` sheet)
- dbPTM static download helper and a parser preview for the methylation export

Key integration logic:
- Preserves exact site support across source families
- Adds fuzzy support clustering within the same UniProt accession and residue, default tolerance ±2 aa
- Assigns confidence tiers rather than using a hard `>=2 exact-match` filter alone
- Preserves original accession alongside a naive canonicalized accession audit field
- Preserves source-provided methyl-state labels as provenance fields, but these are not intended to drive primary downstream splits
- Supports staged canonical FASTA and InterPro interval downloads for post-remap domain-context annotation
- Supports staged disorder interval inputs for IDR/boundary enrichment and cross-PTM context comparison

## Expected inputs
Place or symlink these files:
- `data/base/human_ptm_master.tsv` (existing iPTMnet/UniProt-based master)
- `data/external/Table_S5_Compiled_Methylarginine_Data.xlsx`
- `data/external/mmc2.xlsx`

## Local run order
```bash
python scripts/01_parse_maron_s5.py --input data/external/Table_S5_Compiled_Methylarginine_Data.xlsx --outdir results/parsed_maron
python scripts/01_parse_prometheus_mmc2.py --input data/external/mmc2.xlsx --outdir results/parsed_prometheus
# Stage canonical FASTA and UniProt metadata for remapping
python scripts/00_stage_context_sources.py --outdir data/context --skip-interpro
# Remap all source rows onto canonical reviewed coordinates
python scripts/02_remap_arg_methyl_to_canonical.py --base-master data/base/human_ptm_master.tsv --maron results/parsed_maron/maron_s5_arg_methyl_sites.tsv --prometheus results/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --uniprot-metadata data/context/uniprot_human_reviewed_canonical_metadata.tsv --outdir results/remapped
# Integrate on remapped canonical positions
python scripts/02_integrate_arg_methyl_sources.py --premerged-remapped results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv --outdir results/integrated
python scripts/03_functional_class_enrichment.py --base-master data/base/human_ptm_master.tsv --integrated-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --ontology config/functional_ontology.json --outdir results/functional_class_union
python scripts/04_example_figures.py --integrated-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --outdir results/example_figures
python scripts/08_motif_and_arg_odds.py --integrated-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --outdir results/motif_arg_odds
# If you stage disorder intervals (for example MobiDB-derived intervals) as canonical accession/start/end rows:
python scripts/06_annotate_disorder_context.py --sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --disorder-intervals data/context/mobidb_human_reviewed_disorder_intervals.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --outdir results/disorder_context
# Stage InterPro only for remapped canonical proteins once remapping is done
python scripts/00_stage_context_sources.py --outdir data/context --accessions results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv --interpro-source bulk --interpro-protein2ipr data/context/protein2ipr.dat.gz --bulk-parser pandas
# After canonical remapping adds corrected positions, annotate domain context separately from disorder context
python scripts/05_annotate_domain_context.py --sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --interpro-intervals data/context/interpro_human_reviewed_domain_like_intervals.tsv --position-col corrected_position --outdir results/domain_context
python scripts/07_summarize_domain_enrichment.py --annotated-sites results/domain_context/sites_with_domain_context.tsv --interpro-intervals data/context/interpro_human_reviewed_domain_like_intervals.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --outdir results/domain_context
# Cross-PTM context comparison once disorder and/or domain intervals are available
python scripts/09_compare_ptm_contexts.py --base-master data/base/human_ptm_master.tsv --integrated-arg-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --disorder-intervals data/context/mobidb_human_reviewed_disorder_intervals.tsv --interpro-intervals data/context/interpro_human_reviewed_domain_like_intervals.tsv --outdir results/ptm_compare
```

## HPC
Create a dedicated environment instead of using the shared base env:
```bash
bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml ptm_pipeline
```

The helper defaults to `conda env create`, not `mamba`, because some shared HPC installs hang under `mamba`. If you want to force `mamba`:
```bash
PTM_ENV_SOLVER=mamba bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml ptm_pipeline
```

Stage InterPro locally if you want domain-context outputs:
```bash
curl -L -o data/context/protein2ipr.dat.gz https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/protein2ipr.dat.gz
```

Then submit:
```bash
sbatch hpc/slurm_run_integration.sh "$PWD"
```

If the bulk `protein2ipr.dat.gz` route is unreliable on your cluster, use the InterPro API path instead:
```bash
PTM_INTERPRO_MODE=api PTM_INTERPRO_API_WORKERS=8 sbatch hpc/slurm_run_integration.sh "$PWD"
```

Or stage InterPro/domain context as a separate HPC job after remapping/integration:
```bash
sbatch hpc/slurm_stage_interpro_api.sh "$PWD"
```

If your HPC env uses a different Conda env name or prefix:
```bash
PTM_ENV_NAME=my_env sbatch hpc/slurm_run_integration.sh "$PWD"
# or
PTM_ENV_PREFIX=/path/to/conda/env sbatch hpc/slurm_run_integration.sh "$PWD"
```

## Notes
- `mmc2.xlsx` is the site-level ProMetheusDB supplement to include. `mmc3.xlsx` contains enrichment/cluster summaries and is not used as a primary site source.
- `human_ptm_master.tsv` is a repo-local integration/background table. Treat it as a convenience layer for coverage, not as a standalone authoritative source database.
- Current integration outputs are canonical-site based once you run the remap step. The remapper writes per-row diagnostics so unresolved source rows can be audited separately.
- For large InterPro staging runs, prefer the official download files `protein2ipr.dat.gz` and `entry.list`. If `protein2ipr.dat.gz` is staged locally, the script uses chunked pandas filtering and does not rely on file ordering.
- If you use `PTM_INTERPRO_MODE=api`, the main Slurm runner stages InterPro intervals with the per-accession InterPro API instead of the bulk file. Tune with `PTM_INTERPRO_API_WORKERS`, `PTM_INTERPRO_API_PAGE_SIZE`, `PTM_INTERPRO_API_TIMEOUT`, and `PTM_INTERPRO_API_RETRIES`.
- Disorder inputs are expected as a staged canonical-interval table with at least `canonical_UniProtAC`, `fragment_start`, and `fragment_end`. The current repo does not yet include a downloader for MobiDB, so that file must be staged separately.
- Maron Table S5 currently defaults to excluding `Larsen`-matched references. Override `--exclude-pattern` if you want a different sensitivity run.
- Do not use `Rme1` vs `Rme2` as a primary split for the integrated mass-spec-driven union. Source-provided state labels are retained only as auxiliary provenance.
- Domain context and disorder context should remain separate layers. The domain-context script classifies sites as `in_domain`, `boundary`, `inter_domain_linker`, or `distal` relative to InterPro intervals and does not merge those labels with IDR annotations.
- Figure-producing scripts now emit PDF and PNG outputs.
- The dbPTM downloader script is separated from analysis so static files can be staged first.
- The Einstein HPC `WARNING: overwriting environment variables set in the machine` message during Conda activation is not the failure. The real failure is NumPy/Pandas/Matplotlib binary incompatibility in the shared base env, which is why the pipeline should run in a dedicated Conda environment.
