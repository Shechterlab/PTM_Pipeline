# PTM pipeline v4

This package integrates the existing public human PTM master table with methylarginine-specific external resources.

Clean mental model:
- Raw sources: Maron Table S5, ProMetheusDB `mmc2.xlsx`, and the repo-local `human_ptm_master.tsv`
- Intermediate convenience layer: `human_ptm_master.tsv` is a synthetic public PTM background table, not a single ground-truth export
- Final integration layer: accession-level methylarginine union with per-site provenance, canonical-accession audit fields, and state annotations when present

Included source parsers:
- Maron Table S5 Excel parser with optional Larsen-sensitivity exclusion (default main run keeps Larsen to preserve the total union)
- ProMetheusDB supplementary `mmc2.xlsx` parser (`RK-methylsites` sheet)
- Optional dbPTM experimental methylation parser for extending the methylarginine union when a dbPTM export is staged

Key integration logic:
- Preserves exact site support across source families, with corrected source-family labeling for UniProt / IEDB / iPTMnet provenance
- Adds fuzzy support clustering within the same UniProt accession and residue, default tolerance ±2 aa
- Preserves original accession alongside a naive canonicalized accession audit field
- Includes a source-union audit module that summarizes exact overlaps, source combinations, and review-ready descriptive figures for the full methylarginine union
- Preserves source-provided methyl-state labels as provenance fields, but these are not intended to drive primary downstream splits
- Supports staged canonical FASTA and InterPro interval downloads for post-remap domain-context annotation
- Supports staged disorder interval inputs for IDR/boundary enrichment and cross-PTM context comparison

## Expected inputs
Place or symlink these files:
- `data/base/human_ptm_master.tsv` (existing iPTMnet/UniProt-based master)
- `data/external/Table_S5_Compiled_Methylarginine_Data.xlsx`
- `data/external/mmc2.xlsx`
- Optional: a dbPTM methylation export compatible with `scripts/01_parse_dbptm_methylation.py`
- For condensate analysis, place raw CD-CODE downloads in `data/external/cd_code/`:
  - `proteins_202603181653.csv`
  - `protein2cdcode_v2.2.tsv`
  - `condensates_202603181647.csv`

## Local run order
```bash
python scripts/01_parse_maron_s5.py --input data/external/Table_S5_Compiled_Methylarginine_Data.xlsx --outdir results/parsed_maron
python scripts/01_parse_prometheus_mmc2.py --input data/external/mmc2.xlsx --outdir results/parsed_prometheus
# Optional if you stage a dbPTM experimental methylation export
python scripts/01_parse_dbptm_methylation.py --input data/external/dbptm_methylation.tsv.gz --outdir results/parsed_dbptm
# Stage canonical FASTA and UniProt metadata for remapping
python scripts/00_stage_context_sources.py --outdir data/context --skip-interpro
# Remap all source rows onto canonical reviewed coordinates
python scripts/02_remap_arg_methyl_to_canonical.py --base-master data/base/human_ptm_master.tsv --maron results/parsed_maron/maron_s5_arg_methyl_sites.tsv --prometheus results/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv --dbptm results/parsed_dbptm/dbptm_arg_methyl_sites.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --uniprot-metadata data/context/uniprot_human_reviewed_canonical_metadata.tsv --outdir results/remapped
# Integrate on remapped canonical positions
python scripts/02_integrate_arg_methyl_sources.py --premerged-remapped results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv --outdir results/integrated
# Or integrate directly without remapping, optionally adding dbPTM if staged
python scripts/02_integrate_arg_methyl_sources.py --base-master data/base/human_ptm_master.tsv --maron results/parsed_maron/maron_s5_arg_methyl_sites.tsv --prometheus results/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv --dbptm results/parsed_dbptm/dbptm_arg_methyl_sites.tsv --outdir results/integrated_direct
python scripts/03_functional_class_enrichment.py --base-master data/base/human_ptm_master.tsv --integrated-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --ontology config/functional_ontology.json --interpro-intervals data/context/interpro_human_reviewed_domain_like_intervals.tsv --outdir results/functional_class_union  # multi-label RNA-centric categories, optionally supplemented by InterPro domain text
python scripts/04_example_figures.py --integrated-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --outdir results/example_figures
python scripts/08_motif_and_arg_odds.py --integrated-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --disorder-intervals data/context/mobidb_human_reviewed_disorder_intervals.tsv --outdir results/motif_arg_odds  # shrinkage-adjusted protein ranking, matched-control motif families, de novo centered k-mers, and motif feature model
# Sliding-distance clustering analysis with a within-protein permutation null
python scripts/10_methyl_arg_clustering.py --integrated-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --outdir results/methyl_arg_clustering
# Cross-PTM clustering comparison using the same sliding-distance / permutation framework
python scripts/11_compare_ptm_clustering.py --base-master data/base/human_ptm_master.tsv --integrated-arg-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --outdir results/ptm_clustering
# Directed methyl-Arg neighbor enrichment against homotypic and heterotypic PTM classes
python scripts/12_cross_ptm_neighbor_enrichment.py --base-master data/base/human_ptm_master.tsv --integrated-arg-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --disorder-intervals data/context/mobidb_human_reviewed_disorder_intervals.tsv --outdir results/cross_ptm_neighbors
# You can also point the condensate module directly at raw CD-CODE downloads if you do not want a separate staging step
python scripts/13_condensate_ptm_enrichment.py --base-master data/base/human_ptm_master.tsv --integrated-arg-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --cdcode-proteins data/external/cd_code/proteins_202603181653.csv --cdcode-protein-map data/external/cd_code/protein2cdcode_v2.2.tsv --cdcode-condensates data/external/cd_code/condensates_202603181647.csv --disorder-intervals data/context/mobidb_human_reviewed_disorder_intervals.tsv --outdir results/condensate_enrichment
# Stage raw CD-CODE downloads into the protein table expected by the condensate module
python scripts/00_stage_cd_code_condensates.py --cdcode-proteins data/external/cd_code/proteins_202603181653.csv --cdcode-protein-map data/external/cd_code/protein2cdcode_v2.2.tsv --cdcode-condensates data/external/cd_code/condensates_202603181647.csv --out data/context/cd_code_condensate_proteins.tsv
# Protein-level condensate enrichment from either the staged table or the raw CD-CODE trio
python scripts/13_condensate_ptm_enrichment.py --base-master data/base/human_ptm_master.tsv --integrated-arg-sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --condensate-proteins data/context/cd_code_condensate_proteins.tsv --disorder-intervals data/context/mobidb_human_reviewed_disorder_intervals.tsv --outdir results/condensate_enrichment
# Review-grade union-source audit and figure set
python scripts/14_union_source_audit_and_review_figures.py --all-rows results/integrated_direct/human_arg_methyl_union_all_rows.tsv --dedup-sites results/integrated_direct/human_arg_methyl_union_dedup_by_site.tsv --outdir results/review_audit
# Synthetic smoke test for the repaired clustering / neighbor / condensate modules
python scripts/99_smoke_test_pipeline.py
# If you stage disorder intervals (for example MobiDB-derived intervals) as canonical accession/start/end rows:
python scripts/06_annotate_disorder_context.py --sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --disorder-intervals data/context/mobidb_human_reviewed_disorder_intervals.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --outdir results/disorder_context
# Stage InterPro only for remapped canonical proteins once remapping is done
python scripts/00_stage_context_sources.py --outdir data/context --accessions results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv --interpro-source bulk --interpro-protein2ipr data/context/protein2ipr.dat.gz --bulk-parser pandas
# After canonical remapping adds corrected positions, annotate domain context separately from disorder context
python scripts/05_annotate_domain_context.py --sites results/integrated/human_arg_methyl_union_dedup_by_site.tsv --interpro-intervals data/context/interpro_human_reviewed_domain_like_intervals.tsv --position-col corrected_position --outdir results/domain_context
python scripts/07_summarize_domain_enrichment.py --annotated-sites results/domain_context/sites_with_domain_context.tsv --interpro-intervals data/context/interpro_human_reviewed_domain_like_intervals.tsv --canonical-fasta data/context/uniprot_human_reviewed_canonical.fasta --outdir results/domain_context  # matched non-methyl-Arg controls and domain-boundary plots
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

Then submit the main HPC workflow:
```bash
sbatch hpc/slurm_run_integration.sh "$PWD"
```

Recommended Einstein layout for this review project:
```bash
PTM_PROJECT_HOME=/gs/gsfs0/home/dshecht1/projects/2026-04_ArgReview/ptm_pipeline
PTM_DATA_ROOT=$PTM_PROJECT_HOME/PTM_data
PTM_RESULTS_ROOT=$PTM_PROJECT_HOME/PTM_results
```

The Slurm scripts now source `hpc/ptm_env.sh` and default to that exact layout. Keep the Git repo anywhere convenient, and stage large inputs outside the repo under:
```text
$PTM_DATA_ROOT/base/
$PTM_DATA_ROOT/context/
$PTM_DATA_ROOT/external/
$PTM_DATA_ROOT/external/cd_code/
```

This avoids collisions between Git updates and staged files such as InterPro, MobiDB, UniProt FASTA, and CD-CODE exports. Outputs now default to `$PTM_RESULTS_ROOT/` rather than `results/` inside the repo.

That main Slurm runner now also does the following automatically:
- writes the union-source review audit and figure set,
- stages and runs condensate enrichment automatically if either `data/context/cd_code_condensate_proteins.tsv` already exists or the raw CD-CODE files are present in `data/external/cd_code/`,
- uses 100,000 permutations by default for the core methyl-Arg clustering analysis,
- keeps cross-PTM neighbor enrichment as an optional supplemental step (`PTM_RUN_CROSS_PTM_NEIGHBORS=1`).

If you want condensate analysis as a dedicated separate job, use:
```bash
sbatch hpc/slurm_run_condensate.sh "$PWD"
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

If you want to override the default Einstein project paths, edit `hpc/ptm_env.sh` once or pass overrides at submit time:
```bash
PTM_DATA_ROOT=/some/other/data PTM_RESULTS_ROOT=/some/other/results sbatch hpc/slurm_run_integration.sh "$PWD"
```

## Notes
- `mmc2.xlsx` is the site-level ProMetheusDB supplement to include. `mmc3.xlsx` contains enrichment/cluster summaries and is not used as a primary site source.
- `human_ptm_master.tsv` is a repo-local integration/background table. Treat it as a convenience layer for coverage, not as a standalone authoritative source database.
- Current integration outputs are canonical-site based once you run the remap step. The remapper writes per-row diagnostics so unresolved source rows can be audited separately.
- For large InterPro staging runs, prefer the official download files `protein2ipr.dat.gz` and `entry.list`. If `protein2ipr.dat.gz` is staged locally, the script uses chunked pandas filtering and does not rely on file ordering.
- If you use `PTM_INTERPRO_MODE=api`, the main Slurm runner stages InterPro intervals with the per-accession InterPro API instead of the bulk file. Tune with `PTM_INTERPRO_API_WORKERS`, `PTM_INTERPRO_API_PAGE_SIZE`, `PTM_INTERPRO_API_TIMEOUT`, and `PTM_INTERPRO_API_RETRIES`.
- Disorder inputs are expected as a staged canonical-interval table with at least `canonical_UniProtAC`, `fragment_start`, and `fragment_end`. The current repo does not yet include a downloader for MobiDB, so that file must be staged separately.
- `scripts/10_methyl_arg_clustering.py` models methyl-site proximity as a sliding empirical CDF, `P(nearest methyl-Arg <= X)`, and a local-density curve, `mean other methyl-Args within X`, against a within-protein randomization null. The `3/5/10/20/50 aa` outputs are just checkpoint summaries of that full curve.
- `scripts/11_compare_ptm_clustering.py` applies the same sliding-distance framework across PTM classes and reports both raw nearest-neighbor curves and null-adjusted enrichment curves. It defaults to 100 permutations for runtime reasons; increase that on HPC if you want finer empirical p-value resolution.
- Maron Table S5 now defaults to keeping the Larsen entries so the main run represents the total union. For a sensitivity analysis, set `--exclude-pattern Larsen` or `PTM_MARON_EXCLUDE_PATTERN=Larsen`.
- Do not use `Rme1` vs `Rme2` as a primary split for the integrated mass-spec-driven union. Source-provided state labels are retained only as auxiliary provenance.
- Domain context and disorder context should remain separate layers. The domain-context script classifies sites as `in_domain`, `boundary`, `inter_domain_linker`, or `distal` relative to InterPro intervals and does not merge those labels with IDR annotations.
- Figure-producing scripts now emit PDF and PNG outputs.
- The dbPTM downloader script is separated from analysis so static files can be staged first.
- The Einstein HPC `WARNING: overwriting environment variables set in the machine` message during Conda activation is not the failure. The real failure is NumPy/Pandas/Matplotlib binary incompatibility in the shared base env, which is why the pipeline should run in a dedicated Conda environment.

- `scripts/12_cross_ptm_neighbor_enrichment.py` remains available as a supplemental analysis. It is now off by default in the main Slurm review workflow because the review-focused figures emphasize homotypic methyl-Arg clustering, motif context, domain adjacency, and condensate biology first.
- `scripts/00_stage_cd_code_condensates.py` converts raw CD-CODE downloads into `data/context/cd_code_condensate_proteins.tsv`. Put the downloaded raw files in `data/external/cd_code/` and either stage them first or pass them directly to the condensate module.
- `scripts/13_condensate_ptm_enrichment.py` tests whether PTM-bearing proteins are enriched among condensate-associated proteins from a staged CD-CODE-style export or directly from the raw CD-CODE trio. It reports both unadjusted enrichment and an adjusted logistic model including protein length, candidate-residue count, disorder fraction, and a simple low-complexity fraction.
- The central enrichment helper now treats the second argument as an inclusive universe and removes target rows from the background before building Fisher tables. This affects functional-class, disorder, domain, and cross-PTM context enrichment outputs.
- The direct `02_integrate_arg_methyl_sources.py` path now works again without the remap intermediate, although canonical remapping remains the preferred review-grade route.


### Condensate HPC controls
The main integration Slurm script and the dedicated condensate Slurm script both avoid `set -u` and respect the Einstein Conda activation quirks.

Useful environment overrides:
```bash
# Force condensate analysis on or off inside the main integration job
PTM_RUN_CONDENSATE=1 sbatch hpc/slurm_run_integration.sh "$PWD"
PTM_RUN_CONDENSATE=0 sbatch hpc/slurm_run_integration.sh "$PWD"

# Optional supplemental cross-PTM neighbor step
PTM_RUN_CROSS_PTM_NEIGHBORS=1 sbatch hpc/slurm_run_integration.sh "$PWD"

# Larsen sensitivity run instead of the default total-union run
PTM_MARON_EXCLUDE_PATTERN=Larsen sbatch hpc/slurm_run_integration.sh "$PWD"

# Override the core methyl-Arg clustering permutation count
PTM_METHYL_CLUSTER_PERMUTATIONS=100000 sbatch hpc/slurm_run_integration.sh "$PWD"

# Run experimental-only condensates or impose a confidence threshold
PTM_CONDENSATE_EXPERIMENTAL_ONLY=1 sbatch hpc/slurm_run_condensate.sh "$PWD"
PTM_CONDENSATE_MIN_CONFIDENCE=0.8 sbatch hpc/slurm_run_condensate.sh "$PWD"

# Override raw CD-CODE locations if you keep them somewhere else
PTM_CDCODE_DIR=/path/to/cd_code sbatch hpc/slurm_run_condensate.sh "$PWD"
```
