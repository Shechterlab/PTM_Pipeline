# PTM Pipeline

Residue-centric integration and downstream analysis for the human methylarginine union, with canonical remapping, provenance retention, disorder/domain context, motif enrichment, protein prioritization, and supplemental cross-PTM comparisons.

## Current design
- Default analysis universe: the full remapped methylarginine union
- Source support: retained as metadata/QC (`exact` and `fuzzy` multi-source support), not as a headline confidence-tier framework
- Maron defaults: Larsen included by default; exclusion is a sensitivity analysis via `--exclude-pattern`
- Functional annotation: multi-label, RNA-centric protein annotation
- Protein prioritization: shrinkage-adjusted per-protein methyl-Arg enrichment, with raw Arg odds ratios kept as supplemental
- Statistics: Fisher exact p-values from SciPy, plus BH q-values for test families that are interpreted together
- Domain and disorder are separate annotation layers

## Clean layout
The code repo should stay code-only. Staged data and outputs live outside the tracked source tree and are resolved through [`hpc/ptm_env.sh`](/mnt/m/Codex/PTM_pipeline/hpc/ptm_env.sh).

Local layout in this workspace:
- code: `/mnt/m/Codex/PTM_pipeline`
- staged data: `/mnt/m/Codex/PTM_pipeline/PTM_data`
- outputs: `/mnt/m/Codex/PTM_pipeline/PTM_results`

Einstein target layout:
- project root: `/gs/gsfs0/home/dshecht1/projects/2026-04_ArgReview`
- code repo: `PTM_Pipeline_code`
- staged data: `PTM_data`
- outputs: `PTM_results`

## Required staged inputs
- `PTM_data/base/human_ptm_master.tsv`
- `PTM_data/external/Table_S5_Compiled_Methylarginine_Data.xlsx`
- `PTM_data/external/mmc2.xlsx`
- `PTM_data/context/uniprot_human_reviewed_canonical.fasta`
- `PTM_data/context/uniprot_human_reviewed_canonical_metadata.tsv`

Optional staged inputs:
- `PTM_data/context/protein2ipr.dat.gz`
- `PTM_data/context/interpro_human_reviewed_domain_like_intervals.tsv`
- `PTM_data/context/mobidb_human_reviewed_disorder_intervals.tsv`
- `PTM_data/external/mobidb_human_*.json`
- `PTM_data/condensates/cdcode_staged_membership.tsv`
- `PTM_data/condensates/protein2cdcode_v2.2.tsv`
- `PTM_data/condensates/condensates_*.csv`
- `PTM_data/condensates/proteins_*.csv`

Independent downloads to stage before a full rerun:
- Maron supplemental workbook: `Table_S5_Compiled_Methylarginine_Data.xlsx`
- ProMetheusDB methylation workbook: `mmc2.xlsx`
- Reviewed human canonical UniProt FASTA plus metadata export
- Either InterPro bulk/API-derived domain intervals or `protein2ipr.dat.gz` for staging
- Either MobiDB JSON or a precomputed reviewed-human disorder-interval table
- Optional CD-CODE downloads for condensate analysis: `protein2cdcode_v2.2.tsv`, plus `condensates_*.csv` and `proteins_*.csv` if you want species/name filtering

Condensate enrichment can run from either a staged `cdcode_staged_membership.tsv` or direct CD-CODE downloads (`protein2cdcode`, plus optional `condensates` and `proteins` tables for labels/species filtering). There is no in-repo fallback CD-CODE source bundled here.

## Main local run
All paths below are derived from [`hpc/ptm_env.sh`](/mnt/m/Codex/PTM_pipeline/hpc/ptm_env.sh).

```bash
. hpc/ptm_env.sh
export HOME=/tmp/ptm_home
export XDG_CACHE_HOME=/tmp/ptm_cache
export MPLCONFIGDIR=/tmp/mpl_ptm
mkdir -p "$HOME" "$XDG_CACHE_HOME" "$MPLCONFIGDIR"

python3 scripts/01_parse_maron_s5.py \
  --input "$PTM_MARON_XLSX" \
  --outdir "$PTM_RESULTS_ROOT/parsed_maron"

python3 scripts/01_parse_prometheus_mmc2.py \
  --input "$PTM_PROMETHEUS_MMC2" \
  --outdir "$PTM_RESULTS_ROOT/parsed_prometheus"

python3 scripts/02_remap_arg_methyl_to_canonical.py \
  --base-master "$PTM_MASTER_TSV" \
  --maron "$PTM_RESULTS_ROOT/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$PTM_RESULTS_ROOT/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --uniprot-metadata "$PTM_UNIPROT_METADATA" \
  --outdir "$PTM_RESULTS_ROOT/remapped"

python3 scripts/02_integrate_arg_methyl_sources.py \
  --premerged-remapped "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$PTM_RESULTS_ROOT/integrated"

python3 scripts/03_functional_class_enrichment.py \
  --base-master "$PTM_MASTER_TSV" \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --uniprot-metadata "$PTM_UNIPROT_METADATA" \
  --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --ontology "$PTM_CODE_ROOT/config/functional_ontology.json" \
  --domain-ontology "$PTM_CODE_ROOT/config/domain_class_ontology.json" \
  --outdir "$PTM_RESULTS_ROOT/functional_class_union"

python3 scripts/04_example_figures.py \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --outdir "$PTM_RESULTS_ROOT/example_figures"

python3 scripts/08_motif_and_arg_odds.py \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/motif_arg_odds"

python3 scripts/06_annotate_disorder_context.py \
  --sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/disorder_context"

python3 scripts/05_annotate_domain_context.py \
  --sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --position-col corrected_position \
  --outdir "$PTM_RESULTS_ROOT/domain_context"

python3 scripts/07_summarize_domain_enrichment.py \
  --annotated-sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --ontology "$PTM_CODE_ROOT/config/domain_class_ontology.json" \
  --outdir "$PTM_RESULTS_ROOT/domain_context"

python3 scripts/09_compare_ptm_contexts.py \
  --base-master "$PTM_MASTER_TSV" \
  --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
  --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --outdir "$PTM_RESULTS_ROOT/ptm_compare"
```

Supplemental clustering:
```bash
python3 scripts/10_methyl_arg_clustering.py \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
  --null-match-disorder \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering"

python3 scripts/11_compare_ptm_clustering.py \
  --base-master "$PTM_MASTER_TSV" \
  --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
  --null-match-disorder \
  --outdir "$PTM_RESULTS_ROOT/ptm_clustering"
```

For quick local exploration, override the defaults explicitly, for example `--permutations 500`.

Residual GO annotation with `clusterProfiler`:
```bash
HOME=/tmp/ptm_r_home XDG_CACHE_HOME=/tmp/ptm_r_cache TMPDIR=/tmp \
Rscript scripts/13_clusterprofiler_residual_go.R \
  --labels "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" \
  --outdir "$PTM_RESULTS_ROOT/functional_class_union/clusterprofiler_go" \
  --include_all_methyl true
```

Disease enrichment with DisGeNET via `DOSE`:
```bash
HOME=/tmp/ptm_r_home XDG_CACHE_HOME=/tmp/ptm_r_cache TMPDIR=/tmp \
Rscript scripts/14_disease_enrichment.R \
  --labels "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" \
  --outdir "$PTM_RESULTS_ROOT/functional_class_union/disease_enrichment"
```

Condensate enrichment from direct CD-CODE downloads:
```bash
python3 scripts/12_condensate_enrichment.py \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --staged-membership "$PTM_DATA_ROOT/condensates/protein2cdcode_v2.2.tsv" \
  --proteins-table "$PTM_DATA_ROOT/condensates/proteins_202603181653.csv" \
  --condensates-table "$PTM_DATA_ROOT/condensates/condensates_202603181647.csv" \
  --outdir "$PTM_RESULTS_ROOT/condensates"
```

Assembly-membership clustering from enriched CD-CODE assemblies:
```bash
python3 scripts/15_assembly_membership_clustering.py \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --condensate-enrichment "$PTM_RESULTS_ROOT/condensates/condensate_enrichment.tsv" \
  --condensate-membership "$PTM_RESULTS_ROOT/condensates/condensate_membership.tsv" \
  --functional-labels "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" \
  --outdir "$PTM_RESULTS_ROOT/assembly_clustering"
```

## HPC
Create or clone a dedicated environment, then submit the main runner:

```bash
bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml ptm_pipeline
sbatch hpc/slurm_run_integration.sh "$PWD"
```

The Slurm runner:
- sources [`hpc/ptm_env.sh`](/mnt/m/Codex/PTM_pipeline/hpc/ptm_env.sh)
- avoids `set -u`
- guards empty critical variables and empty staged files
- checks the repo for merge-conflict markers before running
- defaults clustering permutations to `100000` for final HPC runs

Useful overrides:
```bash
PTM_INTERPRO_MODE=api sbatch hpc/slurm_run_integration.sh "$PWD"
PTM_CLUSTER_RESTRICT_CONTEXT=IDR PTM_CROSS_PTM_CLUSTER_RESTRICT_CONTEXT=IDR sbatch hpc/slurm_run_integration.sh "$PWD"
PTM_MARON_EXCLUDE_PATTERN=Larsen sbatch hpc/slurm_run_integration.sh "$PWD"
```

Separate condensate job:
```bash
sbatch hpc/slurm_run_condensate.sh "$PWD"
```

## Output highlights
- `PTM_results/integrated/`: remapped union, support/QC summaries, accession audit
- `PTM_results/functional_class_union/`: multi-label broad/subclass enrichment with expanded depleted/background classes and explicit `Other / unclassified` membership tables
- `PTM_results/functional_class_union/clusterprofiler_go/`: residual GO enrichment for `Other / unclassified` and `Poorly characterized / specialized`, using the reviewed-human functional universe as background
- `PTM_results/functional_class_union/disease_enrichment/`: DisGeNET disease enrichment for the full methylome and residual ontology buckets, including a neurodegenerative-focus view
- `PTM_results/motif_arg_odds/`: motif-family enrichment including explicit `DR`/`RD` tests, reordered positional heatmaps, acidic-context asymmetry, and shrinkage-based protein prioritization
- `PTM_results/disorder_context/`: site-level disorder annotation, binary IDR enrichment, and non-overlapping within-IDR edge-distance bins (`<=5`, `6-10`, `11-20`, `21-40`, `>40`)
- `PTM_results/domain_context/`: site-level domain placement, nearest-domain-class enrichment, and non-overlapping domain-edge bins (`<=5`, `6-10`, `11-20`, `21-40`, `>40`), using strict InterPro `domain/repeat` intervals by default for architectural context
- `PTM_results/ptm_compare/`: cross-PTM IDR `log2(OR)` comparison, matched IDR site-fraction figure, and normalized PTM input QC; detailed disorder/domain context line plots are opt-in
- `PTM_results/assembly_clustering/`: clustering of top methylarginine proteins by membership in enriched CD-CODE assemblies, plus shared-membership edge and community tables
- `PTM_results/methyl_arg_clustering/`, `PTM_results/ptm_clustering/`: supplemental local-density/clustering outputs, including cumulative checkpoint summaries and non-overlapping nearest-neighbor bins (`<=5`, `6-10`, `11-20`, `21-40`, `>40`)

## Notes
- `mmc2.xlsx` is the site-level ProMetheusDB methylation supplement to include. `mmc3.xlsx` is not used as a primary site source.
- `human_ptm_master.tsv` is a convenience background/integration layer, not a ground-truth database export.
- Do not use `Rme1` vs `Rme2` as a primary split for the integrated union.
- Cross-PTM nearest-neighbor analyses are supplemental. The primary structural/context story is disorder plus domain adjacency.
- Condensate odds ratios are calculated against the human non-synthetic CD-CODE condensate universe, not against all reviewed-human proteins.
- Figure scripts emit both PDF and PNG.
