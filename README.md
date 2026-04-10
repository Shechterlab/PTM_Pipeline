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

Optional staged inputs:
- `PTM_data/context/protein2ipr.dat.gz`
- `PTM_data/context/mobidb_human_reviewed_disorder_intervals.tsv`
- `PTM_data/condensates/cdcode_staged_membership.tsv`

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
  --ontology "$PTM_CODE_ROOT/config/functional_ontology.json" \
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
  --outdir "$PTM_RESULTS_ROOT/domain_context"

python3 scripts/09_compare_ptm_contexts.py \
  --base-master "$PTM_MASTER_TSV" \
  --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
  --outdir "$PTM_RESULTS_ROOT/ptm_compare"
```

Supplemental clustering:
```bash
python3 scripts/10_methyl_arg_clustering.py \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
  --null-match-disorder \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering" \
  --permutations 500

python3 scripts/11_compare_ptm_clustering.py \
  --base-master "$PTM_MASTER_TSV" \
  --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
  --null-match-disorder \
  --outdir "$PTM_RESULTS_ROOT/ptm_clustering" \
  --permutations 200
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
- `PTM_results/functional_class_union/`: multi-label broad/subclass enrichment
- `PTM_results/motif_arg_odds/`: motif-family enrichment, positional enrichment, shrinkage-based protein prioritization
- `PTM_results/disorder_context/`: site-level disorder annotation and binary IDR enrichment
- `PTM_results/domain_context/`: site-level domain placement and domain-class enrichment
- `PTM_results/ptm_compare/`: cross-PTM IDR/context comparison
- `PTM_results/methyl_arg_clustering/`, `PTM_results/ptm_clustering/`: supplemental local-density/clustering outputs

## Notes
- `mmc2.xlsx` is the site-level ProMetheusDB methylation supplement to include. `mmc3.xlsx` is not used as a primary site source.
- `human_ptm_master.tsv` is a convenience background/integration layer, not a ground-truth database export.
- Do not use `Rme1` vs `Rme2` as a primary split for the integrated union.
- Cross-PTM nearest-neighbor analyses are supplemental. The primary structural/context story is disorder plus domain adjacency.
- Figure scripts emit both PDF and PNG.
