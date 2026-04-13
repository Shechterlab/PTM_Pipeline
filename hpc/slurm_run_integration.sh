#!/bin/bash
#SBATCH -J ptm_integrate
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 16:00:00

set -eo pipefail

CONDA_BASE=/gs/gsfs0/hpc01/rhel8/apps/conda3
CONDA_EXE="$CONDA_BASE/bin/conda"

if [ -x "$CONDA_EXE" ]; then
  eval "$("$CONDA_EXE" shell.bash hook)"
  CONDA_HOOK_STATUS=$?
else
  . "$CONDA_BASE/bin/activate"
  CONDA_HOOK_STATUS=$?
fi

ENV_NAME=${PTM_ENV_NAME:-ptm_pipeline}
ENV_PREFIX=${PTM_ENV_PREFIX:-}
INTERPRO_MODE=${PTM_INTERPRO_MODE:-auto}
INTERPRO_API_WORKERS=${PTM_INTERPRO_API_WORKERS:-8}
INTERPRO_API_PAGE_SIZE=${PTM_INTERPRO_API_PAGE_SIZE:-200}
INTERPRO_API_TIMEOUT=${PTM_INTERPRO_API_TIMEOUT:-120}
INTERPRO_API_RETRIES=${PTM_INTERPRO_API_RETRIES:-5}
CLUSTER_PERMUTATIONS=${PTM_CLUSTER_PERMUTATIONS:-100000}
CROSS_PTM_CLUSTER_PERMUTATIONS=${PTM_CROSS_PTM_CLUSTER_PERMUTATIONS:-100000}
CLUSTER_RESTRICT_CONTEXT=${PTM_CLUSTER_RESTRICT_CONTEXT:-all}
CROSS_PTM_CLUSTER_RESTRICT_CONTEXT=${PTM_CROSS_PTM_CLUSTER_RESTRICT_CONTEXT:-all}
MARON_EXCLUDE_PATTERN=${PTM_MARON_EXCLUDE_PATTERN:-}

run_step() {
  echo "[ptm_pipeline] $1"
  shift
  "$@"
  status=$?
  if [ "$status" -ne 0 ]; then
    echo "[ptm_pipeline] step failed with exit code $status"
    echo "[ptm_pipeline] command: $*"
    exit "$status"
  fi
}

if [ "$CONDA_HOOK_STATUS" -ne 0 ]; then
  echo "[ptm_pipeline] failed to initialize Conda shell hook"
  exit "$CONDA_HOOK_STATUS"
fi

if [ -n "$ENV_PREFIX" ]; then
  conda activate "$ENV_PREFIX"
else
  conda activate "$ENV_NAME"
fi

BASE_DIR=${1:-$PWD}
. "$BASE_DIR/hpc/ptm_env.sh"

if [ -z "$PTM_DATA_ROOT" ] || [ -z "$PTM_RESULTS_ROOT" ] || [ -z "$PTM_CODE_ROOT" ]; then
  echo "[ptm_pipeline] critical PTM_* paths are empty after sourcing hpc/ptm_env.sh"
  exit 1
fi

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$PTM_RESULTS_ROOT"

run_step "repo hygiene guard" bash "$PTM_CODE_ROOT/hpc/check_repo_hygiene.sh" "$PTM_CODE_ROOT"
echo "[ptm_pipeline] code=$PTM_CODE_ROOT"
echo "[ptm_pipeline] data=$PTM_DATA_ROOT"
echo "[ptm_pipeline] results=$PTM_RESULTS_ROOT"

run_step "stage UniProt context" python "$PTM_CODE_ROOT/scripts/00_stage_context_sources.py" \
  --outdir "$PTM_CONTEXT_ROOT" \
  --skip-interpro

MARON_ARGS=(
  --input "$PTM_MARON_XLSX"
  --outdir "$PTM_RESULTS_ROOT/parsed_maron"
)
if [ -n "$MARON_EXCLUDE_PATTERN" ]; then
  MARON_ARGS+=(--exclude-pattern "$MARON_EXCLUDE_PATTERN")
fi
run_step "parse Maron S5" python "$PTM_CODE_ROOT/scripts/01_parse_maron_s5.py" "${MARON_ARGS[@]}"

run_step "parse ProMetheus mmc2" python "$PTM_CODE_ROOT/scripts/01_parse_prometheus_mmc2.py" \
  --input "$PTM_PROMETHEUS_MMC2" \
  --outdir "$PTM_RESULTS_ROOT/parsed_prometheus"

run_step "remap to canonical coordinates" python "$PTM_CODE_ROOT/scripts/02_remap_arg_methyl_to_canonical.py" \
  --base-master "$PTM_MASTER_TSV" \
  --maron "$PTM_RESULTS_ROOT/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$PTM_RESULTS_ROOT/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --uniprot-metadata "$PTM_UNIPROT_METADATA" \
  --outdir "$PTM_RESULTS_ROOT/remapped"

run_step "integrate remapped sites" python "$PTM_CODE_ROOT/scripts/02_integrate_arg_methyl_sources.py" \
  --premerged-remapped "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$PTM_RESULTS_ROOT/integrated"

run_step "functional enrichment" python "$PTM_CODE_ROOT/scripts/03_functional_class_enrichment.py" \
  --base-master "$PTM_MASTER_TSV" \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --uniprot-metadata "$PTM_UNIPROT_METADATA" \
  --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --ontology "$PTM_CODE_ROOT/config/functional_ontology.json" \
  --domain-ontology "$PTM_CODE_ROOT/config/domain_class_ontology.json" \
  --outdir "$PTM_RESULTS_ROOT/functional_class_union"

run_step "example figures" python "$PTM_CODE_ROOT/scripts/04_example_figures.py" \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --outdir "$PTM_RESULTS_ROOT/example_figures"

run_step "motif enrichment and protein prioritization" python "$PTM_CODE_ROOT/scripts/08_motif_and_arg_odds.py" \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/motif_arg_odds"

CLUSTER_ARGS=(
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
  --canonical-fasta "$PTM_CANONICAL_FASTA"
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering"
  --permutations "$CLUSTER_PERMUTATIONS"
  --restrict-context "$CLUSTER_RESTRICT_CONTEXT"
)
PTM_CLUSTER_ARGS=(
  --base-master "$PTM_MASTER_TSV"
  --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
  --canonical-fasta "$PTM_CANONICAL_FASTA"
  --outdir "$PTM_RESULTS_ROOT/ptm_clustering"
  --permutations "$CROSS_PTM_CLUSTER_PERMUTATIONS"
  --restrict-context "$CROSS_PTM_CLUSTER_RESTRICT_CONTEXT"
)
if [[ -f "$PTM_MOBIDB_INTERVALS" && -s "$PTM_MOBIDB_INTERVALS" ]]; then
  CLUSTER_ARGS+=(--disorder-intervals "$PTM_MOBIDB_INTERVALS" --null-match-disorder)
  PTM_CLUSTER_ARGS+=(--disorder-intervals "$PTM_MOBIDB_INTERVALS" --null-match-disorder)
fi
run_step "methyl-Arg clustering" python "$PTM_CODE_ROOT/scripts/10_methyl_arg_clustering.py" "${CLUSTER_ARGS[@]}"
run_step "cross-PTM clustering" python "$PTM_CODE_ROOT/scripts/11_compare_ptm_clustering.py" "${PTM_CLUSTER_ARGS[@]}"

if [[ -f "$PTM_MOBIDB_INTERVALS" && -s "$PTM_MOBIDB_INTERVALS" ]]; then
  run_step "annotate disorder context" python "$PTM_CODE_ROOT/scripts/06_annotate_disorder_context.py" \
    --sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --disorder-intervals "$PTM_MOBIDB_INTERVALS" \
    --canonical-fasta "$PTM_CANONICAL_FASTA" \
    --outdir "$PTM_RESULTS_ROOT/disorder_context"
else
  echo "[ptm_pipeline] skipping disorder-context annotation: $PTM_MOBIDB_INTERVALS missing or empty"
fi

if [[ "$INTERPRO_MODE" = "api" ]]; then
  run_step "stage InterPro domain intervals via API" python "$PTM_CODE_ROOT/scripts/00_stage_context_sources.py" \
    --outdir "$PTM_CONTEXT_ROOT" \
    --accessions "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
    --interpro-source api \
    --workers "$INTERPRO_API_WORKERS" \
    --page-size "$INTERPRO_API_PAGE_SIZE" \
    --timeout "$INTERPRO_API_TIMEOUT" \
    --retries "$INTERPRO_API_RETRIES"
elif [[ "$INTERPRO_MODE" = "bulk" || ( "$INTERPRO_MODE" = "auto" && -f "$PTM_INTERPRO_BULK" ) ]]; then
  if [[ -f "$PTM_INTERPRO_BULK" ]]; then
    run_step "stage InterPro domain intervals" python "$PTM_CODE_ROOT/scripts/00_stage_context_sources.py" \
      --outdir "$PTM_CONTEXT_ROOT" \
      --accessions "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
      --interpro-source bulk \
      --interpro-protein2ipr "$PTM_INTERPRO_BULK" \
      --bulk-parser pandas
  else
    echo "[ptm_pipeline] skipping InterPro bulk staging: $PTM_INTERPRO_BULK not found"
  fi
fi

if [[ -f "$PTM_INTERPRO_DOMAIN_INTERVALS" && -s "$PTM_INTERPRO_DOMAIN_INTERVALS" ]]; then
  run_step "annotate domain context" python "$PTM_CODE_ROOT/scripts/05_annotate_domain_context.py" \
    --sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
    --position-col corrected_position \
    --outdir "$PTM_RESULTS_ROOT/domain_context"
  run_step "summarize domain enrichment" python "$PTM_CODE_ROOT/scripts/07_summarize_domain_enrichment.py" \
    --annotated-sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
    --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
    --canonical-fasta "$PTM_CANONICAL_FASTA" \
    --ontology "$PTM_CODE_ROOT/config/domain_class_ontology.json" \
    --outdir "$PTM_RESULTS_ROOT/domain_context"
else
  echo "[ptm_pipeline] skipping domain-context annotation: $PTM_INTERPRO_DOMAIN_INTERVALS missing or empty"
fi

if [[ ( -f "$PTM_MOBIDB_INTERVALS" && -s "$PTM_MOBIDB_INTERVALS" ) || ( -f "$PTM_INTERPRO_DOMAIN_INTERVALS" && -s "$PTM_INTERPRO_DOMAIN_INTERVALS" ) ]]; then
  COMPARE_ARGS=(
    --base-master "$PTM_MASTER_TSV"
    --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
    --canonical-fasta "$PTM_CANONICAL_FASTA"
    --outdir "$PTM_RESULTS_ROOT/ptm_compare"
  )
  if [[ -f "$PTM_MOBIDB_INTERVALS" && -s "$PTM_MOBIDB_INTERVALS" ]]; then
    COMPARE_ARGS+=(--disorder-intervals "$PTM_MOBIDB_INTERVALS")
  fi
  if [[ -f "$PTM_INTERPRO_DOMAIN_INTERVALS" && -s "$PTM_INTERPRO_DOMAIN_INTERVALS" ]]; then
    COMPARE_ARGS+=(--interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS")
  fi
  run_step "cross-PTM context comparison" python "$PTM_CODE_ROOT/scripts/09_compare_ptm_contexts.py" "${COMPARE_ARGS[@]}"
fi

if [[ -f "$PTM_CDCODE_STAGED" && -s "$PTM_CDCODE_STAGED" ]]; then
  COND_ARGS=(
    --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
    --staged-membership "$PTM_CDCODE_STAGED"
    --outdir "$PTM_RESULTS_ROOT/condensates"
  )
  if [[ -f "$PTM_CDCODE_PROTEINS" && -s "$PTM_CDCODE_PROTEINS" ]]; then
    COND_ARGS+=(--proteins-table "$PTM_CDCODE_PROTEINS")
  fi
  if [[ -f "$PTM_CDCODE_CONDENSATES" && -s "$PTM_CDCODE_CONDENSATES" ]]; then
    COND_ARGS+=(--condensates-table "$PTM_CDCODE_CONDENSATES")
  fi
  run_step "condensate enrichment" python "$PTM_CODE_ROOT/scripts/12_condensate_enrichment.py" \
    "${COND_ARGS[@]}"
  if [[ -f "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" ]]; then
    run_step "assembly membership clustering" python "$PTM_CODE_ROOT/scripts/15_assembly_membership_clustering.py" \
      --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
      --condensate-enrichment "$PTM_RESULTS_ROOT/condensates/condensate_enrichment.tsv" \
      --condensate-membership "$PTM_RESULTS_ROOT/condensates/condensate_membership.tsv" \
      --functional-labels "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" \
      --outdir "$PTM_RESULTS_ROOT/assembly_clustering"
  fi
fi

if command -v Rscript >/dev/null 2>&1 && [[ -f "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" ]]; then
  export HOME="${TMPDIR:-/tmp}/ptm_r_home_${SLURM_JOB_ID:-manual}"
  export XDG_CACHE_HOME="${TMPDIR:-/tmp}/ptm_r_cache_${SLURM_JOB_ID:-manual}"
  mkdir -p "$HOME" "$XDG_CACHE_HOME"
  run_step "clusterProfiler GO enrichment" Rscript "$PTM_CODE_ROOT/scripts/13_clusterprofiler_residual_go.R" \
    --labels "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" \
    --outdir "$PTM_RESULTS_ROOT/functional_class_union/clusterprofiler_go" \
    --include_all_methyl true
  run_step "disease enrichment" Rscript "$PTM_CODE_ROOT/scripts/14_disease_enrichment.R" \
    --labels "$PTM_RESULTS_ROOT/functional_class_union/protein_functional_labels.tsv" \
    --outdir "$PTM_RESULTS_ROOT/functional_class_union/disease_enrichment"
fi
