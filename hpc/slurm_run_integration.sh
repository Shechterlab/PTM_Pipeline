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

PROJECT_ROOT="${1:-$PWD}"
cd "$PROJECT_ROOT" || exit 1

if [ -x "$CONDA_EXE" ]; then
  eval "$("$CONDA_EXE" shell.bash hook)"
  CONDA_HOOK_STATUS=$?
else
  . "$CONDA_BASE/bin/activate"
  CONDA_HOOK_STATUS=$?
fi

if [ "$CONDA_HOOK_STATUS" -ne 0 ]; then
  echo "[ptm_pipeline] failed to initialize Conda shell hook"
  exit "$CONDA_HOOK_STATUS"
fi

source hpc/ptm_env.sh

if [ -n "$PTM_ENV_PREFIX" ]; then
  conda activate "$PTM_ENV_PREFIX"
  ACTIVATE_STATUS=$?
else
  conda activate "$PTM_ENV_NAME"
  ACTIVATE_STATUS=$?
fi

if [ "$ACTIVATE_STATUS" -ne 0 ]; then
  echo "[ptm_pipeline] failed to activate Conda environment"
  echo "[ptm_pipeline] create one with: bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml $PTM_ENV_NAME"
  exit "$ACTIVATE_STATUS"
fi

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$PTM_DATA_ROOT/context" "$PTM_CDCODE_DIR"
mkdir -p \
  "$PTM_RESULTS_ROOT/parsed_maron" \
  "$PTM_RESULTS_ROOT/parsed_prometheus" \
  "$PTM_RESULTS_ROOT/remapped" \
  "$PTM_RESULTS_ROOT/integrated" \
  "$PTM_RESULTS_ROOT/functional_class_union" \
  "$PTM_RESULTS_ROOT/example_figures" \
  "$PTM_RESULTS_ROOT/motif_arg_odds" \
  "$PTM_RESULTS_ROOT/methyl_arg_clustering" \
  "$PTM_RESULTS_ROOT/ptm_clustering" \
  "$PTM_RESULTS_ROOT/cross_ptm_neighbors" \
  "$PTM_RESULTS_ROOT/disorder_context" \
  "$PTM_RESULTS_ROOT/domain_context" \
  "$PTM_RESULTS_ROOT/ptm_compare" \
  "$PTM_RESULTS_ROOT/review_audit" \
  "$PTM_RESULTS_ROOT/condensate_enrichment"

echo "[ptm_pipeline] code root: $PROJECT_ROOT"
echo "[ptm_pipeline] data root: $PTM_DATA_ROOT"
echo "[ptm_pipeline] results root: $PTM_RESULTS_ROOT"

run_step "stage UniProt context" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/00_stage_context_sources.py" \
  --outdir "$PTM_DATA_ROOT/context" \
  --skip-interpro
run_step "parse Maron S5" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/01_parse_maron_s5.py" \
  --input "$PTM_MARON_XLSX" \
  --exclude-pattern "$PTM_MARON_EXCLUDE_PATTERN" \
  --outdir "$PTM_RESULTS_ROOT/parsed_maron"
run_step "parse ProMetheus mmc2" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/01_parse_prometheus_mmc2.py" \
  --input "$PTM_PROMETHEUS_XLSX" \
  --outdir "$PTM_RESULTS_ROOT/parsed_prometheus"
run_step "remap to canonical coordinates" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/02_remap_arg_methyl_to_canonical.py" \
  --base-master "$PTM_BASE_MASTER" \
  --maron "$PTM_RESULTS_ROOT/parsed_maron/maron_s5_arg_methyl_sites.tsv" \
  --prometheus "$PTM_RESULTS_ROOT/parsed_prometheus/prometheus_mmc2_arg_methyl_sites.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --uniprot-metadata "$PTM_DATA_ROOT/context/uniprot_human_reviewed_canonical_metadata.tsv" \
  --outdir "$PTM_RESULTS_ROOT/remapped"
run_step "integrate remapped sites" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/02_integrate_arg_methyl_sources.py" \
  --premerged-remapped "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
  --outdir "$PTM_RESULTS_ROOT/integrated"

if [[ ! -s "$PTM_INTERPRO_INTERVALS" ]] && [[ "$PTM_INTERPRO_MODE" = "api" ]]; then
  run_step "stage InterPro domain intervals via API" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/00_stage_context_sources.py" \
    --outdir "$PTM_DATA_ROOT/context" \
    --accessions "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
    --interpro-source api \
    --workers "$PTM_INTERPRO_API_WORKERS" \
    --page-size "$PTM_INTERPRO_API_PAGE_SIZE" \
    --timeout "$PTM_INTERPRO_API_TIMEOUT" \
    --retries "$PTM_INTERPRO_API_RETRIES"
elif [[ ! -s "$PTM_INTERPRO_INTERVALS" ]] && [[ "$PTM_INTERPRO_MODE" = "bulk" || ( "$PTM_INTERPRO_MODE" = "auto" && -f "$PTM_INTERPRO_PROTEIN2IPR" ) ]]; then
  if [[ -f "$PTM_INTERPRO_PROTEIN2IPR" ]]; then
    run_step "stage InterPro domain intervals" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/00_stage_context_sources.py" \
      --outdir "$PTM_DATA_ROOT/context" \
      --accessions "$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv" \
      --interpro-source bulk \
      --interpro-protein2ipr "$PTM_INTERPRO_PROTEIN2IPR" \
      --bulk-parser pandas
  else
    echo "[ptm_pipeline] skipping InterPro bulk staging: $PTM_INTERPRO_PROTEIN2IPR not found"
  fi
fi

FUNCTIONAL_ARGS=(
  --base-master "$PTM_BASE_MASTER"
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
  --ontology "$PROJECT_ROOT/config/functional_ontology.json"
  --outdir "$PTM_RESULTS_ROOT/functional_class_union"
)
if [[ -s "$PTM_INTERPRO_INTERVALS" ]]; then
  FUNCTIONAL_ARGS+=(--interpro-intervals "$PTM_INTERPRO_INTERVALS")
fi
run_step "functional enrichment" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/03_functional_class_enrichment.py" "${FUNCTIONAL_ARGS[@]}"

run_step "example figures" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/04_example_figures.py" \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --outdir "$PTM_RESULTS_ROOT/example_figures"

MOTIF_ARGS=(
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
  --canonical-fasta "$PTM_CANONICAL_FASTA"
  --outdir "$PTM_RESULTS_ROOT/motif_arg_odds"
)
if [[ -f "$PTM_DISORDER_INTERVALS" ]]; then
  MOTIF_ARGS+=(--disorder-intervals "$PTM_DISORDER_INTERVALS")
fi
run_step "motif enrichment and protein shrinkage model" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/08_motif_and_arg_odds.py" "${MOTIF_ARGS[@]}"

run_step "methyl-Arg clustering" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/10_methyl_arg_clustering.py" \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --permutations "$PTM_METHYL_CLUSTER_PERMUTATIONS" \
  --progress-every "$PTM_METHYL_CLUSTER_PROGRESS_EVERY" \
  --outdir "$PTM_RESULTS_ROOT/methyl_arg_clustering"

run_step "cross-PTM clustering" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/11_compare_ptm_clustering.py" \
  --base-master "$PTM_BASE_MASTER" \
  --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/ptm_clustering"

if [ "$PTM_RUN_CROSS_PTM_NEIGHBORS" = "1" ]; then
  CROSS_NEIGHBOR_ARGS=(
    --base-master "$PTM_BASE_MASTER"
    --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
    --canonical-fasta "$PTM_CANONICAL_FASTA"
    --outdir "$PTM_RESULTS_ROOT/cross_ptm_neighbors"
  )
  if [[ -f "$PTM_DISORDER_INTERVALS" ]]; then
    CROSS_NEIGHBOR_ARGS+=(--disorder-intervals "$PTM_DISORDER_INTERVALS")
  fi
  run_step "cross-PTM neighbor enrichment" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/12_cross_ptm_neighbor_enrichment.py" "${CROSS_NEIGHBOR_ARGS[@]}"
fi

if [[ -f "$PTM_DISORDER_INTERVALS" ]]; then
  run_step "annotate disorder context" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/06_annotate_disorder_context.py" \
    --sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --disorder-intervals "$PTM_DISORDER_INTERVALS" \
    --canonical-fasta "$PTM_CANONICAL_FASTA" \
    --outdir "$PTM_RESULTS_ROOT/disorder_context"
else
  echo "[ptm_pipeline] skipping disorder-context staging: $PTM_DISORDER_INTERVALS not found"
fi

if [[ -s "$PTM_INTERPRO_INTERVALS" ]]; then
  run_step "annotate domain context" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/05_annotate_domain_context.py" \
    --sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --interpro-intervals "$PTM_INTERPRO_INTERVALS" \
    --position-col corrected_position \
    --outdir "$PTM_RESULTS_ROOT/domain_context"
  run_step "summarize domain enrichment" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/07_summarize_domain_enrichment.py" \
    --annotated-sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
    --interpro-intervals "$PTM_INTERPRO_INTERVALS" \
    --canonical-fasta "$PTM_CANONICAL_FASTA" \
    --outdir "$PTM_RESULTS_ROOT/domain_context"
else
  echo "[ptm_pipeline] skipping domain-context annotation: $PTM_INTERPRO_INTERVALS not found"
fi

if [[ -f "$PTM_DISORDER_INTERVALS" || -s "$PTM_INTERPRO_INTERVALS" ]]; then
  COMPARE_ARGS=(
    --base-master "$PTM_BASE_MASTER"
    --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
    --canonical-fasta "$PTM_CANONICAL_FASTA"
    --outdir "$PTM_RESULTS_ROOT/ptm_compare"
  )
  if [[ -f "$PTM_DISORDER_INTERVALS" ]]; then
    COMPARE_ARGS+=(--disorder-intervals "$PTM_DISORDER_INTERVALS")
  fi
  if [[ -s "$PTM_INTERPRO_INTERVALS" ]]; then
    COMPARE_ARGS+=(--interpro-intervals "$PTM_INTERPRO_INTERVALS")
  fi
  run_step "cross-PTM context comparison" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/09_compare_ptm_contexts.py" "${COMPARE_ARGS[@]}"
fi

if [ "$PTM_RUN_REVIEW_AUDIT" = "1" ] && [[ -f "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_all_rows.tsv" ]] && [[ -f "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" ]]; then
  run_step "review audit figures" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/14_union_source_audit_and_review_figures.py" \
    --all-rows "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_all_rows.tsv" \
    --dedup-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
    --outdir "$PTM_RESULTS_ROOT/review_audit"
fi

if [ "$PTM_RUN_CONDENSATE" = "1" ] || [ "$PTM_RUN_CONDENSATE" = "true" ] || [ "$PTM_RUN_CONDENSATE" = "yes" ] || [ "$PTM_RUN_CONDENSATE" = "auto" ]; then
  HAVE_STAGE=0
  HAVE_RAW=0
  if [[ -f "$PTM_CDCODE_STAGED" ]]; then
    HAVE_STAGE=1
  fi
  if [[ -f "$PTM_CDCODE_PROTEINS" ]] && [[ -f "$PTM_CDCODE_PROTEIN_MAP" ]] && [[ -f "$PTM_CDCODE_CONDENSATES" ]]; then
    HAVE_RAW=1
  fi
  if [ "$PTM_RUN_CONDENSATE" = "1" ] || [ "$PTM_RUN_CONDENSATE" = "true" ] || [ "$PTM_RUN_CONDENSATE" = "yes" ] || [ "$HAVE_STAGE" -eq 1 ] || [ "$HAVE_RAW" -eq 1 ]; then
    if [ "$HAVE_STAGE" -ne 1 ]; then
      run_step "stage CD-CODE condensates" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/00_stage_cd_code_condensates.py" \
        --cdcode-proteins "$PTM_CDCODE_PROTEINS" \
        --cdcode-protein-map "$PTM_CDCODE_PROTEIN_MAP" \
        --cdcode-condensates "$PTM_CDCODE_CONDENSATES" \
        --out "$PTM_CDCODE_STAGED"
    fi
    COND_ARGS=(
      --base-master "$PTM_BASE_MASTER"
      --integrated-arg-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
      --canonical-fasta "$PTM_CANONICAL_FASTA"
      --condensate-proteins "$PTM_CDCODE_STAGED"
      --outdir "$PTM_RESULTS_ROOT/condensate_enrichment"
    )
    if [[ -f "$PTM_DISORDER_INTERVALS" ]]; then
      COND_ARGS+=(--disorder-intervals "$PTM_DISORDER_INTERVALS")
    fi
    if [ "$PTM_CONDENSATE_EXPERIMENTAL_ONLY" = "1" ]; then
      COND_ARGS+=(--experimental-only)
    fi
    if [ -n "$PTM_CONDENSATE_MIN_CONFIDENCE" ]; then
      COND_ARGS+=(--min-confidence-score "$PTM_CONDENSATE_MIN_CONFIDENCE")
    fi
    run_step "condensate PTM enrichment" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/13_condensate_ptm_enrichment.py" "${COND_ARGS[@]}"
  else
    echo "[ptm_pipeline] skipping condensate enrichment: no staged table or raw CD-CODE files found under $PTM_CDCODE_DIR"
  fi
fi

echo "[ptm_pipeline] completed successfully"
