#!/bin/bash
#SBATCH -J ptm_interpro_api
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00

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
INTERPRO_API_WORKERS=${PTM_INTERPRO_API_WORKERS:-8}
INTERPRO_API_PAGE_SIZE=${PTM_INTERPRO_API_PAGE_SIZE:-200}
INTERPRO_API_TIMEOUT=${PTM_INTERPRO_API_TIMEOUT:-120}
INTERPRO_API_RETRIES=${PTM_INTERPRO_API_RETRIES:-5}

run_step() {
  echo "[ptm_interpro_api] $1"
  shift
  "$@"
  status=$?
  if [ "$status" -ne 0 ]; then
    echo "[ptm_interpro_api] step failed with exit code $status"
    echo "[ptm_interpro_api] command: $*"
    exit "$status"
  fi
}

if [ "$CONDA_HOOK_STATUS" -ne 0 ]; then
  echo "[ptm_interpro_api] failed to initialize Conda shell hook"
  exit "$CONDA_HOOK_STATUS"
fi

if [ -n "$ENV_PREFIX" ]; then
  conda activate "$ENV_PREFIX"
  ACTIVATE_STATUS=$?
else
  conda activate "$ENV_NAME"
  ACTIVATE_STATUS=$?
fi

if [ "$ACTIVATE_STATUS" -ne 0 ]; then
  echo "[ptm_interpro_api] failed to activate Conda environment"
  exit "$ACTIVATE_STATUS"
fi

BASE_DIR=${1:-$PWD}
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$BASE_DIR/data/context"
mkdir -p "$BASE_DIR/results/domain_context"

ACCESSIONS_FILE=${PTM_INTERPRO_ACCESSIONS_FILE:-$BASE_DIR/results/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv}
SITES_FILE=${PTM_INTERPRO_SITES_FILE:-$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv}

if [[ ! -f "$ACCESSIONS_FILE" ]]; then
  echo "[ptm_interpro_api] missing accessions file: $ACCESSIONS_FILE"
  exit 1
fi

if [[ ! -f "$SITES_FILE" ]]; then
  echo "[ptm_interpro_api] missing integrated sites file: $SITES_FILE"
  exit 1
fi

run_step "stage InterPro intervals via API" python "$BASE_DIR/scripts/00_stage_context_sources.py" \
  --outdir "$BASE_DIR/data/context" \
  --accessions "$ACCESSIONS_FILE" \
  --interpro-source api \
  --workers "$INTERPRO_API_WORKERS" \
  --page-size "$INTERPRO_API_PAGE_SIZE" \
  --timeout "$INTERPRO_API_TIMEOUT" \
  --retries "$INTERPRO_API_RETRIES"

if [[ -s "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" ]]; then
  run_step "annotate domain context" python "$BASE_DIR/scripts/05_annotate_domain_context.py" \
    --sites "$SITES_FILE" \
    --interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" \
    --position-col corrected_position \
    --outdir "$BASE_DIR/results/domain_context"

  run_step "summarize domain enrichment" python "$BASE_DIR/scripts/07_summarize_domain_enrichment.py" \
    --annotated-sites "$BASE_DIR/results/domain_context/sites_with_domain_context.tsv" \
    --interpro-intervals "$BASE_DIR/data/context/interpro_human_reviewed_domain_like_intervals.tsv" \
    --canonical-fasta "$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta" \
    --outdir "$BASE_DIR/results/domain_context"
else
  echo "[ptm_interpro_api] no non-empty InterPro domain interval file was produced"
fi
