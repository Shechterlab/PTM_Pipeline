#!/bin/bash
#SBATCH -J ptm_interpro_api
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 12:00:00

set -eo pipefail

CONDA_BASE=/gs/gsfs0/hpc01/rhel8/apps/conda3
CONDA_EXE="$CONDA_BASE/bin/conda"

if [ -x "$CONDA_EXE" ]; then
  eval "$("$CONDA_EXE" shell.bash hook)"
else
  . "$CONDA_BASE/bin/activate"
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

if [ -n "$ENV_PREFIX" ]; then
  conda activate "$ENV_PREFIX"
else
  conda activate "$ENV_NAME"
fi

BASE_DIR=${1:-$PWD}
. "$BASE_DIR/hpc/ptm_env.sh"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR" "$PTM_CONTEXT_ROOT" "$PTM_RESULTS_ROOT/domain_context"

ACCESSIONS_FILE=${PTM_INTERPRO_ACCESSIONS_FILE:-$PTM_RESULTS_ROOT/remapped/human_arg_methyl_source_rows_remapped_resolved.tsv}
SITES_FILE=${PTM_INTERPRO_SITES_FILE:-$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv}

if [[ ! -f "$ACCESSIONS_FILE" || ! -s "$ACCESSIONS_FILE" ]]; then
  echo "[ptm_interpro_api] missing accessions file: $ACCESSIONS_FILE"
  exit 1
fi

if [[ ! -f "$SITES_FILE" || ! -s "$SITES_FILE" ]]; then
  echo "[ptm_interpro_api] missing integrated sites file: $SITES_FILE"
  exit 1
fi

run_step "repo hygiene guard" bash "$PTM_CODE_ROOT/hpc/check_repo_hygiene.sh" "$PTM_CODE_ROOT"
run_step "stage InterPro intervals via API" python "$PTM_CODE_ROOT/scripts/00_stage_context_sources.py" \
  --outdir "$PTM_CONTEXT_ROOT" \
  --accessions "$ACCESSIONS_FILE" \
  --interpro-source api \
  --workers "$INTERPRO_API_WORKERS" \
  --page-size "$INTERPRO_API_PAGE_SIZE" \
  --timeout "$INTERPRO_API_TIMEOUT" \
  --retries "$INTERPRO_API_RETRIES"

run_step "annotate domain context" python "$PTM_CODE_ROOT/scripts/05_annotate_domain_context.py" \
  --sites "$SITES_FILE" \
  --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --position-col corrected_position \
  --outdir "$PTM_RESULTS_ROOT/domain_context"

run_step "summarize domain enrichment" python "$PTM_CODE_ROOT/scripts/07_summarize_domain_enrichment.py" \
  --annotated-sites "$PTM_RESULTS_ROOT/domain_context/sites_with_domain_context.tsv" \
  --interpro-intervals "$PTM_INTERPRO_DOMAIN_INTERVALS" \
  --canonical-fasta "$PTM_CANONICAL_FASTA" \
  --outdir "$PTM_RESULTS_ROOT/domain_context"
