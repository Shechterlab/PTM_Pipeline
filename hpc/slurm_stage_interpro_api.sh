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
  echo "[ptm_interpro_api] failed to initialize Conda shell hook"
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
  echo "[ptm_interpro_api] failed to activate Conda environment"
  exit "$ACTIVATE_STATUS"
fi

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$PTM_DATA_ROOT/context" "$PTM_RESULTS_ROOT/domain_context"

if [[ ! -f "$PTM_INTERPRO_ACCESSIONS_FILE" ]]; then
  echo "[ptm_interpro_api] missing accessions file: $PTM_INTERPRO_ACCESSIONS_FILE"
  exit 1
fi

if [[ ! -f "$PTM_INTERPRO_SITES_FILE" ]]; then
  echo "[ptm_interpro_api] missing integrated sites file: $PTM_INTERPRO_SITES_FILE"
  exit 1
fi

run_step "stage InterPro intervals via API" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/00_stage_context_sources.py" \
  --outdir "$PTM_DATA_ROOT/context" \
  --accessions "$PTM_INTERPRO_ACCESSIONS_FILE" \
  --interpro-source api \
  --workers "$PTM_INTERPRO_API_WORKERS" \
  --page-size "$PTM_INTERPRO_API_PAGE_SIZE" \
  --timeout "$PTM_INTERPRO_API_TIMEOUT" \
  --retries "$PTM_INTERPRO_API_RETRIES"

echo "[ptm_interpro_api] completed successfully"
