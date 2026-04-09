#!/bin/bash
#SBATCH -J ptm_condensate
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00

set -eo pipefail

CONDA_BASE=/gs/gsfs0/hpc01/rhel8/apps/conda3
CONDA_EXE="$CONDA_BASE/bin/conda"

run_step() {
  echo "[ptm_condensate] $1"
  shift
  "$@"
  status=$?
  if [ "$status" -ne 0 ]; then
    echo "[ptm_condensate] step failed with exit code $status"
    echo "[ptm_condensate] command: $*"
    exit "$status"
  fi
}

require_file() {
  if [ ! -f "$1" ]; then
    echo "[ptm_condensate] required file not found: $1"
    exit 1
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
  echo "[ptm_condensate] failed to initialize Conda shell hook"
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
  echo "[ptm_condensate] failed to activate Conda environment"
  echo "[ptm_condensate] create one with: bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml $PTM_ENV_NAME"
  exit "$ACTIVATE_STATUS"
fi

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$PTM_RESULTS_ROOT/condensate_enrichment" "$(dirname "$PTM_CDCODE_STAGED")"

require_file "$PTM_BASE_MASTER"
require_file "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv"
require_file "$PTM_CANONICAL_FASTA"

if [ "$PTM_CDCODE_FORCE_RESTAGE" = "1" ] || [ ! -f "$PTM_CDCODE_STAGED" ]; then
  require_file "$PTM_CDCODE_PROTEINS"
  require_file "$PTM_CDCODE_PROTEIN_MAP"
  require_file "$PTM_CDCODE_CONDENSATES"
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

if [ -f "$PTM_DISORDER_INTERVALS" ]; then
  COND_ARGS+=(--disorder-intervals "$PTM_DISORDER_INTERVALS")
fi
if [ "$PTM_CONDENSATE_EXPERIMENTAL_ONLY" = "1" ]; then
  COND_ARGS+=(--experimental-only)
fi
if [ -n "$PTM_CONDENSATE_MIN_CONFIDENCE" ]; then
  COND_ARGS+=(--min-confidence-score "$PTM_CONDENSATE_MIN_CONFIDENCE")
fi

run_step "condensate PTM enrichment" "$PTM_PYTHON_BIN" "$PROJECT_ROOT/scripts/13_condensate_ptm_enrichment.py" "${COND_ARGS[@]}"

echo "[ptm_condensate] completed successfully"
