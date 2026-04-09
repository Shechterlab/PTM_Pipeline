#!/bin/bash
#SBATCH -J ptm_condensate
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 08:00:00

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
BASE_DIR=${1:-$PWD}
CDCODE_DIR=${PTM_CDCODE_DIR:-$BASE_DIR/data/external/cd_code}
CDCODE_PROTEINS=${PTM_CDCODE_PROTEINS:-$CDCODE_DIR/proteins_202603181653.csv}
CDCODE_PROTEIN_MAP=${PTM_CDCODE_PROTEIN_MAP:-$CDCODE_DIR/protein2cdcode_v2.2.tsv}
CDCODE_CONDENSATES=${PTM_CDCODE_CONDENSATES:-$CDCODE_DIR/condensates_202603181647.csv}
STAGED_CONDENSATE=${PTM_CONDENSATE_STAGED:-$BASE_DIR/data/context/cd_code_condensate_proteins.tsv}
CANONICAL_FASTA=${PTM_CANONICAL_FASTA:-$BASE_DIR/data/context/uniprot_human_reviewed_canonical.fasta}
BASE_MASTER=${PTM_BASE_MASTER:-$BASE_DIR/data/base/human_ptm_master.tsv}
INTEGRATED_ARG=${PTM_INTEGRATED_ARG:-$BASE_DIR/results/integrated/human_arg_methyl_union_dedup_by_site.tsv}
OUTDIR=${PTM_CONDENSATE_OUTDIR:-$BASE_DIR/results/condensate_enrichment}
DISORDER_INTERVALS=${PTM_DISORDER_INTERVALS:-$BASE_DIR/data/context/mobidb_human_reviewed_disorder_intervals.tsv}
EXPERIMENTAL_ONLY=${PTM_CONDENSATE_EXPERIMENTAL_ONLY:-0}
MIN_CONFIDENCE=${PTM_CONDENSATE_MIN_CONFIDENCE:-}
FORCE_RESTAGE=${PTM_CDCODE_FORCE_RESTAGE:-0}

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

if [ "$CONDA_HOOK_STATUS" -ne 0 ]; then
  echo "[ptm_condensate] failed to initialize Conda shell hook"
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
  echo "[ptm_condensate] failed to activate Conda environment"
  echo "[ptm_condensate] create one with: bash hpc/create_ptm_env.sh envs/ptm_pipeline_hpc.yml $ENV_NAME"
  exit "$ACTIVATE_STATUS"
fi

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID:-manual}"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$BASE_DIR/data/context" "$OUTDIR"

require_file "$BASE_MASTER"
require_file "$INTEGRATED_ARG"
require_file "$CANONICAL_FASTA"

if [ "$FORCE_RESTAGE" = "1" ] || [ ! -f "$STAGED_CONDENSATE" ]; then
  require_file "$CDCODE_PROTEINS"
  require_file "$CDCODE_PROTEIN_MAP"
  require_file "$CDCODE_CONDENSATES"
  run_step "stage CD-CODE condensate proteins" python "$BASE_DIR/scripts/00_stage_cd_code_condensates.py" \
    --cdcode-proteins "$CDCODE_PROTEINS" \
    --cdcode-protein-map "$CDCODE_PROTEIN_MAP" \
    --cdcode-condensates "$CDCODE_CONDENSATES" \
    --out "$STAGED_CONDENSATE"
else
  echo "[ptm_condensate] using staged condensate table: $STAGED_CONDENSATE"
fi

ANALYSIS_ARGS=(
  --base-master "$BASE_MASTER"
  --integrated-arg-sites "$INTEGRATED_ARG"
  --canonical-fasta "$CANONICAL_FASTA"
  --condensate-proteins "$STAGED_CONDENSATE"
  --outdir "$OUTDIR"
)

if [ -f "$DISORDER_INTERVALS" ]; then
  ANALYSIS_ARGS+=(--disorder-intervals "$DISORDER_INTERVALS")
else
  echo "[ptm_condensate] disorder intervals not found, running without disorder covariate file: $DISORDER_INTERVALS"
fi

if [ "$EXPERIMENTAL_ONLY" = "1" ]; then
  ANALYSIS_ARGS+=(--experimental-only)
fi

if [ -n "$MIN_CONFIDENCE" ]; then
  ANALYSIS_ARGS+=(--min-confidence-score "$MIN_CONFIDENCE")
fi

run_step "run condensate enrichment" python "$BASE_DIR/scripts/13_condensate_ptm_enrichment.py" "${ANALYSIS_ARGS[@]}"

echo "[ptm_condensate] complete"
