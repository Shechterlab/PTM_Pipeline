#!/bin/bash
#SBATCH -J ptm_condensate
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 04:00:00

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

if [ -n "$ENV_PREFIX" ]; then
  conda activate "$ENV_PREFIX"
else
  conda activate "$ENV_NAME"
fi

BASE_DIR=${1:-$PWD}
. "$BASE_DIR/hpc/ptm_env.sh"

if [[ ! -f "$PTM_CDCODE_STAGED" || ! -s "$PTM_CDCODE_STAGED" ]]; then
  echo "[ptm_condensate] missing staged condensate membership: $PTM_CDCODE_STAGED"
  exit 1
fi

bash "$PTM_CODE_ROOT/hpc/check_repo_hygiene.sh" "$PTM_CODE_ROOT"
python "$PTM_CODE_ROOT/scripts/12_condensate_enrichment.py" \
  --integrated-sites "$PTM_RESULTS_ROOT/integrated/human_arg_methyl_union_dedup_by_site.tsv" \
  --staged-membership "$PTM_CDCODE_STAGED" \
  --outdir "$PTM_RESULTS_ROOT/condensates"
