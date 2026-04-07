#!/bin/bash

CONDA_BASE=/gs/gsfs0/hpc01/rhel8/apps/conda3
ENV_FILE=${1:-envs/ptm_pipeline_hpc.yml}
ENV_NAME=${2:-ptm_pipeline}
SOLVER=${PTM_ENV_SOLVER:-conda}

. "$CONDA_BASE/bin/activate"

if [ "$SOLVER" = "mamba" ]; then
  mamba env create -n "$ENV_NAME" -f "$ENV_FILE"
else
  conda env create -n "$ENV_NAME" -f "$ENV_FILE"
fi
