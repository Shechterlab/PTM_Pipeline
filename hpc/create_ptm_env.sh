#!/bin/bash

CONDA_BASE=/gs/gsfs0/hpc01/rhel8/apps/conda3
CONDA_EXE="$CONDA_BASE/bin/conda"
ENV_FILE=${1:-envs/ptm_pipeline_hpc.yml}
ENV_NAME=${2:-ptm_pipeline}
SOLVER=${PTM_ENV_SOLVER:-conda}

if [ -x "$CONDA_EXE" ]; then
  eval "$("$CONDA_EXE" shell.bash hook)"
else
  . "$CONDA_BASE/bin/activate"
fi

if [ "$SOLVER" = "mamba" ]; then
  mamba env create -n "$ENV_NAME" -f "$ENV_FILE"
else
  conda env create -n "$ENV_NAME" -f "$ENV_FILE"
fi
