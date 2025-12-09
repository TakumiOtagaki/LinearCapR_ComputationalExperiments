#!/bin/sh

#------ pjsub option --------#
#PJM -N lincapr_figs
#PJM -L rscgrp=short-a
#PJM -L node=1
#PJM -L elapse=2:00:00
#PJM -g gs58
#PJM -j

# Environment setup
conda_env="lincapr"
export python="/work/gs58/s58007/app/anaconda3/envs/$conda_env/bin/python"

CONFIG_PATH=${1:-graph_config.json}

if [ ! -f "$CONFIG_PATH" ]; then
  echo "[error] Config file not found: $CONFIG_PATH"
  exit 1
fi

echo "Starting publication figure generation..."
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Python interpreter: $python"
echo "Config: $CONFIG_PATH"

$python -m scripts.generate_publication_figs --config "$CONFIG_PATH"
status=$?

echo "Job finished at: $(date)"

exit $status
