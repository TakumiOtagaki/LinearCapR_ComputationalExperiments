#!/bin/sh

#------ pjsub option --------# 
#PJM -N bprna_auc
#PJM -L rscgrp=short-a
#PJM -L node=1
#PJM -L elapse=2:00:00
#PJM -g gs58
#PJM -j

# Environment setup
conda_env="lincapr"
export python="/work/gs58/s58007/app/anaconda3/envs/$conda_env/bin/python"

# Run AUC analysis for bpRNA dataset
echo "Starting bpRNA AUC analysis..."
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Number of CPUs: $(nproc)"

# Execute AUC analysis with distance thresholds
# PROCESSES=${NUM_CPU:-$(nproc)}
PROCESSES=${NUM_CPU:-5}
$python -m scripts.analysis.bpRNAanalysis.aucroc --processes ${PROCESSES}

echo "AUC analysis completed at: $(date)"

echo "Job completed at: $(date)"
