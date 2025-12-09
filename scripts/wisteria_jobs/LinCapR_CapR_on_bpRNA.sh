#!/bin/sh

#------ pjsub option --------# 
#PJM -N ${JOBNAME}
#PJM -L rscgrp=short-a
#PJM -L node=1
#PJM -L elapse=2:00:00
#PJM -g gs58
#PJM -j

#------- Program execution -------#
# cache="/data/scratch/gs58/s58007/"

# HOW TO USE
# $ ./LinCapR 
# Usage: ./LinCapR <input_file> <output_file> <beam_size> [options]

# ------------ rewrite ------------ #
conda_env="lincapr"
export python="/work/gs58/s58007/app/anaconda3/envs/$conda_env/bin/python"
export python_script_m="scripts.analysis.bpRNAanalysis.run_LinCapR_CapR"

cmd_args="--mode ${MODE} --beam-sizes ${BS} --cpu ${NUM_CPU}"
if [ -n "${ENERGY_MODEL}" ]; then
  cmd_args="${cmd_args} --energy-model ${ENERGY_MODEL}"
fi

echo "$python -m $python_script_m ${cmd_args}"
$python -m $python_script_m ${cmd_args}
