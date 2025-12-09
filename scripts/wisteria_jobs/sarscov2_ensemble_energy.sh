#!/bin/sh

#------ pjsub option --------#
#PJM -N ${JOBNAME}
#PJM -L rscgrp=short-a
#PJM -L node=1
#PJM -L elapse=2:00:00
#PJM -g gs58
#PJM -j

#------- Program execution -------#
conda_env="lincapr"
python_bin="/work/gs58/s58007/app/anaconda3/envs/${conda_env}/bin/python"
script_path="/work/gs58/s58007/LinCapR_Experiments/scripts/analysis/sarscov2_ensemble_energy.py"

FASTA="/work/gs58/s58007/LinCapR_Experiments/data/processed/sarscov2/NC_045512.RNA.fa"
OUTPUT_DIR="/work/gs58/s58007/LinCapR_Experiments/result/sarscov2/ensemble_energy/${JOB_SUFFIX}"
BEAMS="50 100 200 300 400 500"
ENERGY_MODELS="turner2004 turner1999"

mkdir -p "${OUTPUT_DIR}"
LOG_PATH="${OUTPUT_DIR}/ensemble_energy.log"
TIME_PATH="${OUTPUT_DIR}/time_report.txt"

CMD="${python_bin} ${script_path} --fasta ${FASTA} --beams ${BEAMS} --energy-models ${ENERGY_MODELS} --output-dir ${OUTPUT_DIR}"

echo "Running: ${CMD}" | tee "${LOG_PATH}"
/usr/bin/time -v ${CMD} >> "${LOG_PATH}" 2>> "${TIME_PATH}"
