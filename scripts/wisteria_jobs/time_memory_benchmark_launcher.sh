#!/bin/bash
set -euo pipefail

REPO_ROOT=${REPO_ROOT:-/work/gs58/s58007/LinCapR_Experiments}
TEMPLATE="${REPO_ROOT}/scripts/wisteria_jobs/time_memory_benchmark.sh"

if [ ! -f "${TEMPLATE}" ]; then
  echo "Template not found: ${TEMPLATE}" >&2
  exit 1
fi

declare -A FASTA_MAP=(
  [bprna]="${REPO_ROOT}/data/processed/bprna/bpRNA_multifasta.fasta"
  [rnacentral]="${REPO_ROOT}/data/processed/rnacentral/rnacentral_active_20.fasta"
)

# datasets=(bprna rnacentral) # 10/15 あたりに bprna は終わった。
datasets=rnacentral
modes=(LinCapR CapR)

# export BEAMS=${BEAMS:-"50 100 200 300 400 500"} # 10/15 あたりに 50, 100 は終えた。
# export BEAMS=${BEAMS:-"300 400"} # 10/21 15:20
export BEAMS=${BEAMS:-"200 500"} # 今からやる 10/21 15:25 ~ 
export JOBS=${JOBS:-32}
export RSCGRP=${RSCGRP:-regular-a}
export NODES=${NODES:-1}
export ELAPSE=${ELAPSE:-06:00:00}
export GROUP_ID=${GROUP_ID:-gs58}
export CONDA_ENV=${CONDA_ENV:-lincapr}
export PYTHON_PREFIX=${PYTHON_PREFIX:-/work/gs58/s58007/app/anaconda3/envs}
export RUN_TAG=${RUN_TAG:-$(date +%Y%m%d_%H%M)}
export MACHINE_TAG_BASE=${MACHINE_TAG_BASE:-wisteria}
export MACHINE_TAG="${MACHINE_TAG_BASE}_${RUN_TAG}"
export ENERGY_MODEL_LINCAPR=${ENERGY_MODEL_LINCAPR:-turner2004}

# Collect submitter node specs so each job can log them.
SUBMITTER_MACHINE_HOSTNAME=${SUBMITTER_MACHINE_HOSTNAME:-$(hostname 2>/dev/null || echo unknown)}
SUBMITTER_MACHINE_ARCH=${SUBMITTER_MACHINE_ARCH:-$(uname -m 2>/dev/null || echo unknown)}
if command -v nproc >/dev/null 2>&1; then
  SUBMITTER_MACHINE_CPU_COUNT=${SUBMITTER_MACHINE_CPU_COUNT:-$(nproc 2>/dev/null || echo "")}
else
  SUBMITTER_MACHINE_CPU_COUNT=${SUBMITTER_MACHINE_CPU_COUNT:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo "")}
fi
if [ -z "${SUBMITTER_MACHINE_CPU_COUNT}" ]; then
  SUBMITTER_MACHINE_CPU_COUNT=unknown
fi
if command -v free >/dev/null 2>&1; then
  SUBMITTER_MACHINE_MEMORY_TOTAL=${SUBMITTER_MACHINE_MEMORY_TOTAL:-$(free -h 2>/dev/null | awk '/^Mem:/ {print $2}' || echo "")}
elif [ -r /proc/meminfo ]; then
  SUBMITTER_MACHINE_MEMORY_TOTAL=${SUBMITTER_MACHINE_MEMORY_TOTAL:-$(awk '/MemTotal/ {printf "%.2fGiB", $2/1024/1024}' /proc/meminfo 2>/dev/null || echo "")}
else
  SUBMITTER_MACHINE_MEMORY_TOTAL=${SUBMITTER_MACHINE_MEMORY_TOTAL:-unknown}
fi
if [ -z "${SUBMITTER_MACHINE_MEMORY_TOTAL}" ]; then
  SUBMITTER_MACHINE_MEMORY_TOTAL=unknown
fi
export SUBMITTER_MACHINE_HOSTNAME SUBMITTER_MACHINE_ARCH SUBMITTER_MACHINE_CPU_COUNT SUBMITTER_MACHINE_MEMORY_TOTAL

job_index=0

for dataset in "${datasets[@]}"; do
  fasta_path="${FASTA_MAP[$dataset]}"
  if [ ! -f "${fasta_path}" ]; then
    echo "FASTA not found for dataset ${dataset}: ${fasta_path}" >&2
    exit 1
  fi
  for mode in "${modes[@]}"; do
    job_index=$((job_index + 1))
    export JOBNAME="tm_${dataset}_${mode}_${job_index}"
    export RSCGRP NODES ELAPSE GROUP_ID CONDA_ENV REPO_ROOT PYTHON_PREFIX
    export MODE="${mode}"
    export DATASET_LABEL="${dataset}"
    export FASTA="${fasta_path}"
    export BEAMS JOBS MACHINE_TAG
    unset RESULT_DIR
    if [ "${mode}" = "LinCapR" ]; then
      export ENERGY_MODEL="${ENERGY_MODEL_LINCAPR}"
    else
      unset ENERGY_MODEL
    fi
    echo "Submitting job ${JOBNAME} (dataset=${dataset}, mode=${mode})"
    envsubst '${JOBNAME} ${RSCGRP} ${NODES} ${ELAPSE} ${GROUP_ID} ${MODE} ${DATASET_LABEL} ${FASTA} ${BEAMS} ${JOBS} ${MACHINE_TAG} ${RESULT_DIR} ${CONDA_ENV} ${REPO_ROOT} ${PYTHON_PREFIX} ${ENERGY_MODEL} ${SUBMITTER_MACHINE_HOSTNAME} ${SUBMITTER_MACHINE_ARCH} ${SUBMITTER_MACHINE_CPU_COUNT} ${SUBMITTER_MACHINE_MEMORY_TOTAL}' < "${TEMPLATE}" | pjsub
  done
done
