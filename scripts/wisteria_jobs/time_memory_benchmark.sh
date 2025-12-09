#!/bin/sh

#------ pjsub option --------#
#PJM -N ${JOBNAME}
#PJM -L rscgrp=${RSCGRP}
#PJM -L node=${NODES}
#PJM -L elapse=${ELAPSE}
#PJM -g ${GROUP_ID}
#PJM -j

set -eu

MODE_FROM_LAUNCHER="${MODE}"
DATASET_LABEL_FROM_LAUNCHER="${DATASET_LABEL}"
FASTA_FROM_LAUNCHER="${FASTA}"
BEAMS_FROM_LAUNCHER="${BEAMS}"

: "${MODE_FROM_LAUNCHER:?MODE must be LinCapR or CapR}"
: "${DATASET_LABEL_FROM_LAUNCHER:?DATASET_LABEL is required}"
: "${FASTA_FROM_LAUNCHER:?FASTA path is required}"
: "${BEAMS_FROM_LAUNCHER:?BEAMS list is required}"

MODE="${MODE_FROM_LAUNCHER}"
DATASET_LABEL="${DATASET_LABEL_FROM_LAUNCHER}"
FASTA="${FASTA_FROM_LAUNCHER}"
BEAMS="${BEAMS_FROM_LAUNCHER}"

CONDA_ENV=${CONDA_ENV:-lincapr}
REPO_ROOT=${REPO_ROOT:-/work/gs58/s58007/LinCapR_Experiments}
PYTHON_PREFIX=${PYTHON_PREFIX:-/work/gs58/s58007/app/anaconda3/envs}
PYTHON_BIN=${PYTHON_BIN:-${PYTHON_PREFIX}/${CONDA_ENV}/bin/python}
SCRIPT_MODULE="scripts.analysis.benchmark_lincapr"
JOBS=${JOBS:-8}
MACHINE_TAG=${MACHINE_TAG:-$(hostname)}
RESULT_ROOT=${RESULT_DIR:-${REPO_ROOT}/result/time_memory/${MACHINE_TAG}}
ENERGY_MODEL=${ENERGY_MODEL:-}

mkdir -p "${RESULT_ROOT}" "${RESULT_ROOT}/logs" "${REPO_ROOT}/result/time_memory/tmp"

SPEC_FILE="${RESULT_ROOT}/machine_info_${DATASET_LABEL}_${MODE}.txt"
if [ ! -f "${SPEC_FILE}" ]; then
  {
    echo "# recorded_at: $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    echo "# hostname: $(hostname)"
    echo "# cpuinfo"
    if command -v lscpu >/dev/null 2>&1; then
      lscpu
    else
      echo "lscpu: unavailable"
    fi
    echo "# memory"
    if command -v free >/dev/null 2>&1; then
      free -h
    else
      echo "free: unavailable"
    fi
    echo "# uname"
    uname -a
  } > "${SPEC_FILE}"
fi

if [ "${MODE}" = "LinCapR" ]; then
  EXEC_PATH="${REPO_ROOT}/LinearCapR/LinCapR"
  ENERGY_VALUE=${ENERGY_MODEL:-turner2004}
else
  EXEC_PATH="${REPO_ROOT}/CapR/CapR"
  ENERGY_VALUE=${ENERGY_MODEL:-na}
fi

if [ ! -x "${EXEC_PATH}" ]; then
  echo "Executable not found or not executable: ${EXEC_PATH}" >&2
  exit 1
fi

LOG_FILE="${RESULT_ROOT}/logs/${DATASET_LABEL}_${MODE}_$(date +%Y%m%d%H%M%S).log"
touch "${LOG_FILE}"

# Record the launcher node characteristics for traceability.
if [ -n "${SUBMITTER_MACHINE_HOSTNAME:-}" ] || [ -n "${SUBMITTER_MACHINE_ARCH:-}" ] || [ -n "${SUBMITTER_MACHINE_CPU_COUNT:-}" ] || [ -n "${SUBMITTER_MACHINE_MEMORY_TOTAL:-}" ]; then
  {
    echo "submitter_machine.hostname=${SUBMITTER_MACHINE_HOSTNAME:-unknown}"
    echo "submitter_machine.arch=${SUBMITTER_MACHINE_ARCH:-unknown}"
    echo "submitter_machine.cpu_count=${SUBMITTER_MACHINE_CPU_COUNT:-unknown}"
    echo "submitter_machine.memory_total=${SUBMITTER_MACHINE_MEMORY_TOTAL:-unknown}"
  } | tee -a "${LOG_FILE}"
fi

cd "${REPO_ROOT}"

for beam in ${BEAMS}; do
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] dataset=${DATASET_LABEL} mode=${MODE} beam=${beam}" | tee -a "${LOG_FILE}"
  CMD_TEMPLATE="${EXEC_PATH} {fasta} /dev/null ${beam}"
  if [ "${MODE}" = "LinCapR" ]; then
    CMD_TEMPLATE="${CMD_TEMPLATE} --energy ${ENERGY_VALUE}"
  fi
  OUTPUT_FILE="${RESULT_ROOT}/${DATASET_LABEL}_${MODE}_beam${beam}.csv"
  echo "Running: ${CMD_TEMPLATE}" | tee -a "${LOG_FILE}"
  "${PYTHON_BIN}" -m "${SCRIPT_MODULE}" \
    --fasta "${FASTA}" \
    --cmd "${CMD_TEMPLATE}" \
    --out "${OUTPUT_FILE}" \
    --jobs "${JOBS}" \
    --metadata dataset=${DATASET_LABEL} tool=${MODE} beam=${beam} energy=${ENERGY_VALUE} machine=${MACHINE_TAG} jobs=${JOBS} \
    >> "${LOG_FILE}" 2>&1
done
