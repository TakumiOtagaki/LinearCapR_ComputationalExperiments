#!/bin/bash

TEMPLATE="/work/gs58/s58007/LinCapR_Experiments/scripts/wisteria_jobs/sarscov2_ensemble_energy.sh"

JOB_SUFFIX=${1:-default}
JOBNAME="sc2_ensemble_${JOB_SUFFIX}"

export JOBNAME
export JOB_SUFFIX

echo "Submitting job: ${JOBNAME}"
envsubst '${JOBNAME} ${JOB_SUFFIX}' < "${TEMPLATE}" | pjsub
