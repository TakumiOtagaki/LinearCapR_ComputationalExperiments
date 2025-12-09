#!/bin/bash

# テンプレートスクリプトのパス
TEMPLATE="/work/gs58/s58007/LinCapR_Experiments/scripts/wisteria_jobs/LinCapR_CapR_on_bpRNA.sh"

# MODE と BS の配列
modes=(LinCapR CapR)
bs_values=(50 100 200 300 400 500)
energy_models_lin=(turner2004 turner1999)
export NUM_CPU=50

for MODE in "${modes[@]}"; do
  if [ "$MODE" = "LinCapR" ]; then
    energies=("${energy_models_lin[@]}")
  else
    energies=("")
  fi

  for BS in "${bs_values[@]}"; do
    for ENERGY_MODEL in "${energies[@]}"; do
      export MODE
      export BS
      export ENERGY_MODEL

      job_suffix="bs${BS}"
      if [ -n "$ENERGY_MODEL" ]; then
        job_suffix+="_${ENERGY_MODEL}"
      fi
      prefix=$(echo "$MODE" | cut -c 1-3)
      export JOBNAME="${prefix}_${job_suffix}"

      echo "Submitting job: $JOBNAME (mode=$MODE, bs=$BS${ENERGY_MODEL:+, energy=$ENERGY_MODEL})"
      envsubst '${JOBNAME} ${MODE} ${BS} ${NUM_CPU} ${ENERGY_MODEL}' < "$TEMPLATE" | pjsub
    done
  done
done
