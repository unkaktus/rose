#!/bin/bash
# This script is to be called from the batch job,
# so that Slurm environment variables are set after scheduling.

STATE_NAME=${1%.*}
TOTAL_TASK_NUMBER=21
echo "State name: $STATE_NAME"
mkdir -p "$STATE_NAME"

apptainer exec --bind /scratch:/scratch \
    /path/to/rose.sif render_state.py \
    --state=$1 \
    --total-task-number=$TOTAL_TASK_NUMBER \
    --task-id=$(($SLURM_NODEID*$SLURM_NTASKS_PER_NODE + $SLURM_LOCALID)) \
    --output-dir=$STATE_NAME