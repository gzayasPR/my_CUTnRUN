#!/bin/bash
# Usage: bash submit_cutrun_pipeline.sh <param_file>

PARAM_FILE=$1

if [[ ! -f "$PARAM_FILE" ]]; then
  echo "❌ Parameter file not found: $PARAM_FILE"
  exit 1
fi

# Source parameter file
source "$PARAM_FILE"

# Validate required variables
if [[ -z "${proj_env:-}" || -z "${NAME:-}" ]]; then
  echo "❌ Missing required variables in $PARAM_FILE (expecting proj_env and NAME)"
  exit 1
fi

# Source project environment
source "$proj_env"

if [[ -z "${my_bash:-}" || -z "${my_bin:-}" ]]; then
  echo "❌ proj_env does not define required variables: my_bash, my_bin"
  exit 1
fi

# Create log directory
LOG_DIR="${my_bash}/${NAME}/master_log"
mkdir -p "$LOG_DIR"

echo "📂 Log directory: $LOG_DIR"

# Submit the pipeline and capture job ID
JOBID=$(sbatch --job-name="${NAME} _cutrun" \
  --output="${LOG_DIR}/cutrun_master_%j.out" \
  --error="${LOG_DIR}/cutrun_master_%j.err" \
  "${my_bin}/cutrun_pipe.sh" "$PARAM_FILE" | awk '{print $4}')

# Save summary info
INFO_FILE="${LOG_DIR}/cutrun_master.info"
{
  echo "🧬 CUT&RUN Pipeline Submission Info"
  echo "----------------------------------"
  echo "Project Name : $NAME"
  echo "Job ID       : $JOBID"
  echo "Log Dir      : $LOG_DIR"
  echo "Out File     : cutrun_master_${JOBID}.out"
  echo "Err File     : cutrun_master_${JOBID}.err"
  echo "Param File   : $PARAM_FILE"
  echo "Submitted At : $(date)"
} > "$INFO_FILE"

# Final echo to terminal
echo "✅ Submitted CUT&RUN pipeline for project: $NAME"
echo "🧾 SLURM Job ID: $JOBID"
echo "📂 Log files:    $LOG_DIR"
echo "📝 Info saved to: $INFO_FILE"
