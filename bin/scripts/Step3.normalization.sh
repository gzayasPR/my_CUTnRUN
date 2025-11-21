#!/bin/bash

PARAM_FILE=$1
# PARAM_FILE="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/parameters/project_env.parameters.sh"
source "${PARAM_FILE}"
# ---- VALIDATE REQUIRED PARAMS ----
REQUIRED_VARS=(NAME CSV account STEP3_CPUS MOUSE_GTF)
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var:-}" ]; then
    echo "❌ Missing required parameter: $var"
    exit 1
  fi
done

#source "${project_env}"
source "${my_bin}/scripts/shell.functions.sh"
ALIGN_DIR="${my_results}/${NAME}/2.Single_Alignment"
NORM_DIR="${my_results}/${NAME}/3.Normalization"
LOG_DIR="${my_bash}/${NAME}/Step3"
BIN_DIR="${my_bin}/scripts/3.Normalize"

mkdir -p "$NORM_DIR" "$LOG_DIR"


echo "Starting Step 3: spike-in normalization and scaling finalized."

scale_script="${BIN_DIR}/calculate.scale.R"
out_scale="${NORM_DIR}/scale_factors.tsv"
touch "$out_scale"
sleep 20
alignment_results="${ALIGN_DIR}/alignment_summary.tsv"
ml R
Rscript "$scale_script" "$alignment_results" "$out_scale"
sleep 20
# --- STEP 4: Normalize and generate scaled bigWigs ---
echo "📈 Submitting normalization jobs..."

norm_job_ids=()

# Avoid subshell to preserve norm_job_ids
while read -r sample; do
  jobid=$(sbatch --account="$account" \
    --job-name="norm_${sample}" \
    --output="${LOG_DIR}/norm_${sample}.out" \
    --error="${LOG_DIR}/norm_${sample}.err" \
    --ntasks=1 --cpus-per-task="$STEP3_CPUS" --mem-per-cpu=8gb \
    "${BIN_DIR}/normalize_signal.sh" "$sample" "$ALIGN_DIR" "$NORM_DIR" "$out_scale" "$STEP3_CPUS" | awk '{print $4}')
  echo "✅ Submitted normalization job for $sample: $jobid"
  norm_job_ids+=("$jobid")
  sleep 1
done < <(tail -n +2 "$CSV" | cut -d',' -f1 | sort | uniq )

echo "🕒 Waiting for normalization jobs to finish..."
wait_for_jobs_array norm_job_ids 1 30

qcid=$(sbatch --account="$account" \
    --job-name="QC_Samples" \
    --output="${LOG_DIR}/QC_Samples.out" \
    --error="${LOG_DIR}/QC_Samples.err" \
    --ntasks=1 --cpus-per-task="$STEP3_CPUS" --mem-per-cpu=8gb \
    "${BIN_DIR}/samples_compare.sh" "$PARAM_FILE" "$CSV" "$ALIGN_DIR" "$NORM_DIR" "$STEP3_CPUS" | awk '{print $4}')
  
wait_for_job_completion_less_output "$qcid" 1 30

# --- STEP 5: Run QC for fragment size and TSS enrichment ---
echo "🧪 Submitting QC jobs..."
QC_DIR="${my_results}/${NAME}/3.Normalization"
mkdir -p "$QC_DIR"
qc_job_ids=()

while read -r sample; do
  jobid=$(sbatch --account="$account" \
    --job-name="qc_${sample}" \
    --output="${LOG_DIR}/qc_${sample}.out" \
    --error="${LOG_DIR}/qc_${sample}.err" \
    --ntasks=1 --cpus-per-task="$STEP3_CPUS" --mem-per-cpu=8gb \
    "${BIN_DIR}/qc_normalization.sh" "$sample" "$ALIGN_DIR" "$NORM_DIR" "$QC_DIR" "$MOUSE_GTF" "$STEP3_CPUS" | awk '{print $4}')
  echo "✅ Submitted QC job for $sample: $jobid"
  qc_job_ids+=("$jobid")
  sleep 1
done < <(tail -n +2 "$CSV" | cut -d',' -f1 | sort | uniq )

echo "🕒 Waiting for  QC jobs to finish..."
wait_for_jobs_array qc_job_ids 1 30

echo "✅ Step 3 complete: spike-in normalization and scaling finalized."
