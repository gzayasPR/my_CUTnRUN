#!/bin/bash
# --- Arguments ---

PARAM_FILE=$1

source "${PARAM_FILE}"
# ---- VALIDATE REQUIRED PARAMS ----
REQUIRED_VARS=(NAME CSV account STEP1_CPUS)
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var:-}" ]; then
    echo "❌ Missing required parameter: $var"
    exit 1
  fi
done

# --- Load environment + functions ---
source "${my_bin}/scripts/shell.functions.sh"

# --- Define paths ---
OUT_DIR="${my_results}/${NAME}/1.KMET_QC"
KMET_DIR="${OUT_DIR}/kmet_analysis"
TRIM_DIR="${OUT_DIR}/trimmed"
LOG_DIR="${my_bash}/${NAME}/Step1"
BIN_DIR="${my_bin}/scripts/1.Kmet_Analysis_QC"

mkdir -p "$OUT_DIR" "$KMET_DIR" "$TRIM_DIR" "$LOG_DIR/kmet" "$LOG_DIR/trim"

echo "✅ Step1: KMET + Trimming + MultiQC for $NAME"

# --- Print parameters ---
echo "📌 Parameters:"
echo "  Project env: $project_env"
echo "  Metadata CSV: $CSV"
echo "  Output name: $NAME"
echo "  KMET_DIR: $KMET_DIR"
echo "  TRIM_DIR: $TRIM_DIR"
echo "  LOG_DIR: $LOG_DIR"

# --- Store job IDs ---
kmet_job_ids=()
trim_job_ids=()
# --- Track how many times a sample appears ---
declare -A sample_counts

while IFS=',' read -r sample run_name condition fastq_r1 fastq_r2; do
  prefix="${fastq_r1%_R1_001.fastq.gz}"

  # Count occurrences and assign unique sample-run label
  sample_counts["$sample"]=$((sample_counts["$sample"] + 1))
  suffix=""
  if [[ ${sample_counts[$sample]} -gt 1 ]]; then
    suffix="_${sample_counts[$sample]}"
  fi
  sample_label="${sample}${suffix}"

  # echo "🚀 Submitting jobs for $sample_label (original run: $run_name, condition: $condition)"

  # --- KMET job ---
  kmet_jobid=$(sbatch --account="$account" \
    --job-name="${sample_label}_kmet" \
    --output="${LOG_DIR}/kmet/${sample_label}.out" \
    --error="${LOG_DIR}/kmet/${sample_label}.err" \
    --ntasks=1 --cpus-per-task="$STEP1_CPUS" --mem-per-cpu=8gb \
    "$BIN_DIR/run_kmetstat_analysis.sh" "$prefix" "$sample_label" "$KMET_DIR" | awk '{print $4}')
  echo "✅ KMET job submitted: $kmet_jobid"
  kmet_job_ids+=("$kmet_jobid")

  # --- Trim job ---
  trim_jobid=$(sbatch --account="$account" \
    --job-name="${sample_label}_trim" \
    --output="${LOG_DIR}/trim/${sample_label}.out" \
    --error="${LOG_DIR}/trim/${sample_label}.err" \
    --ntasks=1 --cpus-per-task="$STEP1_CPUS" --mem-per-cpu=8gb \
    "$BIN_DIR/trim_reads.sh" "$prefix" "$sample_label" "$TRIM_DIR" | awk '{print $4}')
  echo "✅ Trim job submitted: $trim_jobid"
  trim_job_ids+=("$trim_jobid")

  sleep 1
done < <(tail -n +2 "$CSV")


# --- Wait for KMET jobs ---
echo "🕒 Waiting for all KMET jobs to finish..."
wait_for_jobs_array kmet_job_ids 1 30

# --- Wait for trim jobs ---
echo "🕒 Waiting for all trimming jobs to finish..."
wait_for_jobs_array trim_job_ids 1 30
sleep 120
# --- Submit MultiQC after all complete ---
multiqc_jobid=$(sbatch --account="$account" \
  --job-name=multiqc_step1 \
  --output="${LOG_DIR}/multiqc.out" \
  --error="${LOG_DIR}/multiqc.err" \
  --ntasks=1 --cpus-per-task=1 --mem-per-cpu=8gb \
  "$BIN_DIR/multi_QC.sh" "$KMET_DIR" "$TRIM_DIR" "$OUT_DIR" | awk '{print $4}')
echo "✅ MultiQC job submitted: $multiqc_jobid"

echo "🕒 Waiting for MultiQC job $multiqc_jobid to finish..."
wait_for_job_completion_less_output $multiqc_jobid 5

echo "✅ Step1 complete: KMET + Trimming + MultiQC finished for $NAME"

