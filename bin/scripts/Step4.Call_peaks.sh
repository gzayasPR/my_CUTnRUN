#!/bin/bash

PARAM_FILE=$1

source "${PARAM_FILE}"
# ---- VALIDATE REQUIRED PARAMS ----
REQUIRED_VARS=(NAME CSV account STEP4_CPUS MOUSE_GTF)
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var:-}" ]; then
    echo "❌ Missing required parameter: $var"
    exit 1
  fi
done

# source "${proj_env}"
source "${my_bin}/scripts/shell.functions.sh"

NORM_DIR="${my_results}/${NAME}/3.Normalization"
PEAK_DIR="${my_results}/${NAME}/4.PeakCalling"
LOG_DIR="${my_bash}/${NAME}/Step4"
BIN_DIR="${my_bin}/scripts/4.PeakCalling"

mkdir -p "$PEAK_DIR" "$LOG_DIR"
echo "Starting Step 4: Calling peaks with MACS2 and SEACR."

peak_job_ids=()
# --- Submit MACS2 peak calling jobs ---
echo "📡 Submitting MACS2 peak calling jobs..."
while IFS=',' read -r sample run_name condition fastq1 fastq2; do
  [[ "$condition" == *IgG ]] && continue  # Skip IgG controls

  BAM="${NORM_DIR}/${sample}/${sample}_mouse_filtered.bam"
  OUT_DIR="${PEAK_DIR}/${sample}/MACS2"
  mkdir -p "$OUT_DIR"

  jobid=$( sbatch --account="$account" \
                  --job-name="macs2_${sample}" \
                  --output="${LOG_DIR}/macs2_${sample}.out" \
                  --error="${LOG_DIR}/macs2_${sample}.err" \
                  --ntasks=1 --cpus-per-task="$STEP4_CPUS" --mem-per-cpu=8gb \
                  "${BIN_DIR}/macs2_call.sh" \
                  "$sample" "$BAM" "$OUT_DIR" "$STEP4_CPUS" \
                  | awk '{print $4}' )

  peak_job_ids+=("$jobid")
  echo "✅ Submitted MACS2 job for $sample: $jobid"
  sleep 1
done < <(tail -n +2 "$CSV")

# --- Build IgG control map ---
declare -A igg_controls
while IFS=',' read -r sample run_name condition fastq1 fastq2; do
  if [[ "$condition" == *IgG ]]; then
    base_cond="${condition/_IgG/}"
    BW="${NORM_DIR}/${sample}/${sample}_mouse_scaled.bw"
    igg_controls["$base_cond"]="$BW"
  fi
done < <( tail -n +2 "$CSV" )
# --- Submit SEACR peak calling jobs ---
echo "🌊 Submitting SEACR peak calling jobs..."
while IFS=',' read -r sample run_name condition fastq1 fastq2; do
  [[ "$condition" == *IgG ]] && continue  # Skip IgG samples

  target_bw="${NORM_DIR}/${sample}/${sample}_mouse_scaled.bw"
  control_bw="${igg_controls[$condition]}"

  if [[ ! -f "$target_bw" || ! -f "$control_bw" ]]; then
    echo "⚠️ Missing BW for $sample or control ${condition}_IgG – skipping"
    continue
  fi

  OUT_DIR="${PEAK_DIR}/${sample}/SEACR"
  mkdir -p "$OUT_DIR"

  jobid=$( sbatch --account="$account" \
                  --job-name="seacr_${sample}" \
                  --output="${LOG_DIR}/seacr_${sample}.out" \
                  --error="${LOG_DIR}/seacr_${sample}.err" \
                  --ntasks=1 --cpus-per-task="$STEP4_CPUS" --mem-per-cpu=8gb \
                  "${BIN_DIR}/seacr_callv2.sh" \
                  "$sample" "$target_bw" "$control_bw" "$OUT_DIR" "$STEP4_CPUS" \
                  "$MOUSE_GTF" "$proj_env" \
                  | awk '{print $4}' )

  peak_job_ids+=("$jobid")
  echo "✅ Submitted SEACR job for $sample (ctrl=${condition}_IgG): $jobid"
  sleep 1
done < <(tail -n +2 "$CSV")

# --- Debug: list all collected job IDs ---
echo
echo "All peak-calling job IDs:"
printf "  %s\n" "${peak_job_ids[@]}"

echo "🕒 Waiting for MACS2 and SEACR jobs to finish..."
wait_for_jobs_array peak_job_ids 1 30

echo "🎯 Step 4 complete: Peaks called."
sleep 60
