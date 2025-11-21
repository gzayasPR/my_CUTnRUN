#!/bin/bash

PARAM_FILE=$1

source "${PARAM_FILE}"
# ---- VALIDATE REQUIRED PARAMS ----
REQUIRED_VARS=(NAME CSV account STEP5_CPUS MOUSE_GTF MOUSE_BED)
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var:-}" ]; then
    echo "❌ Missing required parameter: $var"
    exit 1
  fi
done

source "$my_bin/scripts/shell.functions.sh"

NORM_DIR="$my_results/$NAME/3.Normalization"
PEAK_DIR="$my_results/$NAME/4.PeakCalling"
LOG_DIR="$my_bash/$NAME/Step5"
BIN_DIR="$my_bin/scripts/5.Annotation"
mkdir -p "$PEAK_DIR" "$LOG_DIR"

echo "Starting Step 5: Enrichment and Peak Comparison"

# --- Parse non-IgG samples and their conditions ---
declare -A condition_to_sample

while IFS=',' read -r sample other condition _; do
  [[ "$condition" =~ IgG ]] && continue
  condition_to_sample["$condition"]="$sample"  # Grab first sample per condition
done < <(tail -n +2 "$CSV")  # Skip header

# --- Get all non-IgG unique conditions ---
H3_conditions=("${!condition_to_sample[@]}")

# --- Run SEACR Peak Comparisons Automatically ---
echo "📊 Running automated SEACR peak comparisons..."

PEAK_COMPARISON_DIR="$my_results/$NAME/5.PeakComparison"
mkdir -p "$PEAK_COMPARISON_DIR"
seacr_job_ids=()

for (( i=0; i<${#H3_conditions[@]}; i++ )); do
  for (( j=i+1; j<${#H3_conditions[@]}; j++ )); do
    COND1="${H3_conditions[i]}"
    COND2="${H3_conditions[j]}"
    SAMPLE1="${condition_to_sample[$COND1]}"
    SAMPLE2="${condition_to_sample[$COND2]}"

    PEAK1="$PEAK_DIR/$SAMPLE1/SEACR/${SAMPLE1}_SEACR_peaks.stringent.bed"
    PEAK2="$PEAK_DIR/$SAMPLE2/SEACR/${SAMPLE2}_SEACR_peaks.stringent.bed"

    OUT_NAME="${COND1}_vs_${COND2}"
    OUT_DIR="$PEAK_COMPARISON_DIR/$OUT_NAME/PEAKS_SEACR"
    mkdir -p "$OUT_DIR"

    if [[ -f "$PEAK1" && -f "$PEAK2" ]]; then
      echo "🔬 Comparing $COND1 ($SAMPLE1) vs $COND2 ($SAMPLE2)..."
      jobid=$(sbatch --account="$account" \
        --ntasks=1 --cpus-per-task="$STEP5_CPUS" --mem-per-cpu=8gb \
        --job-name="Compare_SEACR" \
        --output="$LOG_DIR/Compare_SEACR_${COND1}_vs_${COND2}.out" \
        --error="$LOG_DIR/Compare_SEACR_${COND1}_vs_${COND2}.err" \
        "$BIN_DIR/compare_seacr_peaks.sh" "$PEAK1" "$PEAK2" "$OUT_DIR" | awk '{print $4}')
      seacr_job_ids+=("$jobid")
    else
      echo "⚠️ Missing peak files for $SAMPLE1 or $SAMPLE2 — skipping"
    fi
  done
done

echo "🕒 Waiting for SEACR job to finish..."
wait_for_jobs_array seacr_job_ids 1 30
# --- Enrichment Comparison Using Sample Names ---
enrich_job_ids=()

for (( i=0; i<${#H3_conditions[@]}; i++ )); do
  for (( j=i+1; j<${#H3_conditions[@]}; j++ )); do
    COND1="${H3_conditions[i]}"
    COND2="${H3_conditions[j]}"
    SAMPLE1="${condition_to_sample[$COND1]}"
    SAMPLE2="${condition_to_sample[$COND2]}"

    NORM1="$NORM_DIR/$SAMPLE1/${SAMPLE1}_mouse_scaled.bw"
    NORM2="$NORM_DIR/$SAMPLE2/${SAMPLE2}_mouse_scaled.bw"

    OUT_NAME="${COND1}_vs_${COND2}"
    OUT_DIR="$PEAK_COMPARISON_DIR/$OUT_NAME/Enrichment"
    SCEAR_DIR="$PEAK_COMPARISON_DIR/$OUT_NAME/PEAKS_SEACR/"
    mkdir -p "$OUT_DIR"

    SEACR_SHARED="$SCEAR_DIR/shared_peaks_union.bed"
    SEACR_SAMPLE1_ONLY="$SCEAR_DIR/${SAMPLE1}_SEACR_peaks.stringent_only.bed"
    SEACR_SAMPLE2_ONLY="$SCEAR_DIR/${SAMPLE2}_SEACR_peaks.stringent_only.bed"

    if [[ -f "$NORM1" && -f "$NORM2" ]]; then
      echo "🔬 Enrichment comparison for $COND1 ($SAMPLE1) vs $COND2 ($SAMPLE2)..."
      jobid=$(sbatch --account="$account" \
        --ntasks=1 --cpus-per-task="$STEP5_CPUS" --mem-per-cpu=8gb \
        --job-name="Enrichment_${COND1}_vs_${COND2}" \
        --output="$LOG_DIR/Enrichment_${COND1}_vs_${COND2}.out" \
        --error="$LOG_DIR/Enrichment_${COND1}_vs_${COND2}.err" \
        "$BIN_DIR/fold_enrichment.sh" "$NORM1" "$NORM2" "$OUT_DIR" "$proj_env" "$SEACR_SHARED" "$SEACR_SAMPLE1_ONLY" "$SEACR_SAMPLE2_ONLY" "$MOUSE_GTF" | awk '{print $4}')
      enrich_job_ids+=("$jobid")
    else
      echo "⚠️ Missing NORM files for $SAMPLE1 or $SAMPLE2 — skipping"
    fi
  done
done


echo "🕒 Waiting for Enrichment job to finish..."
wait_for_jobs_array enrich_job_ids 1 30

CSV_PEAKS="$OUT_DIR/combined_peaks_log2FC.csv"
echo "📊 Submitting TSS profiling jobs..."

# Build maps from CSV: condition → sample and corresponding BigWig
declare -A bw_map igg_bw
while IFS=',' read -r sample run_name condition fastq1 fastq2; do
  bw_path="$NORM_DIR/$sample/${sample}_mouse_scaled.bw"
  if [[ "$condition" =~ IgG ]]; then
    igg_bw["$condition"]="$bw_path"
  else
    bw_map["$condition"]="$bw_path"
    condition_to_sample["$condition"]="$sample"
  fi
done < <(tail -n +2 "$CSV")

# Use previously defined H3_conditions
for (( i=0; i<${#H3_conditions[@]}; i++ )); do
  for (( j=i+1; j<${#H3_conditions[@]}; j++ )); do
    COND1="${H3_conditions[i]}"
    COND2="${H3_conditions[j]}"

    SAMPLE1="${condition_to_sample[$COND1]}"
    SAMPLE2="${condition_to_sample[$COND2]}"

    WT_H3_BW="${bw_map[$COND1]}"
    KO_H3_BW="${bw_map[$COND2]}"
    WT_IgG_BW="${igg_bw[${COND1}_IgG]}"
    KO_IgG_BW="${igg_bw[${COND2}_IgG]}"

    OUT_TSS_DIR="$PEAK_COMPARISON_DIR/${COND1}_vs_${COND2}/TSS_Compare"
    mkdir -p "$OUT_TSS_DIR"

    if [[ ! -f "$WT_H3_BW" || ! -f "$KO_H3_BW" || ! -f "$WT_IgG_BW" || ! -f "$KO_IgG_BW" ]]; then
      echo "⚠️ Missing BW for $COND1/$COND2 – skipping TSS profiling"
      continue
    fi

    echo "📈 Submitting TSS job: $COND1 vs $COND2"
    bw_jobid=$(sbatch --account="$account" \
      --ntasks=1 --cpus-per-task="$STEP5_CPUS" --mem-per-cpu=8gb \
      --job-name="TSS_${COND1}_vs_${COND2}" \
      --output="$LOG_DIR/TSS_${COND1}_vs_${COND2}.out" \
      --error="$LOG_DIR/TSS_${COND1}_vs_${COND2}.err" \
      "$BIN_DIR/TSS.Start.sh" "$OUT_TSS_DIR" "$MOUSE_BED" "$CSV_PEAKS" "$WT_H3_BW" "$KO_H3_BW" "$WT_IgG_BW" "$KO_IgG_BW" | awk '{print $4}') 
  done
done


CSV_PEAKS="$OUT_DIR/combined_peaks_log2FC.csv"
echo "📊 Submitting TSS profiling jobs (BAM version)..."

# Build maps from CSV: condition → sample and corresponding BAM
declare -A bam_map igg_bam condition_to_sample
while IFS=',' read -r sample run_name condition fastq1 fastq2; do
  bam_path="$NORM_DIR/$sample/${sample}_mouse_filtered.bam"
  if [[ "$condition" =~ IgG ]]; then
    igg_bam["$condition"]="$bam_path"
  else
    bam_map["$condition"]="$bam_path"
    condition_to_sample["$condition"]="$sample"
  fi
done < <(tail -n +2 "$CSV")

# Iterate through H3 comparisons
for (( i=0; i<${#H3_conditions[@]}; i++ )); do
  for (( j=i+1; j<${#H3_conditions[@]}; j++ )); do
    COND1="${H3_conditions[i]}"
    COND2="${H3_conditions[j]}"

    SAMPLE1="${condition_to_sample[$COND1]}"
    SAMPLE2="${condition_to_sample[$COND2]}"

    WT_H3_BAM="${bam_map[$COND1]}"
    KO_H3_BAM="${bam_map[$COND2]}"
    WT_IgG_BAM="${igg_bam[${COND1}_IgG]}"
    KO_IgG_BAM="${igg_bam[${COND2}_IgG]}"

    OUT_TSS_DIR="$PEAK_COMPARISON_DIR/${COND1}_vs_${COND2}/TSS_Compare_BAM"
    mkdir -p "$OUT_TSS_DIR"

    if [[ ! -f "$WT_H3_BAM" || ! -f "$KO_H3_BAM" || ! -f "$WT_IgG_BAM" || ! -f "$KO_IgG_BAM" ]]; then
      echo "⚠️ Missing BAM for $COND1/$COND2 – skipping TSS profiling"
      continue
    fi

    echo "📈 Submitting TSS job (BAM): $COND1 vs $COND2"
    bam_jobid=$(sbatch --account="$account" \
      --ntasks=1 --cpus-per-task="$STEP5_CPUS" --mem-per-cpu=8gb \
      --job-name="TSS_${COND1}_vs_${COND2}" \
      --output="$LOG_DIR/BAM_TSS_${COND1}_vs_${COND2}.out" \
      --error="$LOG_DIR/BAM_TSS_${COND1}_vs_${COND2}.err" \
      "$BIN_DIR/TSS.Start.BAM.sh" "$OUT_TSS_DIR" "$MOUSE_BED" "$CSV_PEAKS" \
      "$WT_H3_BAM" "$KO_H3_BAM" "$WT_IgG_BAM" "$KO_IgG_BAM" | awk '{print $4}')


  done
done
wait_for_job_completion_less_output $bw_jobid 1 30
wait_for_job_completion_less_output $bam_jobid 1 30