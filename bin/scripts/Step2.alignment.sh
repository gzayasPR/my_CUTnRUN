#!/bin/bash

# Step 1: Read parameter file passed as argument
PARAM_FILE=$1
source "${PARAM_FILE}"  # Load all variables defined in the parameter file
COMB_INDEX="${my_data}/references/combined_ref/mouse_ecoli"
# Step 2: Validate required parameters
REQUIRED_VARS=(NAME CSV account STEP2_CPUS ECOLI_INDEX MOUSE_INDEX BLACKLIST)
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var:-}" ]; then
    echo "❌ Missing required parameter: $var"
    exit 1
  fi
done

# Step 3: Load project-specific environment and utility functions
#source "${project_env}"                          # Typically sets my_results, my_bash, etc.
source "${my_bin}/scripts/shell.functions.sh"    # Includes helper functions like wait_for_jobs_array

# Step 4: Define directories and paths
TRIM_DIR="${my_results}/${NAME}/1.KMET_QC/trimmed"   # Input: Trimmed FASTQs
ALIGN_DIR="${my_results}/${NAME}/2.Alignment"        # Output: Alignment results
LOG_DIR="${my_bash}/${NAME}/Step2"                   # SLURM output logs
BIN_DIR="${my_bin}/scripts/2.Alignment"              # Directory with the alignment script
ALIGN_SCRIPT="${BIN_DIR}/alignment.sh"              # The alignment script to run via sbatch
ALIGN_DIR="${my_results}/${NAME}/2.Single_Alignment" 
# Step 5: Create output directories if not exist
mkdir -p "$ALIGN_DIR" "$LOG_DIR"

# Step 6: Status update to terminal
echo "✅ Step 2: Merging + Aligning FASTQs for $NAME"
echo "📌 Parameters:"
echo "  Metadata CSV: $CSV"
echo "  Trimmed input dir: $TRIM_DIR"
echo "  Output dir: $ALIGN_DIR"
echo "  Log dir: $LOG_DIR"

# Step 7: Track job IDs for SLURM submission
job_ids=()

# Step 8: Loop through unique sample names in the first column of the CSV
while read -r sample; do
  echo "🧬 Merging + aligning for $sample"

  SAMPLE_DIR="${ALIGN_DIR}/${sample}"         # Directory to hold sample-specific alignment
  mkdir -p "$SAMPLE_DIR"

  merged_R1="${SAMPLE_DIR}/${sample}_R1.fastq.gz"
  merged_R2="${SAMPLE_DIR}/${sample}_R2.fastq.gz"

  # Step 9: Clear previous merged files if they exist
  > "$merged_R1"
  > "$merged_R2"

  # Log file to record what files were merged
  echo "  🔗 Merging the following trimmed FASTQs for $sample:" > "${SAMPLE_DIR}/merge_log.txt"

  # Step 10: Identify all trimmed directories for this sample (e.g., sample, sample_2, sample_3)
  MATCHED_DIRS=$(find "$TRIM_DIR" -maxdepth 1 -type d -regex ".*/${sample}\(_[0-9]+\)*" | sort)

  if [[ -z "$MATCHED_DIRS" ]]; then
    echo "⚠️ No matching trimmed directories found for sample: $sample"
    continue
  fi

  # Step 11: Merge all matching runs for R1 and R2 FASTQ files
  for run_path in $MATCHED_DIRS; do
    run_name=$(basename "$run_path")
    r1_trimmed="${run_path}/${run_name}_R1_001_val_1.fq.gz"
    r2_trimmed="${run_path}/${run_name}_R2_001_val_2.fq.gz"

    # Skip if trimmed files are missing
    if [[ ! -f "$r1_trimmed" || ! -f "$r2_trimmed" ]]; then
      echo "⚠️ Missing trimmed files for $run_name, skipping..." | tee -a "${SAMPLE_DIR}/merge_log.txt"
      continue
    fi

    # Log what is being merged
    echo "    R1: $r1_trimmed" >> "${SAMPLE_DIR}/merge_log.txt"
    echo "    R2: $r2_trimmed" >> "${SAMPLE_DIR}/merge_log.txt"

    # Append to merged FASTQs
    cat "$r1_trimmed" >> "$merged_R1"
    cat "$r2_trimmed" >> "$merged_R2"
  done

  # Step 12: Report merged output
  echo "  ✅ Merged FASTQs written to:"
  echo "     $merged_R1"
  echo "     $merged_R2"

  # Step 13: Submit alignment job to SLURM
  # jobid=$(sbatch --account="$account" \
  #   --job-name="align_${sample}" \
  #   --output="${LOG_DIR}/align_${sample}.out" \
  #   --error="${LOG_DIR}/align_${sample}.err" \
  #   --ntasks=1 \
  #   --cpus-per-task="$STEP2_CPUS" \
  #   --mem-per-cpu=8gb \
  #   "$ALIGN_SCRIPT" "$merged_R1" "$merged_R2" "$MOUSE_INDEX" "$ECOLI_INDEX" "$SAMPLE_DIR" "$BLACKLIST" \
  #   | awk '{print $4}')  # Extract job ID from sbatch output
  ALIGN_SCRIPT="${BIN_DIR}/single_alignment.sh"    
  jobid=$(sbatch --account="$account" \
    --job-name="align_${sample}" \
    --output="${LOG_DIR}/align_${sample}.out" \
    --error="${LOG_DIR}/align_${sample}.err" \
    --ntasks=1 \
    --cpus-per-task="$STEP2_CPUS" \
    --mem-per-cpu=8gb \
    "$ALIGN_SCRIPT" "$merged_R1" "$merged_R2" "$COMB_INDEX" "$SAMPLE_DIR" "$BLACKLIST" \
    | awk '{print $4}')  # Extract job ID from sbatch output

  echo "✅ Submitted alignment job for $sample: $jobid"
  job_ids+=("$jobid")
  sleep 1  # Avoid overloading SLURM scheduler
done < <(cut -d',' -f1 "$CSV" | tail -n +2 | sort | uniq)  # Extract unique sample names (skip header)

# Step 14: Wait for all SLURM jobs to complete using helper function
echo "📋 Waiting for submitted alignment jobs to finish..."
wait_for_jobs_array job_ids 1 30  # Check every 30 seconds, max 1 retry per job (custom function)
echo "✅ Step 2: All merging + alignment steps completed successfully."

