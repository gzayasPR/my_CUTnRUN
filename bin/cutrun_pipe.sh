#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10G

# ---- INPUT ARGUMENTS ----
INPARAMS_FILE=$1
if [ ! -f "$INPARAMS_FILE" ]; then
  echo "❌ Parameter file $INPARAMS_FILE not found."
  exit 1
fi

# ---- LOAD PARAMS ----
source "$INPARAMS_FILE"

# ensure proj_env is set before sourcing
if [ -z "${proj_env:-}" ]; then
  echo "❌ Missing required parameter: proj_env"
  exit 1
fi
source "$proj_env"

# ---- VALIDATE REQUIRED PARAMS ----
REQUIRED_VARS=(NAME CSV account ECOLI_INDEX MOUSE_INDEX MOUSE_GTF my_bin my_bash)
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var:-}" ]; then
    echo "❌ Missing required parameter: $var"
    exit 1
  fi
done

echo "Loading environment and helper functions..."
source "${proj_env}"
source "${my_bin}/scripts/shell.functions.sh"

bash "${my_bin}/scripts/create_project_parameters.sh" "$1"

PARAM_DIR="${my_results}/${NAME}/parameters"
mkdir -p "$PARAM_DIR"
PARAM_FILE="${my_results}/${NAME}/parameters/project_env.parameters.sh"

source "${PARAM_FILE}"
echo "Project environment loaded from ${PARAM_FILE}"

# ---- SETUP ----
mkdir -p "${my_bash}/$NAME"

echo "🚀 Starting pipeline for project: $NAME"

# Helper function to check flags (case-insensitive for 1/true/yes)
run_step() {
  local var_name="$1"
  local flag_value="${!var_name,,}"
  [[ "$flag_value" =~ ^(1|true|yes)$ ]]
}

# ---- STEP 1: QC KMET ----
if run_step RUN_STEP1; then
  echo "🔍 Running Step 1: QC KMET"
  STEP1_CPUS=${STEP1_CPUS:-2}
  mkdir -p "${my_bash}/$NAME/Step1"
  bash "${my_bin}/scripts/Step1.QC_KMET.sh" "${PARAM_FILE}"
else
  echo "⏭️ Skipping Step 1: QC KMET"
fi

# ---- STEP 2: Alignment ----
if run_step RUN_STEP2; then
  echo "📡 Running Step 2: Alignment"
  STEP2_CPUS=${STEP2_CPUS:-6}
  mkdir -p "${my_bash}/$NAME/Step2"
  bash "${my_bin}/scripts/Step2.alignment.sh"  "${PARAM_FILE}"
else
  echo "⏭️ Skipping Step 2: Alignment"
fi 

# ---- POST-STEP 2: Alignment Summary Table ----
echo "📄 Running alignment summary..."

ALIGN_SUMMARY_SCRIPT="${my_bin}/scripts/2.Alignment/summarize_alignment_metrics.v2.sh"
QC_FILE="${my_results}/${NAME}/1.KMET_QC/multiqc/CUTRUN-KMET-QC-Report_multiqc_report_data/fastqc_sequence_counts_plot.txt"
ALIGN_DIR="${my_results}/${NAME}/2.Single_Alignment"
SUMMARY_OUT="${ALIGN_DIR}/alignment_summary.tsv"

bash "$ALIGN_SUMMARY_SCRIPT" "$NAME" "$QC_FILE" "$ALIGN_DIR" "$SUMMARY_OUT"



# ---- STEP 3: Normalization ----
if run_step RUN_STEP3; then
  echo "🎛️ Running Step 3: Normalization"
  STEP3_CPUS=${STEP3_CPUS:-8}
  mkdir -p "${my_bash}/$NAME/Step3"
  bash "${my_bin}/scripts/Step3.normalization.sh"  "${PARAM_FILE}"
else
  echo "⏭️ Skipping Step 3: Normalization"
fi

# ---- STEP 4: Peak Calling ----
if run_step RUN_STEP4; then
  echo "⛰️ Running Step 4: Peak Calling"
  STEP4_CPUS=${STEP4_CPUS:-8}
  mkdir -p "${my_bash}/$NAME/Step4"
  bash "${my_bin}/scripts/Step4.Call_peaks.sh"  "${PARAM_FILE}"
else
  echo "⏭️ Skipping Step 4: Peak Calling"
fi

# ---- STEP 5: Annotation ----
if run_step RUN_STEP5; then
  echo "📊 Running Step 5: Annotation"
  STEP5_CPUS=${STEP5_CPUS:-8}
  mkdir -p "${my_bash}/$NAME/Step5"
  bash "${my_bin}/scripts/Step5.Annotation.sh"  "${PARAM_FILE}"
else
  echo "⏭️ Skipping Step 5: Annotation"
fi

echo "✅ Pipeline complete."

echo "📦 Preparing final export and report..."

# Move into results directory
cd "${my_results}/${NAME}/" || exit 1

# Define export paths
EXPORT_DIR="Export_data"
rm -rf $EXPORT_DIR
BIGWIG_DIR="${EXPORT_DIR}/BigWig_Files"
PEAKS_DIR_SRC="5.PeakComparison/WT_vs_KO/Enrichment"
PEAKS_DIR_DEST="${EXPORT_DIR}/Peaks"
SOURCE_DIR="3.Normalization"
TSS_DIR_SRC="5.PeakComparison/WT_vs_KO/TSS_Compare"
TSS_DIR_DEST="${EXPORT_DIR}/TSS_Compare"
MULTIQC_DEST="${EXPORT_DIR}/MultiQC_Report"

# Create export structure
mkdir -p "$BIGWIG_DIR" "$PEAKS_DIR_DEST" "$TSS_DIR_DEST" "$MULTIQC_DEST"

# Export BigWigs
echo "📁 Exporting BigWig files..."
find "$SOURCE_DIR" -type f -name "*.bw" -exec cp -v {} "$BIGWIG_DIR/" \;

# Export peaks
echo "📁 Exporting Peak files..."
cp -v "${PEAKS_DIR_SRC}/combined_peaks_log2FC.csv" "$PEAKS_DIR_DEST/" 2>/dev/null
cp -v "${PEAKS_DIR_SRC}/combined_peaks_log2FC_annotated_merged.csv" "$PEAKS_DIR_DEST/" 2>/dev/null
cp -v "${PEAKS_DIR_SRC}"/*.pdf "$PEAKS_DIR_DEST/" 2>/dev/null

# Export TSS plots
echo "📁 Exporting TSS comparison plots..."
cp -v "${TSS_DIR_SRC}"/*.png "$TSS_DIR_DEST/" 2>/dev/null
cp -v "${TSS_DIR_SRC}"/*.pdf "$TSS_DIR_DEST/" 2>/dev/null

# Copy bash logs/scripts for reproducibility
echo "📂 Copying bash scripts into results..."
cp -r "${my_bash}/$NAME/" "bash_results/"

# Run MultiQC
echo "📊 Running MultiQC report..."
ml multiqc
multiqc . -o "$MULTIQC_DEST"

# Optional: Clean up verbose bash logs from results root
echo "🧹 Cleaning up intermediate logs..."
find . -type f -name "bashout.txt" -delete

rm -rf "bash_results/"
echo "✅ All exports and reports completed for project: $NAME"
