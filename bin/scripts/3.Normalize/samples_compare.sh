#!/bin/bash
#SBATCH --time=48:00:00


# --- Arguments ---
project_env=$1
CSV=$2
ALIGN_DIR=$3
NORM_DIR=$4
cpus=$5

# --- Load environment ---
source "${project_env}"

# --- Load modules ---
module load samtools
module load deeptools

# --- Create output directory ---
OUT_DIR="${NORM_DIR}/all_samples"
mkdir -p "$OUT_DIR"

# --- Build list of BAM files and labels ---
BAM_FILES=()
SAMPLE_LABELS=()

while read -r sample; do
  # BAM_PATH="${ALIGN_DIR}/${sample}/mouse_alignment/${sample}_mouse_filtered.bam"
   BAM_PATH="${ALIGN_DIR}/${sample}/${sample}.mouse.bam"
  if [[ -f "$BAM_PATH" ]]; then
    BAM_FILES+=("$BAM_PATH")
    SAMPLE_LABELS+=("$sample")
  else
    echo "⚠️ Warning: BAM file not found for $sample at $BAM_PATH"
  fi
done < <(tail -n +2 "$CSV" | cut -d',' -f1 | sort | uniq)

# --- Convert arrays to strings ---
BAM_FILES_STR="${BAM_FILES[*]}"
SAMPLE_LABELS_STR="${SAMPLE_LABELS[*]}"

echo $BAM_FILES_STR

# # --- Fragment Size Distribution ---
bamPEFragmentSize \
  -b $BAM_FILES_STR \
  --hist "${OUT_DIR}/fragment_size_distribution.png" \
  --maxFragmentLength 600 \
  --samplesLabel $SAMPLE_LABELS_STR \
  --plotTitle "Fragment Size Distribution - All Samples" \
  --outRawFragmentLengths "${OUT_DIR}/fragment_lengths.txt" \
  -p "$cpus"


BIGWIG_FILES=()
SAMPLE_LABELS=()

while read -r sample; do
  BW_PATH="${NORM_DIR}/${sample}/${sample}_mouse_scaled.bw"
  if [[ -f "$BW_PATH" ]]; then
    BIGWIG_FILES+=("$BW_PATH")
    SAMPLE_LABELS+=("$sample")
  else
    echo "⚠️ Warning: BigWig file not found for $sample at $BW_PATH"
  fi
done < <(tail -n +2 "$CSV" | cut -d',' -f1 | sort | uniq)

BIGWIG_FILES_STR="${BIGWIG_FILES[*]}"
SAMPLE_LABELS_STR="${SAMPLE_LABELS[*]}"

echo $BIGWIG_FILES_STR
# --- Compute Coverage Matrix ---
multiBigwigSummary bins \
  -b $BIGWIG_FILES_STR \
  --outFileName "${OUT_DIR}/matrix.npz" \
  --outRawCounts "${OUT_DIR}/raw_counts.tab" \
  -p "$cpus"

# --- Plot Correlation ---
plotCorrelation \
  --corData "${OUT_DIR}/matrix.npz" \
  --corMethod pearson \
  --whatToPlot heatmap \
  --plotNumbers \
  --removeOutliers \
  --plotFile  "${OUT_DIR}/sample_correlation_heatmap.png" \
  --outFileCorMatrix "${OUT_DIR}/sample_correlation_matrix.tab" 

# --- Plot PCA ---
plotPCA \
  --corData "${OUT_DIR}/matrix.npz" \
  --transpose \
  --plotFile "${OUT_DIR}/PCA_plot.png" \
  --outFileNameData "${OUT_DIR}/PCA_data.tab" \
  --plotTitle "PCA of CUT&RUN Samples" 