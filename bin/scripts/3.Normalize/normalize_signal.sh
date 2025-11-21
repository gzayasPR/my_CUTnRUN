#!/bin/bash
#SBATCH --time=48:00:00

# --- Inputs ---
SAMPLE=$1
INPUT_DIR=$2         # where filtered BAMs and counts are
OUT_DIR=$3           # where to write scaled bigWigs
out_scale=$4    # full path to max_ecoli_read_count.txt
cpus=$5
# --- Paths ---
SAMPLE_DIR="${INPUT_DIR}/${SAMPLE}"
# MOUSE_BAM="${SAMPLE_DIR}/mouse_alignment/${SAMPLE}_mouse_filtered.bam"
MOUSE_BAM="${SAMPLE_DIR}/${SAMPLE}.mouse.bam"
BIGWIG_OUT="${OUT_DIR}/${SAMPLE}/${SAMPLE}_mouse_scaled.bw"
mkdir -p "${OUT_DIR}/${SAMPLE}/"
# ------------------------------------------------
echo "📝 Logging parameters..."
LOGFILE="${OUT_DIR}/${SAMPLE}/Norm.log"
{
  echo "📄 TSS Profiling Log - $(date)"
  echo "SAMPLE=$SAMPLE"
  echo "INPUT_DIR=$INPUT_DIR"
  echo "OUT_DIR=$OUT_DIR"
  echo "out_scale=$out_scale"
  echo "SAMPLE_DIR=$SAMPLE_DIR"
  echo "MOUSE_BAM=$MOUSE_BAM"
  echo "BIGWIG_OUT=$BIGWIG_OUT"
} > "$LOGFILE"

# --- Load modules ---
module load samtools
module load deeptools
SCALE_FACTOR=$(awk -v sample="$SAMPLE" 'BEGIN{FS="\t"} $1 == sample {print $7}' "$out_scale")

echo "[*] $SAMPLE scale factor: $SCALE_FACTOR"
samtools index "$MOUSE_BAM"
# # --- Run bamCoverage ---
# bamCoverage \
#   -b "$MOUSE_BAM" \
#   -o "$BIGWIG_OUT" \
#   --normalizeUsing None \
#   --scaleFactor "$SCALE_FACTOR" \
#   --binSize 10 \
#   --centerReads \
#   --smoothLength 50 \
#   --minMappingQuality 30 \
#   --ignoreForNormalization "chrM" \
#   --skipNonCoveredRegions \
#   -p max


bamCoverage -b "$MOUSE_BAM" -o "$BIGWIG_OUT" \
  --normalizeUsing RPGC --effectiveGenomeSize 2150570000 \
  --binSize 10 --centerReads --smoothLength 50 \
  --minMappingQuality 30 \
  --ignoreForNormalization "chrM" \
  --skipNonCoveredRegions   -p max


# # --- Run bamCoverage ---
# bamCoverage \
#     -b "$MOUSE_BAM" \
#     -o "$BIGWIG_OUT" \
#     --scaleFactor "$SCALE_FACTOR" \
#   --normalizeUsing None \
#   --binSize 10 --centerReads --ignoreDuplicates
  
echo "[✓] Normalized bigWig created for $SAMPLE"
