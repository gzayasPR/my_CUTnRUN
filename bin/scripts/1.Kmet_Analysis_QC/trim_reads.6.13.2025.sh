#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --account=fengyue
# ============================================
# CUT&RUN Trimming + FastQC Script
# ============================================
# Usage: bash trim_reads_with_trimgalore.sh <fastq_prefix> <output_dir>
# Example: bash trim_reads_with_trimgalore.sh /path/to/sample/WT-H3_S1_L007 /path/to/results/trimmed_fastq

# Inputs
PREFIX=$1
RUN_NAME=$2
OUT_DIR=$3

if [[ -z "$PREFIX" || -z "$OUT_DIR" ]]; then
  echo "Usage: bash trim_reads_with_trimgalore.sh <fastq_prefix> <output_dir>"
  exit 1
fi

# Load necessary modules
ml trim_galore
ml fastqc

# Extract sample name from prefix
SAMPLE_NAME=$RUN_NAME

# Create output subdirectory for the sample
SAMPLE_DIR="${OUT_DIR}/${RUN_NAME}"
mkdir -p "$SAMPLE_DIR"

# Define input files
R1="${PREFIX}_R1_001.fastq.gz"
R2="${PREFIX}_R2_001.fastq.gz"

# Run Trim Galore (paired-end)
echo "🔧 Trimming $SAMPLE_NAME..."
trim_galore --length 20 -quality 20 --stringency 3 --trim-n \
            --paired "$R1" "$R2" --gzip -o "$SAMPLE_DIR" -j 4

echo "✅ Done trimming: $SAMPLE_NAME"
TRIMMED_R1="${SAMPLE_DIR}/${SAMPLE_NAME}_R1_001_val_1.fq.gz"
TRIMMED_R2="${SAMPLE_DIR}/${SAMPLE_NAME}_R2_001_val_2.fq.gz"

PREFIX_NAME=$(basename $R1 _R1_001.fastq.gz)
mv "${SAMPLE_DIR}/${PREFIX_NAME}_R1_001_val_1.fq.gz" ${TRIMMED_R1}
mv "${SAMPLE_DIR}/${PREFIX_NAME}_R2_001_val_2.fq.gz" ${TRIMMED_R2}
# Run FastQC on trimmed files

echo "🔬 Running FastQC on trimmed reads..."
fastqc "$TRIMMED_R1" "$TRIMMED_R2" -o "$SAMPLE_DIR"

echo "✅ FastQC complete for $SAMPLE_NAME"
