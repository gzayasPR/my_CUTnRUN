#!/bin/bash
#SBATCH --time=24:00:00


# --- Inputs ---
sample=$1              # sample name (e.g., WT-H3_S1_L007)
ALIGN_DIR=$2
NORM_DIR=$3           # directory with filtered BAMs and BigWigs
OUT_DIR=$4           # where to save QC plots
MOUSE_GTF=$5             # BED file of TSSs (e.g., mm10_TSS.bed)
CPUS=$6                # number of threads to use


# NORM_DIR="${NORM_DIR}"
# OUT_DIR="${QC_DIR}"
# CPUS=1
# --- Paths ---
sample_DIR="${NORM_DIR}/${sample}"
# BAM="${ALIGN_DIR}/${sample}//mouse_alignment/${sample}_mouse_filtered.bam"
BAM="${ALIGN_DIR}/${sample}/${sample}.mouse.bam"
BIGWIG="${sample_DIR}/${sample}_mouse_scaled.bw"

mkdir -p "${OUT_DIR}/${sample}"

# --- Load modules ---
module load samtools
module load deeptools

# --- Fragment Size Distribution ---
bamPEFragmentSize \
  -b "$BAM" \
  --hist "${OUT_DIR}/${sample}/fragment_size_distribution.png" \
  --maxFragmentLength 600 \
  --samplesLabel "$sample" \
  --plotTitle "Fragment Size Distribution - $sample" \
  --outRawFragmentLengths "${OUT_DIR}/${sample}/fragment_lengths.txt" \
  -p "$CPUS"

# --- Signal at TSSs ---
computeMatrix reference-point \
  -S "$BIGWIG" \
  -R "$MOUSE_GTF" \
  -a 2000 -b 2000 \
  -out "${OUT_DIR}/${sample}/matrix_TSS.gz" \
  --referencePoint TSS \
  -p "$CPUS"

plotHeatmap \
  -m "${OUT_DIR}/${sample}/matrix_TSS.gz" \
  -out "${OUT_DIR}/${sample}/TSS_heatmap.png" \
  --plotTitle "Signal around TSSs - $sample"

echo "[✓] QC completed for $sample"
