#!/bin/bash
#SBATCH --account=fengyue
#SBATCH --time=48:00:00
# ==========================================
# Align to E. coli, then align host reads to Mouse Genome
# ==========================================
# Usage: bash align_ecoli_then_mouse.sh <R1.fastq.gz> <R2.fastq.gz> <ecoli_index> <mouse_index> <output_dir>

#!/bin/bash

# --- INPUT ---
SAMPLE=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
cpus=$4
# --- PATHS ---
ECOLI_BAM="${INPUT_DIR}/${SAMPLE}/ecoli_alignment/${SAMPLE}_ecoli_sorted.bam"
MOUSE_BAM="${INPUT_DIR}/${SAMPLE}/mouse_alignment/${SAMPLE}_mouse_sorted.bam"
TOTAL_FASTQ_R1="${INPUT_DIR}/${SAMPLE}/${SAMPLE}_R1.fastq.gz"
TOTAL_FASTQ_R2="${INPUT_DIR}/${SAMPLE}/${SAMPLE}_R2.fastq.gz"
SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
mkdir -p "${SAMPLE_OUT}"

# --- MODULES (adapt for your HPC setup) ---
module load samtools
module load deeptools
if command -v zcat >/dev/null 2>&1; then ZCAT="zcat"
elif command -v gzcat >/dev/null 2>&1; then ZCAT="gzcat"
else
  echo "Need zcat or gzcat in PATH" >&2; exit 1
fi

# --- STEP 0: Index original BAMs ---
echo "[*] Indexing original BAMs..."
samtools index "$ECOLI_BAM"
samtools index "$MOUSE_BAM"

# --- STEP 2: Count spike-in reads ---
echo "[*] Counting E. coli reads for spike-in normalization..."

ECOLI_COUNT=$(samtools view -c "${ECOLI_BAM}")
echo "$ECOLI_COUNT" > "${SAMPLE_OUT}/${SAMPLE}_ecoli_read_count.txt"
echo "[+] E. coli read count for ${SAMPLE}: $ECOLI_COUNT"
MOUSE_COUNT=$(samtools view -c "${MOUSE_BAM}")
echo "$MOUSE_COUNT" > "${SAMPLE_OUT}/${SAMPLE}_mouse_read_count.txt"
echo "[+] E. coli read count for ${SAMPLE}: $MOUSE_COUNT"
# --- STEP 3: Count total reads from FASTQs ---
echo "[*] Counting total reads from FASTQs..."
R1_COUNT=$(( $(zcat "$TOTAL_FASTQ_R1" | wc -l) / 4 ))
R2_COUNT=$(( $(zcat "$TOTAL_FASTQ_R2" | wc -l) / 4 ))
# take the minimum to get number of pairs
TOTAL_COUNT=$(( R1_COUNT < R2_COUNT ? R1_COUNT : R2_COUNT ))

echo "$TOTAL_COUNT" > "${SAMPLE_OUT}/${SAMPLE}_total_read_count.txt"
echo "[+] Total read pairs for ${SAMPLE}: $TOTAL_COUNT"