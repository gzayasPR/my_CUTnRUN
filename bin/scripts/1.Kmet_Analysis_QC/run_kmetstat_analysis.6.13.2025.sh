#!/bin/bash
#SBATCH --job-name=kmetstat
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=48:00:00

# ============================================
# CUT&RUN Kmetstat + FastQC Script
# ============================================
# Usage: sbatch run_kmetstat_analysis.sh <fastq_prefix> <output_dir>

set -e  # Exit on error

# Input arguments
prefix_path=$1
run_name=$2
out_dir=$3

# Check input
if [[ -z "$prefix_path" || -z "$out_dir" ]]; then
  echo "Usage: sbatch run_kmetstat_analysis.sh <fastq_prefix> <output_dir>"
  exit 1
fi

# Load modules
ml fastqc

# Extract sample name
sample_name=$run_name

# Create output directories
sample_out_dir="${out_dir}/${sample_name}"
mkdir -p "$sample_out_dir"
mkdir -p logs

# Define input files
R1="${prefix_path}_R1_001.fastq.gz"
R2="${prefix_path}_R2_001.fastq.gz"

# Unzip fastq files into working directory
echo "📂 Unzipping FASTQ files for $sample_name..."
gunzip -c "$R1" > "${sample_out_dir}/${sample_name}_R1.fastq"
gunzip -c "$R2" > "${sample_out_dir}/${sample_name}_R2.fastq"

R1_FASTQ="${sample_out_dir}/${sample_name}_R1.fastq"
R2_FASTQ="${sample_out_dir}/${sample_name}_R2.fastq"

# 🔬 Run FastQC
echo "🔬 Running FastQC on raw reads for $sample_name..."
fastqc "$R1_FASTQ" "$R2_FASTQ" -o "$sample_out_dir"

# Barcode list
BARCODES=(
  TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG
  CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG
  CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA
  TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG
  TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA
  CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA
  ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG
  CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG
  CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA
  ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA
  CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC
  GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG
  CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG
  GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA
  CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG
  CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA
)

# Output file for barcode counts
RESULT_FILE="${sample_out_dir}/${sample_name}_barcode_counts.csv"
echo "Barcode,R1_Count,R2_Count" > "$RESULT_FILE"

# Count barcodes
echo "🧬 Counting barcodes for $sample_name..."
for barcode in "${BARCODES[@]}"; do
  r1_count=$(grep -c "$barcode" "$R1_FASTQ")
  r2_count=$(grep -c "$barcode" "$R2_FASTQ")
  echo "$barcode,$r1_count,$r2_count" >> "$RESULT_FILE"
done

echo "✅ Done: ${sample_name} → $RESULT_FILE"
