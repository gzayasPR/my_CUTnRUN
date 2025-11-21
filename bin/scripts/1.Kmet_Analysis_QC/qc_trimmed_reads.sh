#!/bin/bash
#SBATCH --account=fengyue
# ============================================
# Run FastQC and MultiQC on Trimmed Reads
# ============================================
# Usage: bash qc_trimmed_reads.sh <trimmed_reads_dir>
# Example: bash qc_trimmed_reads.sh /path/to/results/trimmed_fastq

TRIMMED_DIR=$1

if [[ -z "$TRIMMED_DIR" ]]; then
  echo "Usage: bash qc_trimmed_reads.sh <trimmed_reads_dir>"
  exit 1
fi
ml fastqc
ml multiqc
# Output directory for reports
FASTQC_OUT="${TRIMMED_DIR}/fastqc_reports"
MULTIQC_OUT="${TRIMMED_DIR}/multiqc_report"
mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

echo "🔍 Running FastQC on all trimmed reads in $TRIMMED_DIR..."
find "$TRIMMED_DIR" -type f \( -name "*_val_1.fq.gz" -o -name "*_val_2.fq.gz" \) | while read file; do
  echo "  - FastQC: $(basename "$file")"
  fastqc -o "$FASTQC_OUT" "$file"
done

echo "📊 Running MultiQC..."
multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT"

echo "✅ QC complete. MultiQC report saved to: $MULTIQC_OUT"
