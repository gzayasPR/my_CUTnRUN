#!/bin/bash

# Usage: bash run_fastqc_multiqc.sh <out_dir>
# Example: bash run_fastqc_multiqc.sh /path/to/my_results/1.kmet_analysis

set -e

OUT_DIR=$1

if [ -z "$OUT_DIR" ]; then
  echo "Usage: bash run_fastqc_multiqc.sh <out_dir>"
  exit 1
fi

# Load FastQC and MultiQC if on HPC or module system
# module load fastqc
# module load multiqc

echo "🔍 Running FastQC on all .fastq files in $OUT_DIR..."
ml fastqc
ml multiqc
# Step 1: Run FastQC in each sample directory
find "$OUT_DIR" -type f -name "*.fastq" | while read fastq_file; do
  sample_dir=$(dirname "$fastq_file")
  echo "  - FastQC: $(basename "$fastq_file")"
  fastqc -o "$sample_dir" "$fastq_file"
done

# Step 2: Run MultiQC on the whole directory
echo "📊 Running MultiQC summary..."
mkdir -p "$OUT_DIR/multiqc_report"
multiqc "$OUT_DIR" -o "$OUT_DIR/multiqc_report"

echo "✅ Done! MultiQC report saved to $OUT_DIR/multiqc_report/"
