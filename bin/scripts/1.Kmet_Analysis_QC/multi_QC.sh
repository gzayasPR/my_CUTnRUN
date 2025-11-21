#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --account=fengyue

# Load multiqc module
ml multiqc
kmet_dir=$1
trim_dir=$2
out_dir=$3 

# Run MultiQC across raw and trimmed FastQC folders
multiqc ${kmet_dir} ${trim_dir} \
  --outdir ${out_dir}/multiqc \
  --title "CUT&RUN KMET QC Report" \
  --force

echo "✅ MultiQC summary generated: ${out_dir}/multiqc/multiqc_report.html"
