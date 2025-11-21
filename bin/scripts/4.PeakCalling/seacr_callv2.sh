#!/bin/bash
#SBATCH --time=12:00:00
# seacr_call.sh

# set -euo pipefail

sample=$1
target_bw=$2
control_bw=$3
out_dir=$4
cpus=$5
MOUSE_GTF=$6
proj_env=$7

source "${proj_env}"

# --- LOGGING ---
LOGFILE="${out_dir}/scear_${sample}_run.log"
echo "📄 Logging run to $LOGFILE"
{
  echo "Date: $(date)"
  echo "Sample: ${sample}"
  echo "Sample BigWig: ${target_bw}"
  echo "Control BigWig: ${control_bw}"
  echo "Output Dir: $out_dir"
} > "$LOGFILE"

cat "${LOGFILE}"
# Convert BigWig to BedGraph
target_bg="${out_dir}/${sample}_target.bedgraph"
control_bg="${out_dir}/${sample}_control.bedgraph"
cd "$out_dir"
ml conda
conda activate "${my_softwares}/conda_envs/bigwig"

bigWigToBedGraph "$target_bw" "$target_bg"
bigWigToBedGraph "$control_bw" "$control_bg"

ml bedtools
ml R
cd "$out_dir"
bash "${my_softwares}/SEACR/SEACR_1.3.sh" "$target_bg" "$control_bg" norm stringent "${out_dir}/${sample}_SEACR_peaks"

Rscript ${my_bin}/scripts/4.PeakCalling/annotate_peaks.R  "${out_dir}/${sample}_SEACR_peaks.stringent.bed" \
  "${out_dir}/${sample}_stringent" \
  "${MOUSE_GTF}"


# --- Convert to IGV-compatible BED6 (score = enrichment × 10) ---
awk 'BEGIN{OFS="\t"} {score = int($5*10); if (score > 1000) score = 1000; print "chr"$1, $2, $3, $6, score, "."}' \
  "${out_dir}/${sample}_SEACR_peaks.stringent.bed" > "${out_dir}/${sample}_SEACR_igv.bed"

# --- Add track line to top (optional) ---
sed -i "1itrack name=\"${sample}_SEACR\" description=\"SEACR peaks for ${sample}\" useScore=1" \
  "${out_dir}/${sample}_SEACR_igv.bed"

echo "✅ Created IGV BED for $sample → ${sample}_SEACR_igv.bed"