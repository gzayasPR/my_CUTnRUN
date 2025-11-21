#!/bin/bash
set -euo pipefail

# --- INPUTS ---
KO_BIGWIG=$1
WT_BIGWIG=$2
OUTDIR=$3
proj_env=$4
SEACR_SHARED=$5
SEACR_KO_ONLY=$6
SEACR_WT_ONLY=$7

# KO_BIGWIG="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/3.Normalization/KO_H3/KO_H3_mouse_scaled.bw"
# WT_BIGWIG="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/3.Normalization/WT_H3/WT_H3_mouse_scaled.bw"
# OUTDIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/Enrichment"
# proj_env="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/bin/project_env.sh"
# SEACR_SHARED="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/PEAKS_SEACR/shared_peaks_union.bed"
# SEACR_KO_ONLY="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/PEAKS_SEACR//KO_H3_SEACR_peaks.stringent_only.bed"
# SEACR_WT_ONLY="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/PEAKS_SEACR//WT_H3_SEACR_peaks.stringent_only.bed"



# --- ENVIRONMENT ---
source "$proj_env"
mkdir -p "$OUTDIR"

WT_NAME=$(basename "$WT_BIGWIG" .bw)
KO_NAME=$(basename "$KO_BIGWIG" .bw)

# --- FILE PATHS ---
LOG2BW="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC.bw"s
LOG2BEDG="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC.bedgraph"
LOG2BEDG_SORTED="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC_sorted.bedgraph"

KO_BEDG_SORTED="${OUTDIR}/${KO_NAME}_sorted.bedgraph"
WT_BEDG_SORTED="${OUTDIR}/${WT_NAME}_sorted.bedgraph"

# --- LOGGING ---
LOGFILE="${OUTDIR}/fold_enrichment_run.log"
echo "📄 Logging run to $LOGFILE"
{
  echo "Date: $(date)"
  echo "WT BigWig: $WT_BIGWIG"
  echo "KO BigWig: $KO_BIGWIG"
  echo "SEACR Shared: $SEACR_SHARED"
  echo "SEACR WT-only: $SEACR_WT_ONLY"
  echo "SEACR KO-only: $SEACR_KO_ONLY"
  echo "Output Dir: $OUTDIR"
} > "$LOGFILE"

# --- MODULES ---
ml deeptools
ml bedtools
ml conda

# --- STEP 1: Compute log2FC ---

# --- STEP 2: Convert BigWigs to BedGraphs ---
echo "📁 Converting BigWigs to BedGraphs..."
conda activate "${my_softwares}/conda_envs/bigwig"
bigWigToBedGraph "$KO_BIGWIG" "${OUTDIR}/${KO_NAME}.bedgraph"
bigWigToBedGraph "$WT_BIGWIG" "${OUTDIR}/${WT_NAME}.bedgraph"
conda deactivate

sort -k1,1 -k2,2n "${OUTDIR}/${KO_NAME}.bedgraph" > "$KO_BEDG_SORTED"
sort -k1,1 -k2,2n "${OUTDIR}/${WT_NAME}.bedgraph" > "$WT_BEDG_SORTED"

# --- STEP 3: Sort all peak files ---
for file in "$SEACR_SHARED" "$SEACR_WT_ONLY" "$SEACR_KO_ONLY"; do
  base=$(basename "$file" .bed)
  sort -k1,1 -k2,2n "$file" > "${OUTDIR}/${base}.sorted.bed"
done

# --- STEP 4: Function to map signal stats ---
map_stats () {
  label=$1           # shared, wt_only, ko_only
  peak_bed=$2        # sorted peak file
  echo "📌 Mapping WT/KO stats for: $label"
 echo "📌 Mapping collapsed signal for: $label"
  bedtools map -a "$peak_bed" -b "$WT_BEDG_SORTED" -c 4 -o collapse > "${OUTDIR}/${label}_WT_collapse.bed"
  bedtools map -a "$peak_bed" -b "$KO_BEDG_SORTED" -c 4 -o collapse > "${OUTDIR}/${label}_KO_collapse.bed"

  summary_file="${OUTDIR}/${label}_log2FC_full_summary.bed"
  rm "$summary_file"
  # Define metrics
  metrics=(collapse)

  # Start building the paste command from the peak bed
  paste_cmd="paste \"$peak_bed\""

  # Add WT metric columns
  for m in "${metrics[@]}"; do
    paste_cmd+=" <(awk '{print \$NF}' \"${OUTDIR}/${label}_WT_${m}.bed\")"
  done

  # Add KO metric columns
  for m in "${metrics[@]}"; do
    paste_cmd+=" <(awk '{print \$NF}' \"${OUTDIR}/${label}_KO_${m}.bed\")"
  done
  # Evaluate and write output
  eval "$paste_cmd" >> "$summary_file"
  head -n1 "$summary_file"

}

# --- STEP 5: Run for all categories ---
map_stats "shared_peaks" "${OUTDIR}/$(basename "$SEACR_SHARED" .bed).sorted.bed"
map_stats "wt_only_peaks" "${OUTDIR}/$(basename "$SEACR_WT_ONLY" .bed).sorted.bed"
map_stats "ko_only_peaks" "${OUTDIR}/$(basename "$SEACR_KO_ONLY" .bed).sorted.bed"

echo "✅ All peak summaries generated in: $OUTDIR"

ml R
Rscript ${my_bin}/scripts/5.Annotation/peak_log2fold.R "${OUTDIR}"