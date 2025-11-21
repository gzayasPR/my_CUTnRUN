#!/bin/bash
#SBATCH --time=24:00:00
# --- USAGE ---
# ./fold_enrichment_bigwig.sh path/to/WT.bw path/to/KO.bw path/to/output_dir

set -euo pipefail

# --- INPUTS ---
KO_BIGWIG=$1
WT_BIGWIG=$2
OUTDIR=$3
proj_env=$4
SEACR_SHARED=$5      # e.g. shared_peaks.bed
SEACR_KO_ONLY=$6     # e.g. KO_only.bed
SEACR_WT_ONLY=$7     # e.g. WT_only.bed


# KO_BIGWIG="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/3.Normalization/KO_H3/KO_H3_mouse_scaled.bw"
# WT_BIGWIG="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/3.Normalization/WT_H3/WT_H3_mouse_scaled.bw"
# OUTDIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/Enrichment"
# proj_env="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/bin/project_env.sh"
# SEACR_SHARED="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/PEAKS_SEACR//shared_peaks.bed"
# SEACR_KO_ONLY="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/PEAKS_SEACR//KO_H3_SEACR_peaks.stringent_only.bed"
# SEACR_WT_ONLY="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/PEAKS_SEACR//WT_H3_SEACR_peaks.stringent_only.bed"




source "$proj_env"
# --- SETUP ---
# --- SETUP ---
mkdir -p "$OUTDIR"
WT_NAME=$(basename "$WT_BIGWIG" .bw)
KO_NAME=$(basename "$KO_BIGWIG" .bw)

LOG2BW="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC.bw"
LOG2BEDG="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC.bedgraph"

LOG2_SHARED="${OUTDIR}/shared_peaks_log2FC.bed"
LOG2_WT_ONLY="${OUTDIR}/wt_only_peaks_log2FC.bed"
LOG2_KO_ONLY="${OUTDIR}/ko_only_peaks_log2FC.bed"

SUMMARY_TAB="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_summary.tab"



LOG2BW="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC.bw"
SUMMARY_TAB="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_summary.tab"
LOG2BEDG="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC.bedgraph"

# --- Logging parameters ---
LOGFILE="${OUTDIR}/fold_enrichment_run.log"
echo "📄 Logging run parameters to $LOGFILE"
{
  echo "Date: $(date)"
  echo "WT BigWig: $WT_BIGWIG"
  echo "KO BigWig: $KO_BIGWIG"
  echo "Output Directory: $OUTDIR"
  echo "Project Env: $proj_env"
  echo "SEACR Shared Peaks: $SEACR_SHARED"
  echo "SEACR WT-only Peaks: $SEACR_WT_ONLY"
  echo "SEACR KO-only Peaks: $SEACR_KO_ONLY"
  echo "Log2FC BigWig Output: $LOG2BW"
  echo "Log2FC BedGraph Output: $LOG2BEDG"
  echo ""
} > "$LOGFILE"

# # --- MODULES (adjust as needed for your HPC) ---
ml deeptools

# --- STEP 1: Compute genome-wide log2 fold change ---
echo "🧮 Calculating genome-wide log2 fold change..."
bigwigCompare \
  -b1 "$KO_BIGWIG" \
  -b2 "$WT_BIGWIG" \
  --operation log2 \
  --skipZeroOverZero \
  --pseudocount 1 \
  --binSize 50 \
  --numberOfProcessors 4 \
  -o "$LOG2BW"

ml conda
# --- STEP 2: Convert to BedGraph ---
echo "📁 Converting BigWig to BedGraph..."
conda activate "${my_softwares}/conda_envs/bigwig"
bigWigToBedGraph "$LOG2BW" "$LOG2BEDG"
LOG2BEDG_sorted="${OUTDIR}/${KO_NAME}_vs_${WT_NAME}_log2FC_sorted.bedgraph"
# --- STEP 3: Map log2FC to SEACR peak categories ---
echo "📌 Mapping log2FC to SEACR peak categories..."
ml bedtools
sort -k1,1 -k2,2n "$SEACR_SHARED" > "${OUTDIR}/shared_peaks.sorted.bed"
sort -k1,1 -k2,2n "$SEACR_WT_ONLY" > "${OUTDIR}/wt_only_peaks.sorted.bed"
sort -k1,1 -k2,2n "$SEACR_KO_ONLY" > "${OUTDIR}/ko_only_peaks.sorted.bed"
sort -k1,1 -k2,2n  "$LOG2BEDG" > "$LOG2BEDG_sorted"

bedtools map -a "${OUTDIR}/shared_peaks.sorted.bed"   -b "$LOG2BEDG_sorted" -c 4 -o mean > "$LOG2_SHARED"
bedtools map -a "${OUTDIR}/wt_only_peaks.sorted.bed"  -b "$LOG2BEDG_sorted" -c 4 -o mean > "$LOG2_WT_ONLY"
bedtools map -a "${OUTDIR}/ko_only_peaks.sorted.bed"  -b "$LOG2BEDG_sorted" -c 4 -o mean > "$LOG2_KO_ONLY"

# --- STEP 4: Plot using R ---
echo "📊 Plotting log2FC distributions..."
ml R

Rscript "${my_bin}/scripts/5.Annotation/log2fold_plotting.R" \
  "$LOG2BEDG" \
  "$OUTDIR" \
  "$LOG2_SHARED" \
  "$LOG2_WT_ONLY" \
  "$LOG2_KO_ONLY"

# # --- STEP 5: Optional summary stats over bins ---
# echo "📈 Generating bin summary matrix..."
# multiBigwigSummary bins \
#   --bwfiles "$KO_BIGWIG" "$WT_BIGWIG" \
#   --binSize 50 \
#   -out "${OUTDIR}/signal_summary.npz" \
#   --outRawCounts "$SUMMARY_TAB"

# echo "✅ All steps complete. Results saved to: $OUTDIR"