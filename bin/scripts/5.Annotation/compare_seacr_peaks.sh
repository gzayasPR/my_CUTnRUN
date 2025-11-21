#!/bin/bash
#SBATCH --time=24:00:00
# Usage: compare_seacr_peaks.sh <WT_peak_file> <KO_peak_file> <out_dir>


# --- Inputs ---
WT_PEAK=$1
KO_PEAK=$2
OUT_DIR=$3

# OUT_DIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run1/5.PeakComparison/KO_H3_vs_WT_H3/PEAKS_SEACR"
# WT_PEAK="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run1/4.PeakCalling/WT_H3/SEACR/WT_H3_SEACR_peaks.stringent.bed"
# KO_PEAK="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run1/4.PeakCalling/KO_H3/SEACR/KO_H3_SEACR_peaks.stringent.bed"

mkdir -p "$OUT_DIR"

ml bedtools
ml R
WT_BASENAME=$(basename "$WT_PEAK" .bed)
KO_BASENAME=$(basename "$KO_PEAK" .bed)

# --- Output files ---
SHARED="$OUT_DIR/shared_peaks.bed"
WT_ONLY="$OUT_DIR/${WT_BASENAME}_only.bed"
KO_ONLY="$OUT_DIR/${KO_BASENAME}_only.bed"
MERGED="$OUT_DIR/merged_peaks.bed"
SCORES="$OUT_DIR/WT_vs_KO_enrichment.tsv"
R_PLOT="$OUT_DIR/WT_vs_KO_enrichment_scatter.pdf"

# --- Logging parameters ---
LOGFILE="${OUT_DIR}/compare_seacr_peaks_run.log"
echo "📄 Logging run parameters to $LOGFILE"
{
  echo "Date: $(date)"
  echo "WT SEACR Peak File: $WT_PEAK"
  echo "KO SEACR Peak File: $KO_PEAK"
  echo "Output Directory: $OUT_DIR"
  echo "Shared Peaks Output: $SHARED"
  echo "WT-only Peaks Output: $WT_ONLY"
  echo "KO-only Peaks Output: $KO_ONLY"
  echo "Merged Peaks Output: $MERGED"
  echo "WT vs KO Enrichment Table: $SCORES"
  echo "Enrichment Scatterplot: $R_PLOT"
  echo ""
} > "$LOGFILE"

cat "$LOGFILE"
# --- Step 1: Find overlapping and unique peaks ---
echo "🔍 Comparing peak overlaps..."
bedtools intersect -a "$WT_PEAK" -b "$KO_PEAK" -wa -wb > "$SHARED"

awk 'BEGIN {OFS="\t"}
{
  chrom = $1
  new_start = ($2 > $8) ? $2 : $8
  new_end   = ($3 < $9) ? $3 : $9
  printf "%s\t%d\t%d", chrom, new_start, new_end

  # Print WT columns (4–6)
  for (i = 4; i <= 6; i++) printf "\t%s", $i

  # Print KO metadata columns (after $9)
  for (i = 10; i <= NF; i++) printf "\t%s", $i

  printf "\n"
}' "$SHARED" > "$OUT_DIR/shared_peaks_union.bed"


bedtools intersect -v -a "$WT_PEAK" -b "$KO_PEAK" > "$WT_ONLY"
bedtools intersect -v -a "$KO_PEAK" -b "$WT_PEAK" > "$KO_ONLY"


# # --- Step 2: Merge and prepare enrichment comparison table ---
# echo "📊 Merging enrichment scores for shared peaks..."

# awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' "$WT_PEAK" > "$OUT_DIR/wt_score.tmp"
# awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' "$KO_PEAK" > "$OUT_DIR/ko_score.tmp"

# sort -k1,1 -k2,2n "$OUT_DIR/wt_score.tmp" > "$OUT_DIR/wt_sorted.bed"
# sort -k1,1 -k2,2n "$OUT_DIR/ko_score.tmp" > "$OUT_DIR/ko_sorted.bed"

# bedtools intersect -a "$OUT_DIR/wt_sorted.bed" -b "$OUT_DIR/ko_sorted.bed" -wa -wb > "$OUT_DIR/intersect_scored.tsv"

# # Format: chr  start  end  wt_score  chr  start  end  ko_score
# awk 'BEGIN{OFS="\t"} {print $1":"$2"-"$3, $4, $8}' "$OUT_DIR/intersect_scored.tsv" > "$SCORES"

# # --- Step 3: Generate scatterplot in R ---
# echo "🧬 Plotting WT vs KO enrichment..."
# Rscript - <<EOF
# df <- read.table("$SCORES", header=FALSE)
# colnames(df) <- c("Region", "WT", "KO")

# pdf("$R_PLOT", width=6, height=6)
# plot(df\$WT, df\$KO,
#      xlab = "WT Enrichment",
#      ylab = "KO Enrichment",
#      main = "WT vs KO SEACR Enrichment",
#      col = "darkblue", pch = 16, log = "xy")
# abline(0,1,col="red", lty=2)
# dev.off()
# EOF

# echo "✅ Comparison complete."
# echo " - Shared peaks: $SHARED"
# echo " - WT-only peaks: $WT_ONLY"
# echo " - KO-only peaks: $KO_ONLY"
# echo " - Enrichment scatterplot: $R_PLOT"
