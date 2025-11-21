#!/bin/bash
#SBATCH --time=72:00:00

# Usage:
# bash plot_tss_profiles.sh <OUT_TSS_DIR> <GENE_BED> <CSV_PEAKS> <WT_H3_BW> <KO_H3_BW> <WT_IgG_BW> <KO_IgG_BW>

OUT_TSS_DIR=$1
GENE_BED=$2
CSV_PEAKS=$3
WT_H3_BW=$4
KO_H3_BW=$5
WT_IgG_BW=$6
KO_IgG_BW=$7

# GENE_BED="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/references/Mus_musculus_GRCm39/Ensembl/111/Annotation/Genes/Mus_musculus.GRCm39.111_genes.bed"
# CSV_PEAKS="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project1_CREB1/5.PeakComparison/WT_vs_KO/Enrichment/combined_peaks_log2FC.csv"
# WT_H3_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project1_CREB1/3.Normalization/WT_CREB1/WT_CREB1_mouse_scaled.bw"
# KO_H3_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project1_CREB1/3.Normalization/KO_CREB1/KO_CREB1_mouse_scaled.bw"
# WT_IgG_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project1_CREB1/3.Normalization/WT_IgG/WT_IgG_mouse_scaled.bw"
# KO_IgG_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project1_CREB1/3.Normalization/KO_IgG/KO_IgG_mouse_scaled.bw"
# OUT_TSS_DIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project1_CREB1/5.PeakComparison/WT_vs_KO/TSS_Compare"

mkdir -p "$OUT_TSS_DIR"
cd "$OUT_TSS_DIR"

LOGFILE="${OUT_TSS_DIR}/TSS.log"
{
  echo "📄 TSS Profiling Log - $(date)"
  echo "GENE_BED: $GENE_BED"
  echo "CSV_PEAKS: $CSV_PEAKS"
  echo "WT_H3_BW: $WT_H3_BW"
  echo "KO_H3_BW: $KO_H3_BW"
  echo "WT_IgG_BW: $WT_IgG_BW"
  echo "KO_IgG_BW: $KO_IgG_BW"
} > "$LOGFILE"
cat "$LOGFILE"

ml deeptools
ml bedtools

# -------------------------
# Preprocessing shared peaks
# -------------------------
echo "🔹 Extracting shared peaks from CSV..."
shared_peaks="${OUT_TSS_DIR}/combined_shared_peaks.bed"
awk -F',' 'NR > 1 { OFS = "\t"; print $1, $2, $3 }' "$CSV_PEAKS" > "$shared_peaks"

echo "🔹 Intersecting shared peaks with gene BED..."
bedtools intersect -a "$GENE_BED" -b "$shared_peaks" > TSS_shared.bed

# --------------------------------------------
# Section 1: Raw H3 (WT vs KO) - All genes
# --------------------------------------------
echo "🔹 Section 1: Raw H3 signal around all gene TSS..."
computeMatrix reference-point \
  -R "$GENE_BED" \
  -S "$WT_H3_BW" "$KO_H3_BW" \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros -p 7 \
  -o matrix_allTSS_WTvsKO_raw.gz

plotProfile -m matrix_allTSS_WTvsKO_raw.gz \
  --perGroup --colors blue red \
  --regionsLabel "All Genes" --samplesLabel "WT" "KO" \
  -o metaprofile_allTSS_WTvsKO_raw.png

# --------------------------------------------
# Section 2: Raw H3 (WT vs KO) - Peak genes
# --------------------------------------------
echo "🔹 Section 2: Raw H3 signal around SEACR peak-associated TSS..."
computeMatrix reference-point \
  -R TSS_shared.bed \
  -S "$WT_H3_BW" "$KO_H3_BW" \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 \
  -o matrix_peaks_TSS_WTvsKO_raw.gz

plotProfile -m matrix_peaks_TSS_WTvsKO_raw.gz \
  --perGroup --colors blue red \
  --regionsLabel "Genes w/ SEACR Peaks" --samplesLabel "WT" "KO" \
  -o metaprofile_peaks_TSS_WTvsKO_raw.png

# --------------------------------------------
# Section 3: H3 + IgG - All genes
# --------------------------------------------
echo "🔹 Section 3: H3 + IgG signal around all gene TSS..."
computeMatrix reference-point \
  -R "$GENE_BED" \
  -S "$WT_H3_BW" "$KO_H3_BW" "$WT_IgG_BW" "$KO_IgG_BW" \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros -p 7 \
  -o matrix_allTSS.gz

plotProfile -m matrix_allTSS.gz \
  --perGroup --colors blue red green yellow \
  --regionsLabel "All Genes" \
  -o metaprofile_allTSS.png

# --------------------------------------------
# Section 4: H3 + IgG - Peak genes
# --------------------------------------------
echo "🔹 Section 4: H3 + IgG signal at SEACR peak-associated TSS..."
computeMatrix reference-point \
  -R TSS_shared.bed \
  -S "$WT_H3_BW" "$KO_H3_BW" "$WT_IgG_BW" "$KO_IgG_BW" \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 \
  -o matrix_peaks_TSS.gz

plotProfile -m matrix_peaks_TSS.gz \
  --perGroup --colors blue red green yellow \
  --regionsLabel "Genes w/ SEACR Peaks" \
  -o metaprofile_peaks_TSS.png

# --------------------------------------------
# Section 5: IgG-subtracted signal
# --------------------------------------------
echo "🔹 Section 5: Subtracting IgG from H3..."
bigwigCompare -b1 "$WT_H3_BW" -b2 "$WT_IgG_BW" --operation subtract -p 7 -o WT_H3_minus_IgG.bw
bigwigCompare -b1 "$KO_H3_BW" -b2 "$KO_IgG_BW" --operation subtract -p 7 -o KO_H3_minus_IgG.bw

# All genes
echo "🔹 Section 5A: Subtracted signal around all genes..."
computeMatrix reference-point \
  -R "$GENE_BED" \
  -S WT_H3_minus_IgG.bw KO_H3_minus_IgG.bw \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 -o matrix_allTSS_subtract.gz

plotHeatmap -m matrix_allTSS_subtract.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --heatmapHeight 12 --heatmapWidth 6 \
  --refPointLabel "TSS" --sortRegions descend --sortUsing mean \
  --plotTitle "CUT&RUN Signal Around All Gene TSS (IgG Subtracted)" \
  -out heatmap_allTSS_withProfile.png

plotProfile -m matrix_allTSS_subtract.gz \
  --perGroup --colors blue red \
  --regionsLabel "All Genes" --samplesLabel "WT" "KO" \
  -o metaprofile_allTSS_subtracted.png

# Peak genes
echo "🔹 Section 5B: Subtracted signal at SEACR peak-associated genes..."
computeMatrix reference-point \
  -R TSS_shared.bed \
  -S WT_H3_minus_IgG.bw KO_H3_minus_IgG.bw \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 -o matrix_peaks_TSS_subtract.gz

plotHeatmap -m matrix_peaks_TSS_subtract.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --heatmapHeight 12 --heatmapWidth 6 \
  --refPointLabel "TSS" --sortRegions descend --sortUsing mean \
  --plotTitle "CUT&RUN Signal at SEACR Peak-Associated TSS (IgG Subtracted)" \
  -out heatmap_peaks_TSS_withProfile.png

plotProfile -m matrix_peaks_TSS_subtract.gz \
  --perGroup --colors blue red \
  --regionsLabel "Genes w/ SEACR Peaks" --samplesLabel "WT" "KO" \
  -o metaprofile_peaks_TSS_subtracted.png

# --------------------------------------------
# Section 6: Combined TSS sets (optional bonus)
# --------------------------------------------
echo "🔹 Section 6: Combined view of all and peak genes (optional)..."
computeMatrix reference-point \
  -R "$GENE_BED" TSS_shared.bed \
  -S WT_H3_minus_IgG.bw KO_H3_minus_IgG.bw \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 -o matrix_peaks_subtract.gz

plotHeatmap -m matrix_peaks_subtract.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --heatmapHeight 12 --heatmapWidth 6 \
  --refPointLabel "TSS" --sortRegions descend --sortUsing mean \
  --plotTitle "CUT&RUN Signal: All vs SEACR Peak-Associated TSS" \
  -out heatmap_peaks_genes_TSS_withProfile.png

echo "✅ All steps completed."