#!/bin/bash
#SBATCH --time=72:00:00

# Usage:
# bash plot_tss_profiles.sh <OUT_TSS_DIR> <GENE_BED> <CSV_PEAKS> <WT_H3_BW> <KO_H3_BW> <WT_IgG_BW> <KO_IgG_BW>

# OUT_TSS_DIR=$1
# GENE_BED=$2
# CSV_PEAKS=$3
# WT_H3_BW=$4
# KO_H3_BW=$5
# WT_IgG_BW=$6
# KO_IgG_BW=$7

GENE_BED="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/references/Mus_musculus_GRCm39/Ensembl/111/Annotation/Genes/Mus_musculus.GRCm39.111_genes.bed"
CSV_PEAKS="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/Enrichment/combined_peaks_log2FC.csv"
WT_H3_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/WT_Dpy30/WT_Dpy30_mouse_scaled.bw"
KO_H3_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/KO_Dpy30/KO_Dpy30_mouse_scaled.bw"
WT_IgG_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/WT_IgG/WT_IgG_mouse_scaled.bw"
KO_IgG_BW="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/KO_IgG/KO_IgG_mouse_scaled.bw"
OUT_TSS_DIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/TSS_Compare"


mkdir -p "$OUT_TSS_DIR"
cd "$OUT_TSS_DIR"

# ------------------------------------------------
# Step 1: Logging
# ------------------------------------------------
echo "📝 Logging parameters..."
LOGFILE="$OUT_TSS_DIR/TSS.log"
{
  echo "📄 TSS Profiling Log - $(date)"
  echo "GENE_BED: $GENE_BED"
  echo "CSV_PEAKS: $CSV_PEAKS"
  echo "WT_H3_BW: $WT_H3_BW"
  echo "KO_H3_BW: $KO_H3_BW"
  echo "WT_IgG_BW: $WT_IgG_BW"
  echo "KO_IgG_BW: $KO_IgG_BW"
  echo "OUT_TSS_DIR: $OUT_TSS_DIR"
} > "$LOGFILE"

ml deeptools
ml bedtools
ml R

# ------------------------------------------------
# Step 2: Extract shared peaks
# ------------------------------------------------
echo "📌 Extracting shared SEACR peaks..."
shared_peaks="${OUT_TSS_DIR}/combined_shared_peaks.bed"
awk -F',' 'NR > 1 { OFS = "\t"; print $1, $2, $3 }' "$CSV_PEAKS" > "$shared_peaks"

bedtools intersect -a "$GENE_BED" -b "$shared_peaks" > TSS_shared.bed

# ------------------------------------------------
# Step 2b: Subtracted bigWigs
# ------------------------------------------------
echo "📌 Creating IgG-subtracted bigWig files..."
bigwigCompare -b1 "$WT_H3_BW" -b2 "$WT_IgG_BW" --operation subtract  -p 7 -o WT_H3_minus_IgG.bw
bigwigCompare -b1 "$KO_H3_BW" -b2 "$KO_IgG_BW" --operation subtract  -p 7 -o KO_H3_minus_IgG.bw

# bigwigCompare -b1 "$WT_H3_BW" -b2 "$WT_IgG_BW" --operation subtract  -o WT_H3_minus_IgG.bw
# bigwigCompare -b1 "$KO_H3_BW" -b2 "$KO_IgG_BW" --operation subtract  -o KO_H3_minus_IgG.bw

# ------------------------------------------------
# Step 3: Matrix + Plotting Setup
# ------------------------------------------------
echo "📊 Generating matrices for raw profiles..."
computeMatrix reference-point -R "$GENE_BED" TSS_shared.bed \
  -S "$WT_H3_BW" "$KO_H3_BW" "$WT_IgG_BW" "$KO_IgG_BW" \
  --skipZeros \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 -o matrix_all_peaks_TSS_raw.gz

plotProfile -m matrix_all_peaks_TSS_raw.gz \
  --samplesLabel "WT" "KO" "WT_IgG" "KO_IgG" --regionsLabel "All Genes" "peaks" \
  --perGroup \
  -o profile_raw_all_genes.png

# ------------------------------------------------
echo "📊 Generating matrices for IgG-subtracted signal..."
computeMatrix reference-point -R "$GENE_BED" TSS_shared.bed \
  -S WT_H3_minus_IgG.bw KO_H3_minus_IgG.bw \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
  -p 7 -o matrix_all_peaks_TSS_subtracted.gz

plotProfile -m  matrix_all_peaks_TSS_subtracted.gz \
  --samplesLabel "WT" "KO" --regionsLabel "All Genes" "peaks" \
  --legendLocation "upper-right" \
  --yMin 0 \
  --perGroup \
  -o profile_raw_all_genes_subtracted.png

echo "📊 Generating autoscaled heatmaps..."

# Heatmaps: WT and KO together
plotHeatmap -m matrix_all_peaks_TSS_subtracted.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --regionsLabel "All Genes" "peaks" \
  --legendLocation "upper-right" \
  --yMin 0 \
   --zMin 0 \
  -out heatmap_all_genes_subtracted.png

plotHeatmap -m matrix_all_peaks_TSS_subtracted.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --perGroup \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --regionsLabel "All Genes" "peaks" \
  --legendLocation "upper-right" \
  --plotTitle "SEACR Peaks (WT vs KO)" \
   --yMin 0 \
   --zMin 0 \
  -out heatmap_peaks_subtracted.png

echo "✅ All TSS plots and heatmaps completed!"


echo "📊 Generating matrices for single raw profiles..."
computeMatrix reference-point -R "$GENE_BED" \
  -S "$WT_H3_BW" "$KO_H3_BW" "$WT_IgG_BW" "$KO_IgG_BW" \
  --skipZeros \
    --skipZeros \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 -o matrix_genes_raw.gz

plotProfile -m matrix_genes_raw.gz \
  --samplesLabel "WT" "KO" "WT_IgG" "KO_IgG" --regionsLabel "All Genes" "peaks" \
  --perGroup \
  --yMin 0 \
  --zMin 0 \
  -o profile_raw_genes.png

echo "📊 Generating matrices for IgG-subtracted signal..."

computeMatrix reference-point -R "$GENE_BED"  \
  -S WT_H3_minus_IgG.bw KO_H3_minus_IgG.bw \
    --skipZeros \
      --nanAfterEnd \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 -o matrix_genes_subtracted.gz
ml  deeptools/3.5.2 
plotHeatmap -m matrix_genes_subtracted.gz \
  --colorMap RdBu --zMin auto --zMax 10 \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "none" \
  --regionsLabel "genes" \
    --yMin 0 \
  --zMin 0 \
  --heatmapHeight 14 \
  -out WTvsKO_genes.png

computeMatrix reference-point -R "TSS_shared.bed"  \
  -S WT_H3_minus_IgG.bw KO_H3_minus_IgG.bw \
  --nanAfterEnd \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
  -p 7 -o matrix_peaks_subtracted.gz

plotHeatmap -m matrix_peaks_subtracted.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "none" \
  --regionsLabel "peak genes" \
  --yMin 0 \
  --zMin 0 \
  --heatmapHeight 14 \
  -out WTvsKO_peaks.png


# ------------------------------------------------
echo "📊 Generating matrices for raw signal..."
computeMatrix reference-point -R "$GENE_BED" TSS_shared.bed \
  -S "$WT_H3_BW" "$KO_H3_BW"  \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
    -p 7 -o matrix_h3_peaks_TSS_raw.gz

  
plotProfile -m  matrix_h3_peaks_TSS_raw.gz \
  --samplesLabel "WT" "KO" --regionsLabel "All Genes" "peaks" \
  --legendLocation "upper-right" \
  --perGroup \
  --yMin 0 \
   --zMin 0 \
  -o profile_raw_h3_genes.png




echo "📊 Generating autoscaled heatmaps..."

# Heatmaps: WT and KO together
plotHeatmap -m matrix_h3_peaks_TSS_raw.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --regionsLabel "All Genes" "peaks" \
  --legendLocation "upper-right" \
  -out heatmap_h3_genes_raw.png

plotHeatmap -m matrix_h3_peaks_TSS_raw.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --perGroup \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --regionsLabel "All Genes" "peaks" \
  --legendLocation "upper-right" \
  --plotTitle "SEACR Peaks (WT vs KO)" \
  -out heatmap_h3_peaks_raw.png

echo "✅ All TSS plots and heatmaps completed!"

##################################################################################
##################################################################################
computeMatrix reference-point -R "$GENE_BED" \
  -S "$WT_H3_BW" "$KO_H3_BW"  \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
  -o matrix_genes_TSS_raw.gz


# Heatmaps: WT and KO together
plotHeatmap -m matrix_genes_TSS_raw.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "none" \
  --regionsLabel "all genes" \
  -out KOvsWT_genes.png

  computeMatrix reference-point -R "TSS_shared.bed" \
  -S "$WT_H3_BW" "$KO_H3_BW"  \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
    -p 7 -o matrix_peaks_TSS_raw.gz
# Heatmaps: WT and KO together
plotHeatmap -m  matrix_peaks_TSS_raw.gz\
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "none" \
  --regionsLabel "peaks" \
  -out KOvsWT_peaks.png




##################################################################################
##################################################################################
computeMatrix reference-point -R "$GENE_BED" \
  -S "$WT_H3_BW" "$KO_H3_BW"  \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
  -o matrix_genes_TSS_raw.gz


# Heatmaps: WT and KO together
plotHeatmap -m matrix_genes_TSS_raw.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "none" \
  --regionsLabel "all genes" \
  -out WTvsKO_genes.png

  computeMatrix reference-point -R "TSS_shared.bed" \
  -S "$WT_H3_BW" "$KO_H3_BW"  \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
    -p 7 -o matrix_peaks_TSS_raw.gz
# Heatmaps: WT and KO together
plotHeatmap -m  matrix_peaks_TSS_raw.gz\
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "none" \
  --regionsLabel "peaks" \
  -out  WTvsKO_genes.png
