#!/bin/bash
#SBATCH --time=72:00:00

# Usage:

# ==============================================
# Usage:
# bash plot_tss_profiles_from_bam.sh <OUT_TSS_DIR> <GENE_BED> <CSV_PEAKS> <WT_H3_BAM> <KO_H3_BAM> <WT_IgG_BAM> <KO_IgG_BAM>
# ==============================================

OUT_TSS_DIR=$1
GENE_BED=$2
CSV_PEAKS=$3
WT_H3_BAM=$4
KO_H3_BAM=$5
WT_IgG_BAM=$6
KO_IgG_BAM=$7
# OUT_TSS_DIR=/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/TSS_Compare_BAM
# GENE_BED=/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/references/Mus_musculus_GRCm39/Ensembl/111/Annotation/Genes/Mus_musculus.GRCm39.111_genes.bed
# CSV_PEAKS=/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/Enrichment/combined_peaks_log2FC.csv
# WT_H3_BAM=/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/WT_Dpy30/WT_Dpy30_mouse_filtered.bam
# KO_H3_BAM=/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/KO_Dpy30/KO_Dpy30_mouse_filtered.bam
# WT_IgG_BAM=/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/WT_IgG/WT_IgG_mouse_filtered.bam
# KO_IgG_BAM=/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/3.Normalization/KO_IgG/KO_IgG_mouse_filtered.bam



mkdir -p "$OUT_TSS_DIR"
cd "$OUT_TSS_DIR"

ml deeptools
ml bedtools
ml R

# ------------------------------------------------
# Step 1: Logging
# ------------------------------------------------
LOGFILE="$OUT_TSS_DIR/TSS.log"
{
  echo "📄 TSS Profiling Log - $(date)"
  echo "GENE_BED: $GENE_BED"
  echo "CSV_PEAKS: $CSV_PEAKS"
  echo "WT_H3_BAM: $WT_H3_BAM"
  echo "KO_H3_BAM: $KO_H3_BAM"
  echo "WT_IgG_BAM: $WT_IgG_BAM"
  echo "KO_IgG_BAM: $KO_IgG_BAM"
} > "$LOGFILE"

# ------------------------------------------------
# Step 2: Convert BAM to BigWig (CPM-normalized, .bigwig)
# ------------------------------------------------
echo "🚧 Converting BAM to BigWig (.bigwig extension)..."

bamCoverage -b "$WT_H3_BAM" -o WT_H3.bigwig --normalizeUsing CPM \
  --binSize 25 --extendReads --centerReads -p 6

bamCoverage -b "$KO_H3_BAM" -o KO_H3.bigwig --normalizeUsing CPM \
  --binSize 25 --extendReads --centerReads -p 6

bamCoverage -b "$WT_IgG_BAM" -o WT_IgG.bigwig --normalizeUsing CPM \
  --binSize 25 --extendReads --centerReads -p 6

bamCoverage -b "$KO_IgG_BAM" -o KO_IgG.bigwig --normalizeUsing CPM \
  --binSize 25 --extendReads --centerReads -p 6

# ------------------------------------------------
# Step 3: Extract shared peaks
# ------------------------------------------------
echo "📌 Extracting shared SEACR peaks..."
shared_peaks="${OUT_TSS_DIR}/combined_shared_peaks.bed"
awk -F',' 'NR > 1 { OFS = "\t"; print $1, $2, $3 }' "$CSV_PEAKS" > "$shared_peaks"
bedtools intersect -a "$GENE_BED" -b "$shared_peaks" > TSS_shared.bed

# ------------------------------------------------
# Step 4: IgG-subtracted BigWigs
# ------------------------------------------------
echo "➖ Subtracting IgG from signal..."
bigwigCompare -b1 WT_H3.bigwig -b2 WT_IgG.bigwig --operation subtract -p 6 -o WT_H3_minus_IgG.bigwig
bigwigCompare -b1 KO_H3.bigwig -b2 KO_IgG.bigwig --operation subtract -p 6 -o KO_H3_minus_IgG.bigwig

# ------------------------------------------------
# Step 5: TSS matrices and plots
# ------------------------------------------------

echo "📊 Generating TSS matrix and profile plots..."

computeMatrix reference-point -R "$GENE_BED" TSS_shared.bed \
  -S WT_H3.bigwig KO_H3.bigwig WT_IgG.bigwig KO_IgG.bigwig \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros -p 6 -o matrix_all_peaks_TSS_raw.gz

plotProfile -m matrix_all_peaks_TSS_raw.gz \
  --samplesLabel "WT" "KO" "WT_IgG" "KO_IgG" --regionsLabel "All Genes" "peaks" \
  --perGroup \
  -o profile_raw_all_genes.png

computeMatrix reference-point -R "$GENE_BED" TSS_shared.bed \
  -S WT_H3_minus_IgG.bigwig KO_H3_minus_IgG.bigwig \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros -p 6 -o matrix_all_peaks_TSS_subtracted.gz

plotProfile -m matrix_all_peaks_TSS_subtracted.gz \
  --samplesLabel "WT" "KO" --regionsLabel "All Genes" "peaks" \
  --legendLocation "upper-right" \
  --perGroup \
  -o profile_raw_all_genes_subtracted.png

plotHeatmap -m matrix_all_peaks_TSS_subtracted.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "upper-right" \
  -out heatmap_all_genes_subtracted.png

plotHeatmap -m matrix_all_peaks_TSS_subtracted.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --perGroup \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "upper-right" \
  --plotTitle "SEACR Peaks (WT vs KO)" \
  -out heatmap_peaks_subtracted.png

computeMatrix reference-point -R "$GENE_BED"  \
  -S WT_H3.bigwig KO_H3.bigwig \
      --skipZeros \
      --nanAfterEnd \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  -p 7 -o matrix_genes.gz
ml  deeptools/3.5.2 
plotHeatmap -m matrix_genes.gz \
  --colorMap RdBu --zMin auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "none" \
  --regionsLabel "genes" \
    --yMin 0 \
  --zMin 0 \
  --heatmapHeight 12 \
  -out WTvsKO_genes.png
echo "✅ DONE — All profiles and heatmaps generated from BAMs!"
