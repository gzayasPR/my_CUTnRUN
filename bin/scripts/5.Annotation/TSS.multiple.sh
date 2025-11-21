#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=mateescu
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G

# Usage:
# bash plot_tss_profiles.sh <OUT_TSS_DIR> <GENE_BED> <CSV_PEAKS> <WT_H3_BW> <KO_H3_BW> <WT_IgG_BW> <KO_IgG_BW>


GENE_BED="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/references/Mus_musculus_GRCm39/Ensembl/111/Annotation/Genes/Mus_musculus.GRCm39.111_genes.bed"
CSV_PEAKS="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/Enrichment/combined_peaks_log2FC.csv"
BW_DIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization"
metadata="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/metadata/Project6_H3K27me3.csv"
OUT_TSS_DIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/5.PeakComparison/WT_vs_KO/TSS_Compare"


WT1="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT1_H3K27me3/WT1_H3K27me3_mouse_scaled.bw"
WT2="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT2_H3K27me3/WT2_H3K27me3_mouse_scaled.bw"
WT3="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT3_H3K27me3/WT3_H3K27me3_mouse_scaled.bw"
KO1="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO1_H3K27me3/KO1_H3K27me3_mouse_scaled.bw"
KO2="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO2_H3K27me3/KO2_H3K27me3_mouse_scaled.bw"
KO3="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO3_H3K27me3/KO3_H3K27me3_mouse_scaled.bw"
KO_IgG="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO_IgG/KO_IgG_mouse_scaled.bw"
WT_IgG="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT_IgG/WT_IgG_mouse_scaled.bw"

ml  deeptools/3.5.2 
ml bedtools
ml R

mkdir -p $OUT_TSS_DIR
cd $OUT_TSS_DIR
##################################################################################
##################################################################################

bigwigAverage -b "$WT1" "$WT2" "$WT3" -p 7 -o WT_Average.bw
bigwigAverage -b "$KO1" "$KO2" "$KO3" -p 7 -o KO_Average.bw
computeMatrix reference-point -R "$GENE_BED" \
  -S WT_Average.bw KO_Average.bw   \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
  -p 8 \
  -o matrix_genes_TSS_ave.gz

plotHeatmap -m matrix_genes_TSS_ave.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT" "KO" \
  --legendLocation "upper-right" \
  --regionsLabel "all genes" \
    --heatmapHeight 14 \
  -out Average_KOvsWT_genes.facet.sample.png

plotHeatmap -m matrix_genes_TSS_ave.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS  --samplesLabel "WT" "KO" \
  --legendLocation "none" \
    --perGroup \
  --regionsLabel "all genes" \
    --heatmapHeight 14 \
  -out Average_KOvsWT_genes.no.facet.png



computeMatrix reference-point -R "$GENE_BED" \
  -S "$WT1" "$WT2" "$WT3"  "$KO1" "$KO2" "$KO3" "$WT_IgG" "$KO_IgG"  \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
  -p 8 \
  -o matrix_genes_TSS_raw.gz


# Heatmaps: WT and KO together
plotHeatmap -m matrix_genes_TSS_raw.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT1" "WT2" "WT3" "KO1" "KO2" "KO3" "WT IgG" "KO IgG" \
  --legendLocation "upper-right" \
  --regionsLabel "all genes" \
    --heatmapHeight 14 \
  -out KOvsWT_genes.facet.sample.png

# Heatmaps: WT and KO together
plotHeatmap -m matrix_genes_TSS_raw.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT1" "WT2" "WT3" "KO1" "KO2" "KO3" "WT IgG" "KO IgG" \
  --legendLocation "none" \
    --perGroup \
  --regionsLabel "all genes" \
    --heatmapHeight 14 \
  -out KOvsWT_genes.no.facet.png


computeMatrix reference-point -R "$GENE_BED" \
  -S "$WT1" "$WT2" "$WT3"  "$KO1" "$KO2" "$KO3"  \
  --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 \
  --skipZeros \
  -p 8 \
  -o matrix_genes_TSS_raw_no_IgG.gz


# Heatmaps: WT and KO together
plotHeatmap -m matrix_genes_TSS_raw_no_IgG.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT1" "WT2" "WT3" "KO1" "KO2" "KO3" \
  --legendLocation "upper-right" \
  --regionsLabel "all genes" \
    --heatmapHeight 14 \
  -out KOvsWT_genes.facet.sample_no_IgG.png

# Heatmaps: WT and KO together
plotHeatmap -m matrix_genes_TSS_raw_no_IgG.gz \
  --colorMap RdBu --zMin auto --zMax auto \
  --refPointLabel TSS --samplesLabel "WT1" "WT2" "WT3" "KO1" "KO2" "KO3" \
  --legendLocation "none" \
    --perGroup \
  --regionsLabel "all genes" \
    --heatmapHeight 14 \
  -out KOvsWT_genes.no.facet_no_IgG.png
