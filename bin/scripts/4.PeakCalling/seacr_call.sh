#!/bin/bash
#SBATCH --time=12:00:00
set -euo pipefail

# --- Inputs ---
sample=$1
BW=$2
OUT_DIR=$3
CPUS=$4
proj_env=$5

source "${proj_env}"

# --- Paths ---
SEACR_PATH="${my_softwares}/SEACR/SEACR_1.3.sh"
bigwig_env="${my_softwares}/conda_envs/bigwig"
BG="${OUT_DIR}/${sample}.bedgraph"

# --- Load conda and activate environment ---
ml conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$bigwig_env"

# --- Convert BigWig to BedGraph ---
bigWigToBedGraph "$BW" "$BG"

# --- Load other modules ---
ml bedtools
ml R

# --- Run SEACR ---
bash "$SEACR_PATH" "$BG" 0.01 norm relaxed "${OUT_DIR}/${sample}_SEACR_peaks.bed"

echo "[✓] SEACR peak calling complete for $sample"
