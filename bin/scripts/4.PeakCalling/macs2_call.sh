#!/bin/bash
#SBATCH --time=12:00:00
set -euo pipefail

# --- Inputs ---
sample=$1
BAM=$2
OUT_DIR=$3
CPUS=$4

# --- Load module ---
module load macs

# --- Run MACS2 ---
macs2 callpeak \
  -t "$BAM" \
  -f BAMPE \
  -g mm \
  -n "$sample" \
  --outdir "$OUT_DIR" \
  --keep-dup all \
  --nomodel --shift 0 --extsize 200 \
  --call-summits

echo "[✓] MACS2 peak calling complete for $sample"
