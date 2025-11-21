#!/bin/bash
# motif_analysis.sh
# Usage: bash motif_analysis.sh <WT_ONLY_BED> <OUTPUT_DIR> <GENOME>

set -euo pipefail

WT_ONLY_BED=$1
OUT_DIR=$2
GENOME=$3   # e.g., mm10 or hg38

mkdir -p "$OUT_DIR"

# --- Check HOMER ---
if ! command -v findMotifsGenome.pl &> /dev/null; then
  echo "❌ HOMER is not installed or not in PATH"
  exit 1
fi

# --- Center the peaks around midpoint (optional but preferred) ---
echo "📐 Centering peaks..."
CENTERED_BED="$OUT_DIR/centered_peaks.bed"
awk -F'\t' 'BEGIN{OFS="\t"} {mid=int(($2+$3)/2); print $1, mid-100, mid+100}' "$WT_ONLY_BED" > "$CENTERED_BED"

# --- Run HOMER motif enrichment ---
echo "🔍 Running HOMER motif analysis on WT-specific peaks..."
findMotifsGenome.pl "$CENTERED_BED" "$GENOME" "$OUT_DIR/homer_results" -size 200 -mask

echo "✅ Motif analysis complete. Results saved to: $OUT_DIR/homer_results"