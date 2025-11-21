#!/bin/bash
# Usage: bash summarize_alignment_metrics.sh <project_name> <multiqc_counts_file> <alignment_dir> <output_file>

PROJECT=$1
SEQ_COUNT_FILE=$2
ALIGN_DIR=$3
OUT_FILE=$4

if [[ -z "$PROJECT" || -z "$SEQ_COUNT_FILE" || -z "$ALIGN_DIR" || -z "$OUT_FILE" ]]; then
  echo "Usage: bash summarize_alignment_metrics.sh <project_name> <multiqc_counts_file> <alignment_dir> <output_file>"
  exit 1
fi

LOG_FILE="${OUT_FILE%.tsv}.log"
{
  echo "📄 summarize_alignment_metrics.sh log"
  echo "Date: $(date)"
  echo "User: $(whoami)"
  echo
  echo "Parameters:"
  echo "  Project       = $PROJECT"
  echo "  MultiQC File  = $SEQ_COUNT_FILE"
  echo "  Alignment Dir = $ALIGN_DIR"
  echo "  Output File   = $OUT_FILE"
  echo
} > "$LOG_FILE"

{
  echo "📂 Raw contents of MultiQC fastqc_sequence_counts_plot.txt:"
  echo "----------------------------------------------------------"
  cat "$SEQ_COUNT_FILE"
  echo
} >> "$LOG_FILE"

echo -e "Project\tSample\tRaw Reads\tTrimmed Reads\tMapped Mouse Reads\tMouse Alignment Rate\tClean Reads\tEcoli Reads\tPortion of Ecoli" > "$OUT_FILE"

declare -A TOTAL_COUNTS
while IFS=$'\t' read -r sample unique dup; do
  [[ "$sample" == "Sample" ]] && continue
  total=$(awk "BEGIN { printf \"%.0f\", $unique + $dup }")
  TOTAL_COUNTS["$sample"]=$total
done < "$SEQ_COUNT_FILE"

{
  echo "📊 Parsed MultiQC totals (sample → total reads):"
  for s in "${!TOTAL_COUNTS[@]}"; do
    echo "  $s = ${TOTAL_COUNTS[$s]}"
  done
  echo
} | tee -a "$LOG_FILE"

for sample_dir in "$ALIGN_DIR"/*; do
  [[ -d "$sample_dir" ]] || continue
  SAMPLE=$(basename "$sample_dir")
  echo "🔍 Processing $SAMPLE" | tee -a "$LOG_FILE"

  RAW_R1S=$(printf '%s\n' "${!TOTAL_COUNTS[@]}" | grep -E "^${SAMPLE}(_[0-9]+)?_R1$")
  RAW_R2S=$(printf '%s\n' "${!TOTAL_COUNTS[@]}" | grep -E "^${SAMPLE}(_[0-9]+)?_R2$")
  TRIM_R1S=$(printf '%s\n' "${!TOTAL_COUNTS[@]}" | grep -E "^${SAMPLE}(_[0-9]+)?_R1_001_val_1$")
  TRIM_R2S=$(printf '%s\n' "${!TOTAL_COUNTS[@]}" | grep -E "^${SAMPLE}(_[0-9]+)?_R2_001_val_2$")

  {
    echo "🧪 FASTQ match summary for $SAMPLE:"
    echo "  RAW_R1S:   $RAW_R1S"
    echo "  RAW_R2S:   $RAW_R2S"
    echo "  TRIM_R1S:  $TRIM_R1S"
    echo "  TRIM_R2S:  $TRIM_R2S"
    echo
  } | tee -a "$LOG_FILE"

  RAW_TOTAL=0
  for s in $RAW_R1S; do
    base=${s%_R1}
    r1=${TOTAL_COUNTS["${base}_R1"]}
    r2=${TOTAL_COUNTS["${base}_R2"]}
    if [[ "$r1" != "$r2" ]]; then
      echo "⚠️ Mismatch: ${base}_R1 and _R2 counts do not match → $r1 vs $r2" | tee -a "$LOG_FILE"
    fi
    RAW_TOTAL=$(awk -v a="$RAW_TOTAL" -v b="$r1" 'BEGIN { print a + b }')
  done

  TRIM_TOTAL=0
  for s in $TRIM_R1S; do
    base=${s%_R1_001_val_1}
    r1=${TOTAL_COUNTS["${base}_R1_001_val_1"]}
    r2=${TOTAL_COUNTS["${base}_R2_001_val_2"]}
    if [[ "$r1" != "$r2" ]]; then
      echo "⚠️ Mismatch: ${base}_val_1 and _val_2 counts do not match → $r1 vs $r2" | tee -a "$LOG_FILE"
    fi
    TRIM_TOTAL=$(awk -v a="$TRIM_TOTAL" -v b="$r1" 'BEGIN { print a + b }')
  done

  FLAGSTAT="${sample_dir}/mouse_alignment/${SAMPLE}_mouse_total.flagstat.txt"
  UNIQUE_FILE="${sample_dir}/mouse_alignment/${SAMPLE}_mouse_unique_counts.txt"
  ECOLI_FILE="${sample_dir}/ecoli_alignment/${SAMPLE}_ecoli_counts.txt"

  MAPPED_MOUSE_RAW=0
  if [[ -f "$FLAGSTAT" ]]; then
    MAPPED_MOUSE_RAW=$(grep " mapped (" "$FLAGSTAT" | head -1 | awk '{print $1}')
  else
    echo "⚠️ Missing mouse_total.flagstat.txt for $SAMPLE" | tee -a "$LOG_FILE"
  fi

  CLEAN_READS=0
  [[ -f "$UNIQUE_FILE" ]] && CLEAN_READS=$(cut -f2 "$UNIQUE_FILE")

  ECOLI_READS=0
  [[ -f "$ECOLI_FILE" ]] && ECOLI_READS=$(cut -f2 "$ECOLI_FILE")

  # Convert and calculate
  RAW_READS=$(awk -v x="$RAW_TOTAL" 'BEGIN { printf "%d", x * 2 }')
  TRIM_READS=$(awk -v x="$TRIM_TOTAL" 'BEGIN { printf "%d", x * 2 }')
  MAPPED_MOUSE=$(awk -v x="$MAPPED_MOUSE_RAW" 'BEGIN { printf "%d", x / 2 }')
  CLEAN_READS=$(awk -v x="$CLEAN_READS" 'BEGIN { printf "%d", x / 2 }')
  ECOLI_READS=$(awk -v x="$ECOLI_READS" 'BEGIN { printf "%d", x / 2 }')

  if [[ "$TRIM_TOTAL" -gt 0 ]]; then
    MOUSE_ALIGN_RATE=$(awk "BEGIN { printf \"%.2f\", ($MAPPED_MOUSE / $TRIM_TOTAL) * 100 }")
  else
    MOUSE_ALIGN_RATE="NA"
  fi

  if [[ "$CLEAN_READS" -gt 0 ]]; then
    PORTION_ECOLI=$(awk "BEGIN { printf \"%.2f\", ($ECOLI_READS / $CLEAN_READS) * 100 }")
  else
    PORTION_ECOLI="NA"
  fi

  echo -e "${PROJECT}\t${SAMPLE}\t${RAW_READS}\t${TRIM_READS}\t${MAPPED_MOUSE}\t${MOUSE_ALIGN_RATE}%\t${CLEAN_READS}\t${ECOLI_READS}\t${PORTION_ECOLI}" >> "$OUT_FILE"
done

echo "✅ Summary written to: $OUT_FILE" | tee -a "$LOG_FILE"
