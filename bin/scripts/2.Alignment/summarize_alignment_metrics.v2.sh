#!/usr/bin/env bash
# summarize_alignment_metrics.sh
# Usage:
#   bash summarize_alignment_metrics.sh <project_name> <multiqc_counts_file> <alignment_dir> <output_file>
#
# Inputs:
#   project_name            Label for output
#   multiqc_counts_file     MultiQC "fastqc_sequence_counts_plot.txt" (tab-delimited: Sample,Unique,Duplicates)
#   alignment_dir           Directory with sample subdirs (e.g., KO1, WT1, ...)
#   output_file             Path to write summary TSV
#
# Output TSV columns:
#   Project  Sample  RawReads  TrimmedReads  MouseAligned  EcoliAligned  TotalAligned  MouseClean  EcoliClean  TotalClean  MouseAlignRate  PortionOfEcoli
#
# Notes:
# - RawReads/TrimmedReads are computed from MultiQC (Unique+Duplicates per FASTQ), summed over lanes.
#   We sum R1 entries across lanes, sanity-check the matching R2, and then multiply by 2 to get total reads.
#   Raw pattern:   <SAMPLE>[_<lane>]?_R1 and _R2
#   Trim pattern:  <SAMPLE>[_<lane>]?_R1_001_val_1 and _R2_001_val_2
# - Aligned/Clean metrics are read primarily from <sample>/<sample>.reads.tsv (header-aware).
#   If missing, we fallback to:
#     mouse_alignment/<SAMPLE>_mouse_total.flagstat.txt              -> first "mapped (" field (reads; convert to pairs)
#     mouse_alignment/<SAMPLE>_mouse_unique_counts.txt (tab: name N) -> N (reads; convert to pairs)  [MouseClean]
#     ecoli_alignment/<SAMPLE>_ecoli_counts.txt (tab: name N)        -> N (reads; convert to pairs)  [EcoliClean]
#   TotalAligned = MouseAligned + EcoliAligned (on fallback)
#   TotalClean   = MouseClean + EcoliClean  (on fallback)
# - MouseAlignRate (%) = (MouseAligned / TrimmedReadPairs) * 100
# - PortionOfEcoli     =  EcoliClean / TotalClean
#

PROJECT="${1:-}"
SEQ_COUNT_FILE="${2:-}"
ALIGN_DIR="${3:-}"
OUT_FILE="${4:-}"

usage() {
  echo "Usage: bash $(basename "$0") <project_name> <multiqc_counts_file> <alignment_dir> <output_file>" >&2
  exit 1
}

[[ -z "$PROJECT" || -z "$SEQ_COUNT_FILE" || -z "$ALIGN_DIR" || -z "$OUT_FILE" ]] && usage
[[ ! -f "$SEQ_COUNT_FILE" ]] && { echo "ERROR: MultiQC counts file not found: $SEQ_COUNT_FILE" >&2; exit 2; }
[[ ! -d "$ALIGN_DIR" ]] && { echo "ERROR: Alignment directory not found: $ALIGN_DIR" >&2; exit 2; }

LOG_FILE="${OUT_FILE%.tsv}.log"
mkdir -p "$(dirname "$OUT_FILE")"

# Single log capturing stdout and stderr
exec > >(tee -a "$LOG_FILE") 2>&1

echo "📄 summarize_alignment_metrics.sh"
echo "Date: $(date)"
echo "User: $(whoami)"
echo "Project: $PROJECT"
echo "MultiQC: $SEQ_COUNT_FILE"
echo "AlignDir: $ALIGN_DIR"
echo "OutFile: $OUT_FILE"
echo

# === Output header ===
printf "Project\tSample\tRawReads\tTrimmedReads\tMouseAligned\tEcoliAligned\tTotalAligned\tMouseClean\tEcoliClean\tTotalClean\tMouseAlignRate\tPortionOfEcoli\n" > "$OUT_FILE"

# === Parse MultiQC counts into an associative array ===
# Expecting columns: Sample<TAB>Unique<TAB>Duplicates
declare -A TOTAL_COUNTS  # key=FASTQ sample id (e.g., KO1_R1, KO1_R1_001_val_1), value=(Unique+Duplicates)
while IFS=$'\t' read -r s u d rest; do
  # Skip header lines
  [[ "$s" == "Sample" ]] && continue
  # Defensive: require at least two numeric fields
  [[ -z "${u:-}" || -z "${d:-}" ]] && continue
  # Sum Unique+Duplicates (round to integer)
  total=$(awk -v a="${u:-0}" -v b="${d:-0}" 'BEGIN{printf "%d", a+b}')
  TOTAL_COUNTS["$s"]="$total"
done < "$SEQ_COUNT_FILE"

echo "✅ Loaded $((${#TOTAL_COUNTS[@]})) FASTQ count entries from MultiQC"
echo

# === Helper: sum read-pair counts across lanes by pattern ===
# We sum over R1 entries (by lane), check matching R2 exists and matches, then double to get total reads.
sum_pairs_by_pattern() {
  local sample="$1" r1_regex="$2" r2_swap_from="$3" r2_swap_to="$4"
  local sum_r1=0
  local mismatches=0

  # list all keys matching the r1 regex
  local k
  for k in "${!TOTAL_COUNTS[@]}"; do
    if [[ "$k" =~ $r1_regex ]]; then
      local r1="${TOTAL_COUNTS[$k]}"
      local mate="${k/$r2_swap_from/$r2_swap_to}"
      local r2="${TOTAL_COUNTS[$mate]:-}"
      if [[ -z "$r2" ]]; then
        echo "⚠️  [$sample] Missing mate for $k -> expected $mate" >&2
        mismatches=$((mismatches+1))
      elif [[ "$r1" != "$r2" ]]; then
        echo "⚠️  [$sample] R1/R2 mismatch for base ${k%$r2_swap_from}: $r1 vs $r2" >&2
        mismatches=$((mismatches+1))
      fi
      sum_r1=$(( sum_r1 + r1 ))
    fi
  done

  # Return two numbers via stdout: <read_pairs> <mismatches>
  # read_pairs = sum of R1 counts (each == number of pairs)
  echo "$sum_r1 $mismatches"
}

# === Helper: find and parse <sample>.reads.tsv if present ===
find_reads_tsv() {
  local sdir="$1"; local sname="$2"
  if [[ -f "$sdir/${sname}.reads.tsv" ]]; then
    echo "$sdir/${sname}.reads.tsv"
    return 0
  fi
  shopt -s nullglob
  local c=("$sdir"/*.reads.tsv)
  shopt -u nullglob
  if (( ${#c[@]} > 0 )); then
    echo "${c[0]}"; return 0
  fi
  return 1
}

# Parse .reads.tsv → one data line with all aligned/clean fields
parse_reads_tsv() {
  local fpath="$1"
  awk -v OFS="\t" '
    BEGIN { FS="\t"; have=0 }
    NR==1 {
      for(i=1;i<=NF;i++){ lc=tolower($i); gsub(/ /,"",lc); h[lc]=i }
      next
    }
    NR>1 {
      s  = (h["sample"]? $h["sample"] : $1)
      ma = (h["mousealigned"]? $h["mousealigned"] : $2)
      ea = (h["ecolialigned"]? $h["ecolialigned"] : $3)
      ta = (h["totalaligned"]? $h["totalaligned"] : $4)
      mc = (h["mouseclean"]? $h["mouseclean"] : $5)
      ec = (h["ecoliclean"]? $h["ecoliclean"] : $6)
      tc = (h["totalclean"]? $h["totalclean"] : $7)
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", s, ma, ea, ta, mc, ec, tc
      exit
    }
  ' "$fpath"
}

# Fallback loaders (reads→pairs conversions inside)
load_fallback_metrics() {
  local sdir="$1"; local sname="$2"
  local flag="$sdir/mouse_alignment/${sname}_mouse_total.flagstat.txt"
  local uniq="$sdir/mouse_alignment/${sname}_mouse_unique_counts.txt"
  local eco ="$sdir/ecoli_alignment/${sname}_ecoli_counts.txt"  # (tab: name N)

  local ma=0 mc=0 ec=0
  if [[ -f "$flag" ]]; then
    # first line with " mapped (" → first field = mapped reads
    ma=$(grep -m1 " mapped (" "$flag" | awk '{print $1}')
  else
    echo "⚠️  [$sname] Missing flagstat: $flag" >&2
  fi
  if [[ -f "$uniq" ]]; then
    mc=$(cut -f2 "$uniq")
  else
    echo "⚠️  [$sname] Missing unique counts: $uniq" >&2
  fi
  if [[ -f "$sdir/ecoli_alignment/${sname}_ecoli_counts.txt" ]]; then
    ec=$(cut -f2 "$sdir/ecoli_alignment/${sname}_ecoli_counts.txt")
  else
    echo "⚠️  [$sname] Missing ecoli counts." >&2
  fi

  # Convert reads→pairs (divide by 2; round)
  local ma_p mc_p ec_p
  ma_p=$(awk -v x="${ma:-0}" 'BEGIN{printf "%d", x/2}')
  mc_p=$(awk -v x="${mc:-0}" 'BEGIN{printf "%d", x/2}')
  ec_p=$(awk -v x="${ec:-0}" 'BEGIN{printf "%d", x/2}')
  local ta_p tc_p
  ta_p=$(( ma_p + 0 ))          # no ecoli_aligned in fallback; keep as mouse-only aligned
  tc_p=$(( mc_p + ec_p ))

  echo -e "${ma_p}\t0\t${ta_p}\t${mc_p}\t${ec_p}\t${tc_p}"
}

echo "▶️  Scanning samples under $ALIGN_DIR"
processed=0
skipped=0

shopt -s nullglob
for sample_dir in "$ALIGN_DIR"/*/; do
  [[ -d "$sample_dir" ]] || continue
  SAMPLE="$(basename "$sample_dir")"
  echo "— Processing $SAMPLE —"

  # 1) Raw & Trimmed read counts from MultiQC
  # Sum over R1 lanes for raw (pattern: ^SAMPLE(_\d+)?_R1$), mate must be _R2
  read RAW_PAIRS RAW_MIS <<< "$(sum_pairs_by_pattern "$SAMPLE" "^${SAMPLE}(_[0-9]+)?_R1$" "_R1" "_R2")"
  # Trimmed (pattern: ^SAMPLE(_\d+)?_R1_001_val_1$), mate _R2_001_val_2
  read TRIM_PAIRS TRIM_MIS <<< "$(sum_pairs_by_pattern "$SAMPLE" "^${SAMPLE}(_[0-9]+)?_R1_001_val_1$" "_R1_001_val_1" "_R2_001_val_2")"

  RAW_READS=$(( RAW_PAIRS * 2 ))
  TRIM_READS=$(( TRIM_PAIRS * 2 ))

  if (( RAW_PAIRS == 0 )); then
    echo "⚠️  [$SAMPLE] No raw FASTQ entries found in MultiQC. RawReads=0" >&2
  fi
  if (( TRIM_PAIRS == 0 )); then
    echo "⚠️  [$SAMPLE] No trimmed FASTQ entries found in MultiQC. TrimmedReads=0" >&2
  fi
  (( RAW_MIS > 0 ))  && echo "⚠️  [$SAMPLE] Found $RAW_MIS raw R1/R2 pairing mismatches" >&2
  (( TRIM_MIS > 0 )) && echo "⚠️  [$SAMPLE] Found $TRIM_MIS trimmed R1/R2 pairing mismatches" >&2

  # 2) Aligned/Clean metrics (prefer .reads.tsv, else fallback)
  rtsv=""
  if rtsv="$(find_reads_tsv "$sample_dir" "$SAMPLE")"; then
    line="$(parse_reads_tsv "$rtsv" || true)"
    if [[ -n "$line" ]]; then
      # fields: s ma ea ta mc ec tc
      IFS=$'\t' read -r _S MA EA TA MC EC TC <<< "$line"
    else
      echo "⚠️  [$SAMPLE] Failed to parse $rtsv, using fallback." >&2
      IFS=$'\t' read -r MA EA TA MC EC TC <<< "$(load_fallback_metrics "$sample_dir" "$SAMPLE")"
    fi
  else
    echo "⚠️  [$SAMPLE] No .reads.tsv found, using fallback metrics." >&2
    IFS=$'\t' read -r MA EA TA MC EC TC <<< "$(load_fallback_metrics "$sample_dir" "$SAMPLE")"
  fi

  # 3) Derived metrics
  # MouseAlignRate (%) = MouseAligned / TrimmedReadPairs * 100
  if (( TRIM_PAIRS > 0 )); then
    MOUSE_ALIGN_RATE=$(awk -v a="${MA:-0}" -v b="${TRIM_PAIRS}" 'BEGIN{printf "%.2f", (b==0)?0:(a/b*100)}')
  else
    MOUSE_ALIGN_RATE="NA"
  fi

  # PortionOfEcoli = EcoliClean / TotalClean
  if [[ -n "${TC:-}" && "${TC:-0}" -gt 0 ]]; then
    PORTION_ECOLI=$(awk -v ec="${EC:-0}" -v tc="${TC:-0}" 'BEGIN{printf "%.6f", (tc==0)?0:(ec/tc)}')
  else
    PORTION_ECOLI="NA"
  fi

  # 4) Emit row
  printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" \
    "$PROJECT" "$SAMPLE" \
    "$RAW_READS" "$TRIM_READS" \
    "${MA:-0}" "${EA:-0}" "${TA:-0}" \
    "${MC:-0}" "${EC:-0}" "${TC:-0}" \
    "$MOUSE_ALIGN_RATE" "$PORTION_ECOLI" >> "$OUT_FILE"

  echo "✅ $SAMPLE done."
  ((processed++))
done
shopt -u nullglob

echo
echo "Summary:"
echo "  Samples processed: $processed"
echo "  Output: $OUT_FILE"
echo "  Log   : $LOG_FILE"
