#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --time=48:00:00

# ==========================================
# Single-pass alignment to combined mouse+E.coli index
# Split to mouse-only & ecoli-only BAMs
# Write organism counts (pre-filter & post-filter)
# Clean up intermediates at the end
# Usage:
#   bash align_mouse_ecoli.sh R1.fastq.gz R2.fastq.gz COMB_INDEX OUT_DIR BLACKLIST.bed
# ==========================================

# -------- arguments --------
R1=${1:?R1 fastq.gz}
R2=${2:?R2 fastq.gz}
COMB_INDEX=${3:?Combined Bowtie2 index prefix}
OUT=${4:?Output directory}
BLACKLIST=${5:?Mouse blacklist BED}

SAMPLE=$(basename "$R1" | sed 's/_R1.*//')
CPU=${SLURM_CPUS_PER_TASK:-8}

# -------- toggles (set to "true"/"false") --------
KEEP_COMBINED=false       # keep ${SAMPLE}.markdup.bam ?
KEEP_ECOLI_BAM=true       # keep ${SAMPLE}.ecoli.bam for QC?
KEEP_CHR_LISTS=false      # keep mouse.chroms.txt / ecoli.chroms.txt ?

# -------- setup / logging --------
mkdir -p "$OUT"
LOGFILE="${OUT}/alignment.log"

{
  echo "=== CUT&RUN combined align ===  $(date)"
  echo "SAMPLE: $SAMPLE"
  echo "R1: $R1"
  echo "R2: $R2"
  echo "COMB_INDEX: $COMB_INDEX"
  echo "BLACKLIST: $BLACKLIST"
  echo "OUT: $OUT"
  echo "CPU: $CPU"
} > "$LOGFILE"

module load bowtie2 samtools bedtools picard

echo "🐭🧬 Aligning ${SAMPLE} to combined (mouse+ecoli)..."
bowtie2 -x "$COMB_INDEX" \
  -1 "$R1" -2 "$R2" \
  --very-sensitive-local --no-mixed --no-discordant -I 10 -X 2000 \
  --threads "$CPU" \
  --rg-id "$SAMPLE" --rg SM:"$SAMPLE" \
| samtools sort -@ "$CPU" -o "${OUT}/${SAMPLE}.sorted.bam" \
  2>> "$LOGFILE"
samtools index -@ "$CPU" "${OUT}/${SAMPLE}.sorted.bam"

# ---------------- pre-filter organism counts (from sorted.bam) ----------------
# Sum mapped reads by organism using contig prefix 'ecoli|'
read MOUSE_PRE ECOLI_PRE TOTAL_PRE < <(
  samtools idxstats "${OUT}/${SAMPLE}.sorted.bam" \
  | awk '
      $1=="*" {next}
      $1 ~ /^ecoli\|/ {ec+=$3; next}
      {mm+=$3}
      END{printf "%d %d %d", mm+0, ec+0, mm+ec+0}
    '
)

# ---------------- filter (same rules for both organisms) ----------------
# Properly paired, MAPQ>=30; drop secondary(256),supp(2048),QC-fail(512) → 2816; KEEP duplicates
samtools view -@ "$CPU" -b -q 30 -f 2 -F 2816 \
  "${OUT}/${SAMPLE}.sorted.bam" > "${OUT}/${SAMPLE}.flt.bam"

# Remove blacklist regions (preserve BAM header)
bedtools intersect -v -ubam -a "${OUT}/${SAMPLE}.flt.bam" -b "$BLACKLIST" \
  > "${OUT}/${SAMPLE}.flt.bl.bam"

# Mark (do not remove) duplicates
picard MarkDuplicates \
  I="${OUT}/${SAMPLE}.flt.bl.bam" \
  O="${OUT}/${SAMPLE}.markdup.bam" \
  M="${OUT}/${SAMPLE}.dup.txt" \
  REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT \
  >> "$LOGFILE" 2>&1
samtools index -@ "$CPU" "${OUT}/${SAMPLE}.markdup.bam"

FINAL="${OUT}/${SAMPLE}.markdup.bam"

# ---------------- split by organism ----------------
samtools idxstats "$FINAL" | awk '$1!="*" && $1 !~ /^ecoli\|/ {print $1}' > "${OUT}/mouse.chroms.txt"
samtools idxstats "$FINAL" | awk '$1 ~  /^ecoli\|/                {print $1}' > "${OUT}/ecoli.chroms.txt"

# Mouse-only BAM
samtools view -@ "$CPU" -b "$FINAL" $(cat "${OUT}/mouse.chroms.txt") > "${OUT}/${SAMPLE}.mouse.bam"
samtools index -@ "$CPU" "${OUT}/${SAMPLE}.mouse.bam"

# E. coli-only BAM (optional)
samtools view -@ "$CPU" -b "$FINAL" $(cat "${OUT}/ecoli.chroms.txt") > "${OUT}/${SAMPLE}.ecoli.bam"
samtools index -@ "$CPU" "${OUT}/${SAMPLE}.ecoli.bam"

# ---------------- post-filter organism counts (from split BAMs) ----------------
MOUSE_POST=$(samtools idxstats "${OUT}/${SAMPLE}.mouse.bam" | awk '$1!="*"{n+=$3} END{print n+0}')
ECOLI_POST=$(samtools idxstats "${OUT}/${SAMPLE}.ecoli.bam" | awk '$1!="*"{n+=$3} END{print n+0}')
TOTAL_POST=$(( MOUSE_POST + ECOLI_POST ))


# ---------------- single per-sample combined file ----------------
SREAD="${OUT}/${SAMPLE}.reads.tsv"
echo -e "Sample\tMouseAligned\tEcoliAligned\tTotalAligned\tMouseClean\tEcoliClean\tTotalClean" > "$SREAD"
echo -e "${SAMPLE}\t${MOUSE_PRE}\t${ECOLI_PRE}\t${TOTAL_PRE}\t${MOUSE_POST}\t${ECOLI_POST}\t${TOTAL_POST}" >> "$SREAD"


echo -e "${SAMPLE}\t${MOUSE_POST}\t${ECOLI_POST}\t${TOTAL_POST}" >> "${OUT}/combined_counts.postfilter.tsv"

echo "✅ Alignment, split, and counting complete for $SAMPLE" | tee -a "$LOGFILE"

# ---------------- cleanup ----------------
echo "🧹 Cleaning up intermediates..." | tee -a "$LOGFILE"

rm -f "${OUT}/${SAMPLE}.sorted.bam" "${OUT}/${SAMPLE}.sorted.bam.bai" || true
rm -f "${OUT}/${SAMPLE}.flt.bam" || true
rm -f "${OUT}/${SAMPLE}.flt.bl.bam" || true

if [[ "${KEEP_COMBINED}" != "true" ]]; then
  rm -f "${OUT}/${SAMPLE}.markdup.bam" "${OUT}/${SAMPLE}.markdup.bam.bai" || true
fi

if [[ "${KEEP_ECOLI_BAM}" != "true" ]]; then
  rm -f "${OUT}/${SAMPLE}.ecoli.bam" "${OUT}/${SAMPLE}.ecoli.bam.bai" || true
fi

if [[ "${KEEP_CHR_LISTS}" != "true" ]]; then
  rm -f "${OUT}/mouse.chroms.txt" "${OUT}/ecoli.chroms.txt" || true
fi

echo "Done. $(date)" >> "$LOGFILE"
