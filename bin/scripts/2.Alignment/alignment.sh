#!/bin/bash
#SBATCH --time=72:00:00
#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --time=48:00:00
# ==========================================
# Align to Mouse (primary), then to E. coli (spike-in)
# ==========================================
# Usage: bash align_mouse_then_ecoli.sh <R1.fastq.gz> <R2.fastq.gz> \
#       <mouse_index> <ecoli_index> <mouse_blacklist.bed> <output_dir>


# ------------------- arguments -------------------
R1=$1; R2=$2
MOUSE_INDEX=$3
ECOLI_INDEX=$4
OUT=$5
BLACKLIST=$6

[[ -z "$R1" || -z "$R2" || -z "$MOUSE_INDEX" || -z "$ECOLI_INDEX" || -z "$BLACKLIST" || -z "$OUT" ]] && {
  echo "Usage: bash align_mouse_then_ecoli.sh <R1> <R2> <mouse_index> <ecoli_index> <mouse_blacklist.bed> <outdir>"
  exit 1
}

SAMPLE=$(basename "$R1" | sed 's/_R1.*//')
CPU=${SLURM_CPUS_PER_TASK:-1}
echo "🔍 Aligning sample: $SAMPLE with $CPU threads"
# output folders
MOUSE_OUT="${OUT}/mouse_alignment"
ECOLI_OUT="${OUT}/ecoli_alignment"
mkdir -p "$MOUSE_OUT" "$ECOLI_OUT"


# --- LOGGING ---
LOGFILE="${OUT}/alignment.log"
echo "📄 Logging run to $LOGFILE"
{
  echo "Date: $(date)"
  echo "R1: $R1"
  echo "R2: $R2"
  echo "Blacklist: $BLACKLIST"
  echo "Output Dir: $OUT"
} > "$LOGFILE"

cat ${LOGFILE}



module load bowtie2 samtools picard bedtools

# =========================== STEP 1 – Mouse ================================
echo "🐭 Aligning $SAMPLE to mouse genome..."

bowtie2 -x "$MOUSE_INDEX" \
        -1 "$R1" -2 "$R2" \
        --local --very-sensitive --no-mixed --no-discordant -I 10 -X 2000 \
        --threads $CPU \
        --un-conc-gz "${MOUSE_OUT}/${SAMPLE}_unaligned_R%.fastq.gz" \
        --rg-id "$SAMPLE" --rg SM:"$SAMPLE" \
  | samtools sort -@ $CPU -o "${MOUSE_OUT}/${SAMPLE}_mouse_sorted.bam" \
  2> "${MOUSE_OUT}/${SAMPLE}_mouse_alignment.log"

samtools flagstat "${MOUSE_OUT}/${SAMPLE}_mouse_sorted.bam" > "${MOUSE_OUT}/${SAMPLE}_mouse_total.flagstat.txt"

# 1a) keep only uniquely‐mapped pairs (MAPQ ≥ 30)
samtools view -@ $CPU -b -q 30 "${MOUSE_OUT}/${SAMPLE}_mouse_sorted.bam" \
  > "${MOUSE_OUT}/${SAMPLE}_mouse_mapq30.bam"

# 1a.5) strip /1 and /2 from read names before blacklist filtering
samtools view -h "${MOUSE_OUT}/${SAMPLE}_mouse_mapq30.bam" \
  | sed 's/\/[12]$//' \
  | samtools view -b -@ $CPU -o "${MOUSE_OUT}/${SAMPLE}_mouse_mapq30_clean.bam"

# 1b) remove ENCODE DAC exclusion (blacklist) regions
bedtools intersect -v \
  -a "${MOUSE_OUT}/${SAMPLE}_mouse_mapq30_clean.bam" \
  -b "$BLACKLIST" \
  > "${MOUSE_OUT}/${SAMPLE}_mouse_filtered.bam"

# 1c) add read groups (if not already)
samtools addreplacerg  -w \
  -@ $CPU -r ID:${SAMPLE} -r LB:lib1 -r PL:ILLUMINA -r PU:unit1 -r SM:${SAMPLE} \
  -o "${MOUSE_OUT}/${SAMPLE}_mouse_rg.bam" \
  "${MOUSE_OUT}/${SAMPLE}_mouse_filtered.bam"


# 1d) mark (but do not remove) duplicates
picard MarkDuplicates \
  I="${MOUSE_OUT}/${SAMPLE}_mouse_rg.bam" \
  O="${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam" \
  M="${MOUSE_OUT}/${SAMPLE}_mouse_dup_metrics.txt" \
  REMOVE_DUPLICATES=false \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT

samtools index "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam"

# 1e) final count & stats
mouse_reads=$(samtools view -c -F 4 "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam")
echo -e "${SAMPLE}\t${mouse_reads}" > "${MOUSE_OUT}/${SAMPLE}_mouse_counts.txt"
samtools flagstat "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam" \
  > "${MOUSE_OUT}/${SAMPLE}_mouse.flagstat.txt"



mouse_reads=$(samtools view -c -F 4 "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam")
echo -e "${SAMPLE}\t${mouse_reads}" > "${MOUSE_OUT}/${SAMPLE}_mouse_unique_counts.txt"
samtools flagstat "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam" \
  > "${MOUSE_OUT}/${SAMPLE}_mouse_unique.flagstat.txt"




# =========================== STEP 2 – E. coli ==============================
echo "🧬 Aligning host‐unaligned reads to E. coli..."

UN_R1="${MOUSE_OUT}/${SAMPLE}_unaligned_R1.fastq.gz"
UN_R2="${MOUSE_OUT}/${SAMPLE}_unaligned_R2.fastq.gz"

bowtie2 -x "$ECOLI_INDEX" \
        -1 "$UN_R1" -2 "$UN_R2" \
        --end-to-end --very-sensitive --no-mixed --no-discordant -I 10 -X 700 \
        --threads $CPU \
        --un-conc-gz "${ECOLI_OUT}/${SAMPLE}_ecoli_unaligned_R%.fastq.gz" \
  | samtools sort -@ $CPU -o "${ECOLI_OUT}/${SAMPLE}_ecoli_sorted.bam" \
  2> "${ECOLI_OUT}/${SAMPLE}_ecoli_alignment.log"

# keep only uniquely‐mapped reads (MAPQ ≥ 10)
samtools view -@ $CPU -b -q 30 "${ECOLI_OUT}/${SAMPLE}_ecoli_sorted.bam" \
  > "${ECOLI_OUT}/${SAMPLE}_ecoli_unique.bam"

# counts & stats
ecoli_reads=$(samtools view -c -F 4 "${ECOLI_OUT}/${SAMPLE}_ecoli_unique.bam")
echo -e "${SAMPLE}\t${ecoli_reads}" > "${ECOLI_OUT}/${SAMPLE}_ecoli_counts.txt"
samtools flagstat "${ECOLI_OUT}/${SAMPLE}_ecoli_unique.bam" \
  > "${ECOLI_OUT}/${SAMPLE}_ecoli.flagstat.txt"

echo "✅ Alignment complete for $SAMPLE"
