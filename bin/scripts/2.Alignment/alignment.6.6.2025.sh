#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --time=48:00:00
# ==========================================
# Align to E. coli, then align host reads to Mouse Genome
# ==========================================
# Usage: bash align_ecoli_then_mouse.sh <R1.fastq.gz> <R2.fastq.gz> <ecoli_index> <mouse_index> <output_dir>

set -euo pipefail

# ------------------- keep your original argument order -------------------
R1=$1; R2=$2; ECOLI_INDEX=$3; MOUSE_INDEX=$4; OUT=$5
[[ -z "$R1" || -z "$R2" || -z "$ECOLI_INDEX" || -z "$MOUSE_INDEX" || -z "$OUT" ]] && {
  echo "Usage: bash align_ecoli_then_mouse.sh <R1> <R2> <ecoli_index> <mouse_index> <outdir>"; exit 1; }

SAMPLE=$(basename "$R1" | sed 's/_R1.*//')
CPU=${SLURM_CPUS_PER_TASK:-4}
echo $CPU
# output folders (unchanged)
ECOLI_OUT="${OUT}/ecoli_alignment"
MOUSE_OUT="${OUT}/mouse_alignment"
mkdir -p "$ECOLI_OUT" "$MOUSE_OUT"

module load bowtie2 samtools picard

# =========================== STEP 1 – E. coli =============================
echo "🧬 Aligning $SAMPLE to E. coli..."
bowtie2 -x "$ECOLI_INDEX" -1 "$R1" -2 "$R2" \
        --very-sensitive --no-mixed --no-discordant \
        --threads $CPU \
        --un-conc-gz "${ECOLI_OUT}/${SAMPLE}_host_unaligned_R%.fastq.gz" \
  | samtools sort -@ $CPU -o "${ECOLI_OUT}/${SAMPLE}_ecoli_sorted.bam" \
  2> "${ECOLI_OUT}/${SAMPLE}_ecoli_alignment.log"

# keep uniquely mapped spike-in reads (MAPQ ≥ 10)
samtools view -h -b -q 10 "${ECOLI_OUT}/${SAMPLE}_ecoli_sorted.bam" \
  > "${ECOLI_OUT}/${SAMPLE}_ecoli_unique.bam"

ecoli_reads=$(samtools view -c -F 4 "${ECOLI_OUT}/${SAMPLE}_ecoli_unique.bam")
echo -e "${SAMPLE}\t${ecoli_reads}" > "${ECOLI_OUT}/${SAMPLE}_ecoli_unique_counts.txt"
samtools flagstat "${ECOLI_OUT}/${SAMPLE}_ecoli_unique.bam" \
  > "${ECOLI_OUT}/${SAMPLE}_ecoli_unique.flagstat.txt"

# =========================== STEP 2 – Mouse ===============================
echo "🐭 Aligning host (mouse) reads for $SAMPLE..."

UN_R1="${ECOLI_OUT}/${SAMPLE}_host_unaligned_R1.fastq.gz"
UN_R2="${ECOLI_OUT}/${SAMPLE}_host_unaligned_R2.fastq.gz"

bowtie2 -x "$MOUSE_INDEX" -1 "$UN_R1" -2 "$UN_R2" \
        --local --very-sensitive --no-mixed --no-discordant \
        -I 10 -X 2000 \
        --threads $CPU \
        --un-conc-gz "${MOUSE_OUT}/${SAMPLE}_final_unaligned_R%.fastq.gz" \
        --rg-id "$SAMPLE" --rg SM:"$SAMPLE" \
  | samtools sort -@ $CPU -o "${MOUSE_OUT}/${SAMPLE}_mouse_sorted.bam" \
  2> "${MOUSE_OUT}/${SAMPLE}_mouse_alignment.log"

samtools addreplacerg \
  -w \
  -r ID:${SAMPLE} -r LB:lib1 -r PL:ILLUMINA -r PU:unit1 -r SM:${SAMPLE} \
  -o "${MOUSE_OUT}/${SAMPLE}_mouse_rg.bam" \
  "${MOUSE_OUT}/${SAMPLE}_mouse_sorted.bam"

# mark **and remove** duplicates
picard MarkDuplicates \
  I="${MOUSE_OUT}/${SAMPLE}_mouse_rg.bam" \
  O="${MOUSE_OUT}/${SAMPLE}_mouse_dedup.bam" \
  M="${MOUSE_OUT}/${SAMPLE}_mouse_dup_metrics.txt" \
  REMOVE_DUPLICATES=true \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=SILENT
samtools index "${MOUSE_OUT}/${SAMPLE}_mouse_dedup.bam"

# keep uniquely mapped pairs (MAPQ ≥ 30) for downstream analysis
samtools view -@ $CPU -b -q 30 "${MOUSE_OUT}/${SAMPLE}_mouse_dedup.bam" \
  > "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam"
samtools index "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam"

mouse_reads=$(samtools view -c -F 4 "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam")
echo -e "${SAMPLE}\t${mouse_reads}" > "${MOUSE_OUT}/${SAMPLE}_mouse_unique_counts.txt"
samtools flagstat "${MOUSE_OUT}/${SAMPLE}_mouse_unique.bam" \
  > "${MOUSE_OUT}/${SAMPLE}_mouse_unique.flagstat.txt"

echo "✅ Alignment complete for $SAMPLE"

# =========================== STEP 3 – QC ================================
# echo "📊 Generating QC summary for $SAMPLE..."

# flagstat_file="${MOUSE_OUT}/${SAMPLE}_mouse_unique.flagstat.txt"
# total_reads=$(grep "in total" "$flagstat_file" | awk '{print $1}')
# mapped_reads=$(grep "mapped (" "$flagstat_file" | awk '{print $1}')
# proper_pairs=$(grep "properly paired (" "$flagstat_file" | awk '{print $1}')
# singletons=$(grep "singletons" "$flagstat_file" | awk '{print $1}')
# mate_diff_chr=$(grep "with mate mapped to a different chr$" "$flagstat_file" | awk '{print $1}')
# dup_count=$(grep -v "^#" "${MOUSE_OUT}/${SAMPLE}_mouse_dup_metrics.txt" | awk '{sum+=$3} END{print sum+0}')

# pct_mapped=$(awk -v m=$mapped_reads -v t=$total_reads 'BEGIN{printf "%.2f",(m/t)*100}')
# pct_proper=$(awk -v p=$proper_pairs -v t=$total_reads 'BEGIN{printf "%.2f",(p/t)*100}')
# pct_single=$(awk -v s=$singletons -v t=$total_reads 'BEGIN{printf "%.2f",(s/t)*100}')
# pct_diffchr=$(awk -v d=$mate_diff_chr -v t=$total_reads 'BEGIN{printf "%.2f",(d/t)*100}')

# qc_summary="${MOUSE_OUT}/${SAMPLE}_alignment_qc_summary.csv"
# echo -e "Sample,Total_Reads,Host_Mapped_Reads(%),Properly_Paired(%),Singletons(%),Mate_on_Diff_Chr(%),Duplicates" > "$qc_summary"
# echo -e "${SAMPLE},${total_reads},${pct_mapped},${pct_proper},${pct_single},${pct_diffchr},${dup_count}" >> "$qc_summary"

# echo "🧾 Alignment QC summary written to $qc_summary"