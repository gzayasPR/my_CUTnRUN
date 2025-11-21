#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10gb
#SBATCH --account=mateescu

# Usage: compare_seacr_peaks.sh <WT_peak_file> <KO_peak_file> <out_dir>


# --- Inputs ---
WT1="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/4.PeakCalling/WT1_H3K27me3/SEACR/WT1_H3K27me3_SEACR_peaks.stringent.bed"
WT2="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/4.PeakCalling/WT2_H3K27me3/SEACR/WT2_H3K27me3_SEACR_peaks.stringent.bed"
WT3="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/4.PeakCalling/WT3_H3K27me3/SEACR/WT3_H3K27me3_SEACR_peaks.stringent.bed"
KO1="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/4.PeakCalling/KO1_H3K27me3/SEACR/KO1_H3K27me3_SEACR_peaks.stringent.bed"
KO2="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/4.PeakCalling/KO2_H3K27me3/SEACR/KO2_H3K27me3_SEACR_peaks.stringent.bed"
KO3="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/4.PeakCalling/KO3_H3K27me3/SEACR/KO3_H3K27me3_SEACR_peaks.stringent.bed"
OUT_DIR="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/5.PeakComparison/WT_vs_KO_multiple"
proj_env="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/bin/project_env.sh"
source "$proj_env"
MERGE_DIST=50  # if >0, we’ll stitch close intervals after building union
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merge-dist) MERGE_DIST="$2"; shift 2;;
    *) echo "Unknown option: $1"; exit 1;;
  esac
done

mkdir -p "$OUT_DIR"

# --- Modules (adapt to your cluster) ---
ml bedtools || module load bedtools || true

# --- Names for readability in outputs (derived from basenames) ---
WT1_NAME=$(basename "${WT1%.*}")
WT2_NAME=$(basename "${WT2%.*}")
WT3_NAME=$(basename "${WT3%.*}")
KO1_NAME=$(basename "${KO1%.*}")
KO2_NAME=$(basename "${KO2%.*}")
KO3_NAME=$(basename "${KO3%.*}")

# --- Logging ---
LOGFILE="${OUT_DIR}/compare_seacr_peaks_union_matrix_run.log"
{
  echo "Date: $(date)"
  echo "WT files:"
  echo "  1) $WT1  ($WT1_NAME)"
  echo "  2) $WT2  ($WT2_NAME)"
  echo "  3) $WT3  ($WT3_NAME)"
  echo "KO files:"
  echo "  1) $KO1  ($KO1_NAME)"
  echo "  2) $KO2  ($KO2_NAME)"
  echo "  3) $KO3  ($KO3_NAME)"
  echo "Output Directory: $OUT_DIR"
  echo "Parameters:"
  echo "  MERGE_DIST = $MERGE_DIST"
  echo ""
} > "$LOGFILE"
echo "📄 Logging run parameters to $LOGFILE"
cat "$LOGFILE"

sort_bed() { bedtools sort -i "$1"; }

UNION_RAW="${OUT_DIR}/multiinter_raw.tsv"
UNION_TSV="${OUT_DIR}/peaks_union_presence.tsv"

echo "sample,bed
WT1,$WT1
WT2,$WT2
WT3,$WT3
KO1,$KO1
KO2,$KO2
KO3,$KO3
" > "${OUT_DIR}/input.csv"
# --- Step 1: Multi-intersection to create union of all six ---
# bedtools multiinter columns:
# 1:chrom 2:start 3:end 4:num_files_support 5:list_of_file_indexes
echo "🔧 Computing union across all six replicates with R enhanced script similar to bedtools multiinter..."
# bedtools multiinter \
#   -i <(sort_bed "$WT1") <(sort_bed "$WT2") <(sort_bed "$WT3") \
#      <(sort_bed "$KO1") <(sort_bed "$KO2") <(sort_bed "$KO3") \
#   > "$UNION_RAW"
ml R
Rscript ${my_bin}/scripts/5.Annotation/merge_peaks_multi.R \
        "${OUT_DIR}/input.csv" \
        "$UNION_TSV"


# # --- Step 2: Build presence/absence matrix + group labels ---
# # File index mapping from multiinter to sample columns:
# # 1->WT1, 2->WT2, 3->WT3, 4->KO1, 5->KO2, 6->KO3
# echo "🧮 Annotating per-sample presence matrix..."
# {
#   printf "chrom\tstart\tend\tn_support\tfiles\t%s\t%s\t%s\t%s\t%s\t%s\tWT_any\tKO_any\tGroup\n" \
#     "$WT1_NAME" "$WT2_NAME" "$WT3_NAME" "$KO1_NAME" "$KO2_NAME" "$KO3_NAME"

#   awk -v OFS="\t" \
#       -v w1="$WT1_NAME" -v w2="$WT2_NAME" -v w3="$WT3_NAME" \
#       -v k1="$KO1_NAME" -v k2="$KO2_NAME" -v k3="$KO3_NAME" \
#   '
#   BEGIN{
#     # nothing
#   }
#   NR>=1{
#     chrom=$1; start=$2; end=$3; nsupp=$4; list=$5
#     # initialize flags
#     WT1=WT2=WT3=KO1=KO2=KO3=0

#     # list is comma-separated file indexes (1..6)
#     n=split(list, arr, ",")
#     for(i=1;i<=n;i++){
#       if(arr[i]==1){WT1=1}
#       else if(arr[i]==2){WT2=1}
#       else if(arr[i]==3){WT3=1}
#       else if(arr[i]==4){KO1=1}
#       else if(arr[i]==5){KO2=1}
#       else if(arr[i]==6){KO3=1}
#     }

#     WT_any = (WT1||WT2||WT3)?1:0
#     KO_any = (KO1||KO2||KO3)?1:0

#     group="NA"
#     if (WT_any==1 && KO_any==1) {group="Shared"}
#     else if (WT_any==1 && KO_any==0) {group="WT-only"}
#     else if (WT_any==0 && KO_any==1) {group="KO-only"}

#     print chrom, start, end, nsupp, list, WT1, WT2, WT3, KO1, KO2, KO3, WT_any, KO_any, group
#   }' "$UNION_RAW"
# } > "$UNION_TSV"

# echo "✅ Presence matrix written: $UNION_TSV"

# --- Step 3: Export BEDs for Shared / WT-only / KO-only ---
SHARED_BED="${OUT_DIR}/shared_peaks.bed"
WT_ONLY_BED="${OUT_DIR}/WT_only.bed"
KO_ONLY_BED="${OUT_DIR}/KO_only.bed"

awk -v OFS="\t" 'NR>1 && $14=="Shared"{print $1,$2,$3}' "$UNION_TSV" \
  | bedtools sort -i - > "$SHARED_BED"
awk -v OFS="\t" 'NR>1 && $14=="WT-only"{print $1,$2,$3}' "$UNION_TSV" \
  | bedtools sort -i - > "$WT_ONLY_BED"
awk -v OFS="\t" 'NR>1 && $14=="KO-only"{print $1,$2,$3}' "$UNION_TSV" \
  | bedtools sort -i - > "$KO_ONLY_BED"

# --- Step 4: Summaries ---
SH_N=$(wc -l < "$SHARED_BED" || echo 0)
WO_N=$(wc -l < "$WT_ONLY_BED" || echo 0)
KO_N=$(wc -l < "$KO_ONLY_BED" || echo 0)
TOTAL_N=$(($(wc -l < "$UNION_TSV") - 1))

{
  echo ""
  echo "Summary:"
  echo "  Total union intervals: $TOTAL_N"
  echo "  Shared:                $SH_N"
  echo "  WT-only:               $WO_N"
  echo "  KO-only:               $KO_N"
  echo ""
  echo "Outputs:"
  echo "  Matrix TSV:   $UNION_TSV"
  echo "  Shared BED:   $SHARED_BED"
  echo "  WT-only BED:  $WT_ONLY_BED"
  echo "  KO-only BED:  $KO_ONLY_BED"
} | tee -a "$LOGFILE"

echo "🎉 Done."

###############################################################################
# Extend: annotate union peaks with per-sample signal stats (bigWig -> bedGraph)
# Assumptions:
#   1) You already ran the union/matrix steps and have $UNION_TSV, $OUT_DIR set.
#   2) You have SIX bigWig signal files exported as environment variables:
#        SIGNAL_WT1, SIGNAL_WT2, SIGNAL_WT3, SIGNAL_KO1, SIGNAL_KO2, SIGNAL_KO3
#      e.g.,
#        export SIGNAL_WT1=/path/WT1.bw
#        export SIGNAL_WT2=/path/WT2.bw
#        ...
#   3) UCSC bigWigToBedGraph is available in PATH, and bedtools is loaded.
###############################################################################

echo "➕ Starting signal annotation (bigWig → bedGraph → stats)..."

SIGNAL_IgG_KO="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO_IgG/KO_IgG_mouse_scaled.bw"
SIGNAL_KO1="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO1_H3K27me3/KO1_H3K27me3_mouse_scaled.bw"
SIGNAL_KO2="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO2_H3K27me3/KO2_H3K27me3_mouse_scaled.bw"
SIGNAL_KO3="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/KO3_H3K27me3/KO3_H3K27me3_mouse_scaled.bw"

SIGNAL_IgG_WT="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT_IgG/WT_IgG_mouse_scaled.bw"
SIGNAL_WT1="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT1_H3K27me3/WT1_H3K27me3_mouse_scaled.bw"
SIGNAL_WT2="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT2_H3K27me3/WT2_H3K27me3_mouse_scaled.bw"
SIGNAL_WT3="/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project6_H3K27me3/3.Normalization/WT3_H3K27me3/WT3_H3K27me3_mouse_scaled.bw"


# --- Validate required env vars for signals ---
REQS=(SIGNAL_WT1 SIGNAL_WT2 SIGNAL_WT3 SIGNAL_KO1 SIGNAL_KO2 SIGNAL_KO3)
for v in "${REQS[@]}"; do
  if [[ -z "${!v:-}" ]]; then
    echo "ERROR: Environment variable $v is not set. Please export paths to your .bw files." >&2
  fi
done
ml deeptools
ml conda
conda activate "${my_softwares}/conda_envs/bigwig"

# --- Helper functions ---
sort_bed_like() { sort -k1,1 -k2,2n "$1"; }

name_from() {
  local p="$1"; p=$(basename "$p")
  # strip double extensions like .bw or .bigWig
  p="${p%.bigWig}"; p="${p%.bw}"
  echo "$p"
}

bigwig_to_bedgraph() {
  local in="$1" out="$2"
  if [[ ! -s "$out" ]]; then
    echo "🔄 Converting bigWig → bedGraph: $(basename "$in")"
    bigWigToBedGraph "$in" "$out"
  fi
}



compute_stats() {
  # compute_stats <sample_name> <bedgraph_in> -> writes ${OUT_DIR}/${sample_name}.stats.tsv
  local sample="$1" bg_in="$2"
  local bg_sorted="${OUT_DIR}/${sample}.sorted.bedgraph"
  local out_tab="${OUT_DIR}/${sample}.stats.tsv"

  sort_bed_like "$bg_in" > "$bg_sorted"

  # Map signal (col 4) onto union peaks and compute mean, median, sum, max (fill 0 if no overlap)
  bedtools map -a "$UNION_BED" -b "$bg_sorted" -c 4 -o mean,median,sum,max -null 0 > "${out_tab}.raw"

  # Attach key and sanitize '.' to 0
  paste <(cat "$PEAK_KEY") "${out_tab}.raw" \
    | awk -v OFS="\t" '{
        mean=$5; median=$6; sum=$7; max=$8;
        if(mean==".") mean=0; if(median==".") median=0; if(sum==".") sum=0; if(max==".") max=0;
        print $1, mean, median, sum, max
      }' > "$out_tab"

  rm -f "$bg_sorted" "${out_tab}.raw"
}

# --- Derive union peaks BED from $UNION_TSV (already produced earlier) ---
UNION_BED="${OUT_DIR}/peaks_union_nonoverlap.bed"
awk -v OFS="\t" 'NR>1{print $1,$2,$3}' "$UNION_TSV" > "${UNION_BED}.unsorted"
sort_bed_like "${UNION_BED}.unsorted" > "$UNION_BED"
rm -f "${UNION_BED}.unsorted"

# Key to join on (chrom:start-end)
PEAK_KEY="${OUT_DIR}/peaks_union.key"
awk -v OFS="\t" '{print $1 ":" $2 "-" $3}' "$UNION_BED" > "$PEAK_KEY"

# --- Convert six bigWigs to bedGraph ---
WT1_NAME=$(name_from "$SIGNAL_WT1"); WT1_BG="${OUT_DIR}/${WT1_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_WT1" "$WT1_BG"
WT2_NAME=$(name_from "$SIGNAL_WT2"); WT2_BG="${OUT_DIR}/${WT2_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_WT2" "$WT2_BG"
WT3_NAME=$(name_from "$SIGNAL_WT3"); WT3_BG="${OUT_DIR}/${WT3_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_WT3" "$WT3_BG"
KO1_NAME=$(name_from "$SIGNAL_KO1"); KO1_BG="${OUT_DIR}/${KO1_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_KO1" "$KO1_BG"
KO2_NAME=$(name_from "$SIGNAL_KO2"); KO2_BG="${OUT_DIR}/${KO2_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_KO2" "$KO2_BG"
KO3_NAME=$(name_from "$SIGNAL_KO3"); KO3_BG="${OUT_DIR}/${KO3_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_KO3" "$KO3_BG"
IgG_WT_NAME=$(name_from "$SIGNAL_IgG_WT"); IgG_WT_BG="${OUT_DIR}/${IgG_WT_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_IgG_WT" "$IgG_WT_BG"
IgG_KO_NAME=$(name_from "$SIGNAL_IgG_KO"); IgG_KO_BG="${OUT_DIR}/${IgG_KO_NAME}.bedgraph"; bigwig_to_bedgraph "$SIGNAL_IgG_KO" "$IgG_KO_BG"

# --- Compute stats per sample ---
compute_stats "$WT1_NAME" "$WT1_BG"
compute_stats "$WT2_NAME" "$WT2_BG"
compute_stats "$WT3_NAME" "$WT3_BG"
compute_stats "$KO1_NAME" "$KO1_BG"
compute_stats "$KO2_NAME" "$KO2_BG"
compute_stats "$KO3_NAME" "$KO3_BG"

compute_stats "$IgG_WT_NAME" "$IgG_WT_BG"
compute_stats "$IgG_KO_NAME" "$IgG_KO_BG"

# --- Build wide table: original presence matrix + per-sample stats ---
UNION_WITH_KEY="${OUT_DIR}/peaks_union_presence.withkey.tsv"
{
  read -r header < <(head -n1 "$UNION_TSV")
  echo -e "key\t${header}"
  awk -v OFS="\t" 'NR>1{print $1 ":" $2 "-" $3, $0}' "$UNION_TSV"
} > "$UNION_WITH_KEY"


# Sanitize helper (nukes CRLF if any)
sanitize_tsv() { sed -e 's/\r$//' "$1" > "$1.san" && mv "$1.san" "$1"; }
add_stats_cols() {
  # add_stats_cols <main_in> <stats_tsv> <sample_name> <out_tsv>
  local main_in="$1"
  local stats="$2"
  local sample="$3"
  local out="$4"

  sanitize_tsv "$main_in"
  sanitize_tsv "$stats"

  # header
  local hdr
  hdr=$(head -n1 "$main_in")
  printf "%s\t%s_mean\t%s_median\t%s_sum\t%s_max\n" "$hdr" "$sample" "$sample" "$sample" "$sample" > "$out"

  # left-join using awk map; fill zeros when missing
  awk -v OFS="\t" '
    FNR==1 && NR==1 { next }               # skip header of stats if it’s first file
    NR==FNR && FNR>1 {                     # pass 1: stats file -> map
      m[$1]=$2"\t"$3"\t"$4"\t"$5; next
    }
    NR>FNR {
      if (FNR==1) { next }                 # skip header of main
      key=$1
      vals = (key in m) ? m[key] : "0\t0\t0\t0"
      print $0, vals
    }
  ' "$stats" "$main_in" >> "$out"
}


TMP1="${OUT_DIR}/step1.tmp.tsv"
TMP2="${OUT_DIR}/step2.tmp.tsv"
TMP3="${OUT_DIR}/step3.tmp.tsv"
TMP4="${OUT_DIR}/step4.tmp.tsv"
TMP5="${OUT_DIR}/step5.tmp.tsv"
TMP6="${OUT_DIR}/step6.tmp.tsv"
TMP7="${OUT_DIR}/step7.tmp.tsv"
TMP8="${OUT_DIR}/step8.tmp.tsv"

add_stats_cols "$UNION_WITH_KEY" "${OUT_DIR}/${WT1_NAME}.stats.tsv" "$WT1_NAME" "$TMP1"
add_stats_cols "$TMP1"          "${OUT_DIR}/${WT2_NAME}.stats.tsv" "$WT2_NAME" "$TMP2"
add_stats_cols "$TMP2"          "${OUT_DIR}/${WT3_NAME}.stats.tsv" "$WT3_NAME" "$TMP3"
add_stats_cols "$TMP3"          "${OUT_DIR}/${KO1_NAME}.stats.tsv" "$KO1_NAME" "$TMP4"
add_stats_cols "$TMP4"          "${OUT_DIR}/${KO2_NAME}.stats.tsv" "$KO2_NAME" "$TMP5"
add_stats_cols "$TMP5"          "${OUT_DIR}/${KO3_NAME}.stats.tsv" "$KO3_NAME" "$TMP6"
add_stats_cols "$TMP6"          "${OUT_DIR}/${IgG_WT_NAME}.stats.tsv" "$IgG_WT_NAME" "$TMP7"
add_stats_cols "$TMP7"          "${OUT_DIR}/${IgG_KO_NAME}.stats.tsv" "$IgG_KO_NAME" "$TMP8"

FINAL_STATS="${OUT_DIR}/peaks_union_signal_stats.tsv"
mv "$TMP8" "$FINAL_STATS"
rm -f "$TMP1" "$TMP2" "$TMP3" "$TMP4" "$TMP5" "$TMP5" "$TMP6" "$TMP7" "$UNION_WITH_KEY"

echo "✅ Signal-annotated table: $FINAL_STATS" | tee -a "$LOGFILE"

# Brief summary
TOTAL_ROWS=$(( $(wc -l < "$FINAL_STATS") - 1 ))
{
  echo "Signal annotation summary:"
  echo "  Peaks annotated: $TOTAL_ROWS"
  echo "  Per-sample columns added: *_mean, *_median, *_sum, *_max"
  echo "  File: $FINAL_STATS"
} | tee -a "$LOGFILE"

echo "🏁 Signal annotation complete."
###############################################################################

ANNO_FINAL_STATS="${OUT_DIR}/peaks_union_signal_stats"
MOUSE_GTF="${my_data}/references/Mus_musculus_GRCm39/Ensembl/111/Annotation/Genes/Mus_musculus.GRCm39.111.gtf"

ml R
# |log2FC| >= 0.58 (≈1.5-fold) and FDR <= 0.1
Rscript ${my_bin}/scripts/5.Annotation/run_diff_CUTnRUN.R "$FINAL_STATS" \
                                                            "${OUT_DIR}/peaks_union_signal_log2fc" \
                                                            sum \
                                                            0 \
                                                            0 \
                                                            0.58 \
                                                            0.1


Rscript ${my_bin}/scripts/5.Annotation/annotate_peaks.multiple.R "${OUT_DIR}/peaks_union_signal_log2fc.differential_all_peaks.tsv" \
                                                                 "${OUT_DIR}/annotated_peaks_union_signal_log2fc" \
                                                                 "${MOUSE_GTF}"


