#!/usr/bin/env Rscript

## ------------------------------------------------------------
## Differential CUT&RUN peak calling: KO vs WT
## ------------------------------------------------------------
## INPUT:
##   - Wide peak table like `new peaks.txt` with:
##       * chr, start, end, key
##       * WT*_H3K27me3_SEACR_peaks.stringent (0/1)
##       * KO*_H3K27me3_SEACR_peaks.stringent (0/1)
##       * WT*_H3K27me3_mouse_scaled_{mean,median,sum}
##       * KO*_H3K27me3_mouse_scaled_{mean,median,sum}
##
## METHOD:
##   - Filter peaks by minimum number of WT / KO replicates with a peak
##   - Build WT/KO matrix using chosen metric (mean / median / sum)
##   - log2-transform and run limma (KO vs WT)
##   - Return log2FC, p-value, FDR, direction, significance flag
##
## USAGE (example):
##   Rscript run_diff_CUTnRUN.R \
##       new_peaks.txt \
##       H3K27me3_diff \
##       mean \
##       2 \
##       2 \
##       0.58 \
##       0.05
##
## ARGS:
##   1) input_file      (required)  e.g., new_peaks.txt
##   2) out_prefix      (optional)  default: basename(input_file) w/o extension
##   3) metric          (optional)  one of: mean, median, sum   [default: mean]
##   4) min_WT_calls    (optional)  minimum WT replicates with peak [default: 1]
##   5) min_KO_calls    (optional)  minimum KO replicates with peak [default: 1]
##   6) log2fc_thresh   (optional)  |log2FC| threshold for sig     [default: 1]
##   7) fdr_thresh      (optional)  FDR threshold for sig          [default: 0.05]
##
## OUTPUT:
##   <out_prefix>.differential_all_peaks.tsv
##   <out_prefix>.differential_significant.tsv
##
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("
Usage:
  Rscript run_diff_CUTnRUN.R <input_file> [out_prefix] [metric] [min_WT_calls] [min_KO_calls] [log2fc_thresh] [fdr_thresh]

See header comments for details.\n")
}

## -----------------------------
## USER-EDITABLE DEFAULTS
## -----------------------------
input_file    <- args[1]
out_prefix    <- ifelse(length(args) >= 2 && nzchar(args[2]),
                        args[2],
                        tools::file_path_sans_ext(basename(input_file)))
metric        <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], "mean")
min_WT_calls  <- ifelse(length(args) >= 4 && nzchar(args[4]), as.integer(args[4]), 1L)
min_KO_calls  <- ifelse(length(args) >= 5 && nzchar(args[5]), as.integer(args[5]), 1L)
log2fc_thresh <- ifelse(length(args) >= 6 && nzchar(args[6]), as.numeric(args[6]), 1)
fdr_thresh    <- ifelse(length(args) >= 7 && nzchar(args[7]), as.numeric(args[7]), 0.05)


message("=== Parameters ===")
message("Input file     : ", input_file)
message("Output prefix  : ", out_prefix)
message("Metric         : ", metric, " (mean/median/sum)")
message("min_WT_calls   : ", min_WT_calls)
message("min_KO_calls   : ", min_KO_calls)
message("log2FC thresh  : ", log2fc_thresh)
message("FDR thresh     : ", fdr_thresh)
message("===================")

if (!metric %in% c("mean","median","sum")) {
  stop("Metric must be one of: mean, median, sum")
}

## Choose the metric suffix used in column names
metric_suffix <- switch(
  metric,
  mean   = "mouse_scaled_mean",
  median = "mouse_scaled_median",
  sum    = "mouse_scaled_sum"
)

## ------------------------------------------------------------
## Step 1: Load data
## ------------------------------------------------------------
dt <- fread(input_file)
if (!"key" %in% names(dt)) {
  ## If no 'key' column, create one from chr:start-end
  dt[, key := paste0(chr, ":", start, "-", end)]
}

## ------------------------------------------------------------
## Step 2: Identify peak-call and signal columns
## ------------------------------------------------------------

## SEACR 0/1 peak call columns (for support filters)
wt_call_cols <- grep("^WT[0-9]$", names(dt), value = TRUE)
ko_call_cols <- grep("^KO[0-9]$", names(dt), value = TRUE)

if (length(wt_call_cols) == 0L || length(ko_call_cols) == 0L) {
  stop("Could not find WT/KO SEACR_peaks.stringent columns. Check input column names.")
}

## Metric signal columns
wt_metric_cols <- grep(
  paste0("^WT[0-9]+_H3K27me3_.*", metric_suffix, "$"),
  names(dt),
  value = TRUE
)
ko_metric_cols <- grep(
  paste0("^KO[0-9]+_H3K27me3_.*", metric_suffix, "$"),
  names(dt),
  value = TRUE
)

if (length(wt_metric_cols) == 0L || length(ko_metric_cols) == 0L) {
  stop("Could not find WT/KO metric columns for suffix: ", metric_suffix,
       "\nExample expected patterns: WT1_H3K27me3_", metric_suffix, " etc.")
}

message("Found ", length(wt_metric_cols), " WT metric columns: ",
        paste(wt_metric_cols, collapse = ", "))
message("Found ", length(ko_metric_cols), " KO metric columns: ",
        paste(ko_metric_cols, collapse = ", "))

## ------------------------------------------------------------
## Step 3: Filter peaks by replicate support
## ------------------------------------------------------------

## how many WT/KO replicates have a peak (SEACR call = 1)
dt[, n_WT_calls := rowSums(.SD == 1, na.rm = TRUE), .SDcols = wt_call_cols]
dt[, n_KO_calls := rowSums(.SD == 1, na.rm = TRUE), .SDcols = ko_call_cols]

dt[, keep_for_DE := (n_WT_calls >= min_WT_calls & n_KO_calls >= min_KO_calls)]

message("Total peaks            : ", nrow(dt))
message("Peaks kept for DE test : ", sum(dt$keep_for_DE), " (",
        round(100 * mean(dt$keep_for_DE), 1), "%)")

## ------------------------------------------------------------
## Step 4: Build expression matrix for limma
## ------------------------------------------------------------

## subset only rows to be used in limma
dt_de <- dt[keep_for_DE == TRUE]

if (nrow(dt_de) == 0L) {
  stop("After applying min_WT_calls and min_KO_calls, no peaks remain for DE.")
}

## Extract metric matrix (WT + KO)
expr_mat <- as.matrix(dt_de[, c(wt_metric_cols, ko_metric_cols), with = FALSE])

## Make sure everything is numeric
if (!is.numeric(expr_mat)) {
  expr_mat <- apply(expr_mat, 2, as.numeric)
  expr_mat <- as.matrix(expr_mat)
}

## Add a small pseudocount to avoid log2(0)
pseudocount <- 1e-3
expr_log2   <- log2(expr_mat + pseudocount)

## Condition design
condition <- c(rep("WT", length(wt_metric_cols)),
               rep("KO", length(ko_metric_cols)))
condition <- factor(condition, levels = c("WT", "KO"))
design <- model.matrix(~ condition)

## ------------------------------------------------------------
## Step 5: Run limma (KO vs WT)
## ------------------------------------------------------------
fit <- lmFit(expr_log2, design)
fit <- eBayes(fit)

## coef for KO vs WT is 'conditionKO'
tt <- topTable(fit,
               coef = "conditionKO",
               number = Inf,
               sort.by = "none")

## tt rows should match dt_de rows in order
if (nrow(tt) != nrow(dt_de)) {
  stop("Mismatch between limma result rows and dt_de rows.")
}

## ------------------------------------------------------------
## Step 6: Add statistics back to full table
## ------------------------------------------------------------

## Calculate group-level means (on original scale) for chosen metric
mean_WT <- rowMeans(dt_de[, ..wt_metric_cols], na.rm = TRUE)
mean_KO <- rowMeans(dt_de[, ..ko_metric_cols], na.rm = TRUE)

dt[keep_for_DE == TRUE,
   c("metric_mean_WT", "metric_mean_KO", "log2FC", "pvalue", "FDR") :=
     .(mean_WT,
       mean_KO,
       tt$logFC,        # KO vs WT (log2)
       tt$P.Value,
       tt$adj.P.Val)]

## Peaks not used for DE will have NA for stats
dt[keep_for_DE == FALSE,
   c("metric_mean_WT", "metric_mean_KO", "log2FC", "pvalue", "FDR") :=
     .(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)]

## ------------------------------------------------------------
## Step 7: Call significance & direction
## ------------------------------------------------------------

dt[, is_sig := (!is.na(FDR) &
                FDR <= fdr_thresh &
                abs(log2FC) >= log2fc_thresh)]

dt[, direction := fifelse(
  is.na(log2FC), NA_character_,
  fifelse(
    is_sig & log2FC > 0, "KO_enriched",
    fifelse(is_sig & log2FC < 0, "WT_enriched", "not_significant")
  )
)]

## Optional: keep the original Group (WT-only, KO-only, Shared) if present
if (!"Group" %in% names(dt)) {
  dt[, Group := NA_character_]
}

## ------------------------------------------------------------
## Step 8: Write outputs
## ------------------------------------------------------------
out_all <- paste0(out_prefix, ".differential_all_peaks.tsv")
out_sig <- paste0(out_prefix, ".differential_significant.tsv")

fwrite(dt, out_all, sep = "\t", quote = FALSE)

dt_sig <- dt[is_sig == TRUE]
fwrite(dt_sig, out_sig, sep = "\t", quote = FALSE)

message("Wrote all peaks with statistics   : ", out_all)
message("Wrote significant differential peaks: ", out_sig)
message("  (n = ", nrow(dt_sig), " peaks with |log2FC| >= ",
        log2fc_thresh, " and FDR <= ", fdr_thresh, ")")

message("Done.")
