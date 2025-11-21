#!/usr/bin/env Rscript

## ------------------------------------------------------------
## annotate_peaks.R
## - Input: custom peak file (e.g. "new peaks.tsv") with columns:
##   key, chr, start, end, ..., log2FC, pvalue, FDR, is_sig, direction
## - Output:
##   1) All annotated peaks table
##   2) Significant-only annotated peaks table (is_sig == TRUE)
##   3) 3-panel annotation plots:
##        - All peaks
##        - KO_enriched (significant only)
##        - WT_enriched (significant only)
## ------------------------------------------------------------

suppressMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  library(tools)
})

## -----------------------
## Input arguments
## -----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("\nUsage: Rscript annotate_peaks.R <new_peaks.tsv> <out_prefix> <gtf_file>\n",
       " - <new_peaks.tsv>: tab-delimited peak file with columns key, chr, start, end, log2FC, pvalue, FDR, is_sig, direction\n",
       " - <out_prefix>: prefix for output files\n",
       " - <gtf_file>: GTF annotation used to build TxDb\n")
}

peak_file <- args[1]
out_prefix <- args[2]
gtf_file <- args[3]

## -----------------------
## Build / load TxDb
## -----------------------
txdb_file    <- file_path_sans_ext(basename(gtf_file))
txdb_sqlite  <- paste0(txdb_file, "_ensembl.sqlite")

if (!file.exists(txdb_sqlite)) {
  message("🔧 Building TxDb from GTF: ", gtf_file)
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  saveDb(txdb, file = txdb_sqlite)
} else {
  message("📂 Loading TxDb from: ", txdb_sqlite)
  txdb <- loadDb(txdb_sqlite)
}

## -----------------------
## Load custom peak file
## -----------------------
message("📂 Reading peak file: ", peak_file)

peaks_df <- read.table(
  peak_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c("key", "chr", "start", "end", "is_sig", "direction")
missing_cols  <- setdiff(required_cols, colnames(peaks_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in peak file: ",
       paste(missing_cols, collapse = ", "))
}

## Make a 'name' column for GRanges
peaks_df$name <- peaks_df$key

## -----------------------
## Convert to GRanges
## -----------------------
message("🧬 Converting to GRanges...")

peak_gr <- makeGRangesFromDataFrame(
  peaks_df,
  seqnames.field     = "chr",
  start.field        = "start",
  end.field          = "end",
  keep.extra.columns = TRUE,
  ignore.strand      = TRUE
)

## -----------------------
## Annotate ALL peaks
## -----------------------
message("🧬 Annotating ALL peaks...")

peakAnno_all <- annotatePeak(
  peak_gr,
  TxDb      = txdb,
  tssRegion = c(-10000, 10000),
  annoDb    = "org.Mm.eg.db"
)

anno_all_df <- as.data.frame(peakAnno_all)

## Save all annotated peaks
out_all_txt <- paste0(out_prefix, "_annotated_all.tsv")
write.table(
  anno_all_df,
  file      = out_all_txt,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("📄 Saved ALL annotated peaks to: ", out_all_txt)

## -----------------------
## Annotate SIGNIFICANT peaks (is_sig == TRUE)
## -----------------------
message("🎯 Filtering significant peaks (is_sig == TRUE)...")

# is_sig may be logical or character; coerce safely
is_sig_vec <- anno_all_df$is_sig
if (is.logical(is_sig_vec)) {
  sig_idx <- is_sig_vec
} else {
  sig_idx <- tolower(as.character(is_sig_vec)) %in% c("true", "1", "yes")
}

anno_sig_df <- anno_all_df[sig_idx, ]

out_sig_txt <- paste0(out_prefix, "_annotated_significant.tsv")
write.table(
  anno_sig_df,
  file      = out_sig_txt,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("📄 Saved SIGNIFICANT annotated peaks to: ", out_sig_txt)

## -----------------------
## 3-panel plot for ALL peaks (global)
## -----------------------
out_pdf_all <- paste0(out_prefix, "_ALL_annotation_plots_3panel.pdf")
pdf(out_pdf_all, width = 14, height = 5)
par(mfrow = c(1, 3))
plotAnnoPie(peakAnno_all, main = "All peaks: Annotation Pie")
plotAnnoBar(peakAnno_all, main = "All peaks: Annotation Bar")
plotDistToTSS(peakAnno_all, title = "All peaks: Distance to TSS")
dev.off()
message("🖼️ Saved ALL-peaks 3-panel plot to: ", out_pdf_all)

## -----------------------
## Group-specific plots:
##   KO_enriched and WT_enriched
##   (significant peaks only)
## -----------------------
groups_to_plot <- c("KO_enriched", "WT_enriched")

for (grp in groups_to_plot) {

  message("🎯 Preparing group: ", grp)

  # Subset GRanges by direction + significance
  dir_vec <- mcols(peak_gr)$direction
  sig_vec <- mcols(peak_gr)$is_sig

  # Coerce sig_vec robustly (in case it's not logical)
  if (!is.logical(sig_vec)) {
    sig_vec <- tolower(as.character(sig_vec)) %in% c("true", "1", "yes")
  }

  idx_grp <- (dir_vec == grp) & sig_vec

  if (!any(idx_grp)) {
    message("⚠️ No significant peaks found for group: ", grp, ". Skipping plots.")
    next
  }

  peak_gr_grp <- peak_gr[idx_grp]

  message("🧬 Annotating ", sum(idx_grp), " ", grp, " peaks...")

  peakAnno_grp <- annotatePeak(
    peak_gr_grp,
    TxDb      = txdb,
    tssRegion = c(-10000, 10000),
    annoDb    = "org.Mm.eg.db"
  )

  # Custom title and filename
  nice_name <- switch(
    grp,
    "KO_enriched" = "KO-enriched peaks",
    "WT_enriched" = "WT-enriched peaks",
    grp
  )

  out_pdf_grp <- paste0(out_prefix, "_", grp, "_annotation_plots_3panel.pdf")

  pdf(out_pdf_grp, width = 14, height = 5)
  par(mfrow = c(1, 3))
  plotAnnoPie(peakAnno_grp, main = paste0(nice_name, ": Annotation Pie"))
  plotAnnoBar(peakAnno_grp, main = paste0(nice_name, ": Annotation Bar"))
  plotDistToTSS(peakAnno_grp, title = paste0(nice_name, ": Distance to TSS"))
  dev.off()

  message("🖼️ Saved 3-panel plot for ", grp, " to: ", out_pdf_grp)
}

cat("✅ Annotation complete for ", peak_file, "\n", sep = "")
