#!/usr/bin/env Rscript

# --- Load libraries ---
suppressMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(tools)
  library(data.table)
})

# # --- Input arguments ---
args <- commandArgs(trailingOnly = TRUE)

peak_file <- args[1]
gtf_file <- args[2]
fold_change_file <- args[3]
out_dir <- args[4]

# peak_file <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/Enrichment/combined_peaks_log2FC.bed"
# gtf_file <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/data/references/Mus_musculus_GRCm39/Ensembl/111/Annotation/Genes/Mus_musculus.GRCm39.111.gtf"
# out_dir <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/Enrichment"
# fold_change_file <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project2_Dpy30/5.PeakComparison/WT_vs_KO/Enrichment/combined_peaks_log2FC.csv"
# --- Set TxDb path ---
txdb_file <- file_path_sans_ext(basename(gtf_file))
txdb_sqlite <- paste0(txdb_file, "_ensembl.sqlite")

# --- Build or load TxDb ---
if (!file.exists(txdb_sqlite)) {
  message("🔧 Building TxDb from GTF...")
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  saveDb(txdb, file = txdb_sqlite)
} else {
  txdb <- loadDb(txdb_sqlite)
}

# --- Load and annotate peaks ---
message("📂 Annotating peaks from: ", peak_file)
peak <- readPeakFile(peak_file)
suppressWarnings({
  peakAnno <- annotatePeak(peak, TxDb = txdb,
                           tssRegion = c(-10000, 10000),
                           annoDb = "org.Mm.eg.db")
})

# --- Convert annotated peaks to data.frame ---
# Convert annotation to data.table
# --- Load and prepare fc_df ---
fc_df <- fread(fold_change_file)
fc_df[, Chr := as.character(Chr)]
fc_df[, Start := as.integer(Start)]
fc_df[, End := as.integer(End)]

# --- Prepare anno_df ---
anno_df <- as.data.table(as.data.frame(peakAnno))
anno_df[, Chr := as.character(seqnames)]
anno_df[, Start := as.integer(start) - 1L]
anno_df[, End := as.integer(end)]

# --- Merge on coordinates ---
merged_df <- merge(fc_df, anno_df, by = c("Chr", "Start", "End"), all= TRUE)
keep_cols <- c(names(fc_df), "geneStart",
                "geneEnd", "geneLength",
                 "geneStrand", "geneId",
                 "transcriptId", "distanceToTSS",
                  "ENTREZID", "SYMBOL", "GENENAME"  )
clean_df <- merged_df[, intersect(keep_cols, names(merged_df)), with = FALSE]

# --- Save merged file ---
output_file <- file.path(out_dir, "combined_peaks_log2FC_annotated_merged.csv")
fwrite(clean_df, output_file)

cat("✅ Merged annotation saved to:", output_file, "\n")
