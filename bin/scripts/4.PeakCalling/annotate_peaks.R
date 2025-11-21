#!/usr/bin/env Rscript

# --- Load libraries ---
suppressMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(tools)
})

# --- Input arguments ---
args <- commandArgs(trailingOnly = TRUE)
peak_file <- args[1]
out_prefix <- args[2]
gtf_file <- args[3]

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
peakAnno <- annotatePeak(peak, TxDb = txdb, tssRegion = c(-10000,10000), annoDb = "org.Mm.eg.db")

# --- Save annotation table ---
out_txt <- paste0(out_prefix, "_annotated.txt")
write.table(as.data.frame(peakAnno), out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
message("📄 Saved annotation table to: ", out_txt)

# --- 3-panel base R plot ---
out_pdf <- paste0(out_prefix, "_annotation_plots_3panel.pdf")
pdf(out_pdf, width = 14, height = 5)
par(mfrow = c(1, 3))  # 3 columns, adjust margins

plotAnnoPie(peakAnno, main = "Annotation Pie")
plotAnnoBar(peakAnno, main = "Annotation Bar")
plotDistToTSS(peakAnno, title = "Distance to TSS")
dev.off()
message("🖼️ Saved 3-panel annotation plot to: ", out_pdf)

cat("✅ Annotation complete for", peak_file, "\n")
