suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

# --- ARGUMENT ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript peak_log2fold.R <out_dir>")
out_dir  <- args[1]

setwd(out_dir)
out_dir = getwd()
# --- FILE PATHS ---
shared_file <- file.path(out_dir , "shared_peaks_log2FC_full_summary.bed")
wt_file     <- file.path(out_dir , "wt_only_peaks_log2FC_full_summary.bed")
ko_file     <- file.path(out_dir , "ko_only_peaks_log2FC_full_summary.bed")

# --- Safe reader ---
safe_read <- function(path, colnames, type) {
  if (!file.exists(path)) {
    warning(sprintf("⚠️  File not found: %s — skipping %s peaks", path, type))
    return(NULL)
  }
  if (file.info(path)$size == 0) {
    warning(sprintf("⚠️  File is empty: %s — skipping %s peaks", path, type))
    return(NULL)
  }
  df <- fread(path, header = FALSE)
  setnames(df, colnames)
  df$Peak <- type
  return(df)
}

# --- Define column names ---
shared_cols <- c("Chr", "Start", "End", 
  "WT_Total_Signal", "WT_Max_Signal", "WT_Max_Pos",
  "KO_Total_Signal", "KO_Max_Signal", "KO_Max_Pos",
  "WT_Mean", "WT_Median" , "WT_Sum", "WT_SD", "WT_Max", "WT_Min",
  "KO_Mean", "KO_Median", "KO_Sum", "KO_SD", "KO_Max", "KO_Min", "log2fc")

exclusive_cols <- c("Chr", "Start", "End",
  "Total_Signal", "Max_Signal", "Max_Pos",
  "WT_Mean", "WT_Median", "WT_Sum", "WT_SD", "WT_Max", "WT_Min",
  "KO_Mean", "KO_Median", "KO_Sum", "KO_SD", "KO_Max", "KO_Min", "log2fc")

# --- LOAD DATASETS SAFELY ---
shared  <- safe_read(shared_file, shared_cols, "shared")
wt_only <- safe_read(wt_file, exclusive_cols, "WT-only")
ko_only <- safe_read(ko_file, exclusive_cols, "KO-only")

# --- Combine if any data exists ---
data_list <- list(shared, wt_only, ko_only)
data_list <- Filter(Negate(is.null), data_list)

if (length(data_list) == 0) stop("❌ No valid peak files found. Exiting.")

combined <- rbindlist(data_list, fill = TRUE)

# --- SAVE TABLE ---
fwrite(combined, file = file.path(out_dir , "combined_peaks_log2FC.csv"), sep = ",")
# --- SAVE BED FILE ---
bed_file <- file.path(out_dir, "combined_peaks_log2FC.bed")

# Ensure required BED fields exist; if not, create dummy values
combined_bed <- combined[, .(Chr, Start, End)]
combined_bed[, name := paste0("peak_", .I)]
combined_bed[, score := 0]
combined_bed[, strand := "."]

# Write BED file
fwrite(combined_bed, bed_file, sep = "\t", col.names = FALSE)



# --- ORDER GROUPS ---
# combined[, group := fifelse(!is.na(WT_Total_Signal) & !is.na(KO_Total_Signal), "Shared",
#                             fifelse(!is.na(WT_Total_Signal), "WT-only", "KO-only"))]
# group_counts <- combined[, .N, by = group][order(N)]
# combined[, group := factor(group, levels = group_counts$group)]

# --- PLOT: DENSITY ---
p1 <- ggplot(combined, aes(x = log2fc)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(title = "log2FC Density (All Peaks)", x = "log2(KO / WT)", y = "Density") +
  theme_minimal()

p2 <- ggplot(combined, aes(x = log2fc, fill = Peak, color = Peak)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(title = "log2FC Density by Peak Group", x = "log2(KO / WT)", y = "Density") +
  theme_minimal()

ggsave(file.path(out_dir , "log2FC_two_panel_density.pdf"), plot = p1 / p2, width = 6, height = 6)

# --- PLOT: HISTOGRAM ---
h1 <- ggplot(combined, aes(x = log2fc)) +
  geom_histogram(binwidth = 0.25, fill = "steelblue", alpha = 0.6, color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(title = "log2FC Histogram (All Peaks)", x = "log2(KO / WT)", y = "Count") +
  theme_minimal()

h2 <- ggplot(combined, aes(x = log2fc, fill = Peak, color = Peak)) +
  geom_histogram(binwidth = 0.25, alpha = 0.4, position = "stack") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(title = "log2FC Histogram by Group", x = "log2(KO / WT)", y = "Count") +
  theme_minimal()

ggsave(file.path(out_dir , "log2FC_two_panel_histogram.pdf"), plot = h1 / h2, width = 6, height = 6)

# --- DUAL-PANEL: Density & Histogram ---
final <- p2 / h2
ggsave(file.path(out_dir , "log2FC_distribution.pdf"), plot = final, width = 6, height = 6)

cat("✅ Plots and table saved to:", out_dir , "\n")
