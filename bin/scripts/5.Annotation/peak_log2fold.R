#!/usr/bin/env Rscript

# --- DESCRIPTION ---
# Combines KO-only, WT-only, and Shared peaks, and calculates log2FC from `collapse` values
# Generates multiple log2FC visualizations

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
    library(stringr)
})

# --- ARGUMENT ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript peak_log2fold.R <out_dir >")
out_dir  <- args[1]


# out_dir  = "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/5.PeakComparison/KO_H3_vs_WT_H3/Enrichment"
# --- FILE PATHS ---

shared_file = file.path(out_dir , "shared_peaks_log2FC_full_summary.bed")
wt_file  = file.path(out_dir , "wt_only_peaks_log2FC_full_summary.bed")
ko_file  = file.path(out_dir , "ko_only_peaks_log2FC_full_summary.bed")

# --- LOAD ---
shared <- fread(shared_file, header = FALSE)
wt_only     <- fread(wt_file, header = FALSE)
ko_only     <- fread(ko_file, header = FALSE)


# --- Define column names ---
colnames(shared) <- c(
  "Chr", "Start", "End" , 
  "WT_Total_Signal", "WT_Max_Signal", "WT_Max_Pos",
  "KO_Total_Signal", "KO_Max_Signal", "KO_Max_Pos",
  "Collapsed_WT","Collapsed_KO")

colnames(ko_only  ) <- c(
  "Chr", "Start", "End" ,
  "KO_Total_Signal", "KO_Max_Signal", "KO_Max_Pos",
  "Collapsed_WT","Collapsed_KO")

colnames(wt_only  ) <- c(
  "Chr", "Start", "End" ,
  "WT_Total_Signal", "WT_Max_Signal", "WT_Max_Pos",
  "Collapsed_WT","Collapsed_KO")

shared$Peak <- "shared"
ko_only$Peak <- "KO-only"
wt_only$Peak <- "WT-only"
# --- COMBINE ---
combined <- rbind(shared, wt_only, ko_only, fill = TRUE)

# --- Helper: split collapsed vector and convert to numeric ---
parse_vector <- function(x) as.numeric(str_split(x, ",", simplify = TRUE))

# --- Custom stats ---
safe_stat <- function(vec, stat = "mean") {
  vec <- as.numeric(vec)
#   if (stat %in% c("mean", "median")) {
#     vec <- vec[vec != 0]  # Remove zeros
#   }
  if (length(vec) == 0) return(NA_real_)
  switch(stat,
         mean   = mean(vec),
         median = median(vec),
         sum    = sum(vec),
         sd     = sd(vec),
         max    = max(vec),
         min    = min(vec),
         NA_real_)
}

# --- Apply to both WT and KO ---
metrics <- c("mean", "median", "sum", "sd", "max", "min")

for (m in metrics) {
  combined[[paste0("KO_", m)]] <- sapply(combined$Collapsed_KO, function(x) {
    safe_stat(str_split(x, ",", simplify = TRUE), m)
  })
  
  combined[[paste0("WT_", m)]] <- sapply(combined$Collapsed_WT, function(x) {
    safe_stat(str_split(x, ",", simplify = TRUE), m)
  })
}

# --- Calculate log2 fold change using mean ---
combined[, log2fc := log2((KO_sum) / (WT_sum))]

# --- DROP collapsed columns ---
combined[, c("Collapsed_WT", "Collapsed_KO") := NULL]

# --- SAVE TABLE ---
fwrite(combined, file = file.path(out_dir , "combined_peaks_log2FC.csv"), sep = ",")
cat("✅ Saved combined log2FC table to combined_peaks_log2FC.csv\n")

# --- ORDER GROUPS ---
combined[, group := fifelse(!is.na(WT_Total_Signal) & !is.na(KO_Total_Signal), "Shared",
                            fifelse(!is.na(WT_Total_Signal), "WT-only", "KO-only"))]
group_counts <- combined[, .N, by = group][order(N)]
combined[, group := factor(group, levels = group_counts$group)]

# --- PLOT: DENSITY ---
p1 <- ggplot(combined, aes(x = log2fc)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(title = "log2FC Density (All Peaks)", x = "log2(KO / WT)", y = "Density") +
  theme_minimal()

p2 <- ggplot(combined, aes(x = log2fc, fill = group, color = group)) +
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

h2 <- ggplot(combined, aes(x = log2fc, fill = group, color = group)) +
  geom_histogram(binwidth = 0.25, alpha = 0.4, position = "stack") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(title = "log2FC Histogram by Group", x = "log2(KO / WT)", y = "Count") +
  theme_minimal()

ggsave(file.path(out_dir , "log2FC_two_panel_histogram.pdf"), plot = h1 / h2, width = 6, height = 6)

# --- DUAL-PANEL: Density & Histogram ---
final <- p2 / h2
ggsave(file.path(out_dir , "log2FC_distribution.pdf"), plot = final, width = 6, height = 6)

cat("✅ Plots saved to:", out_dir , "\n")
