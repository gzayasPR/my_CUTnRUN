# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

file_path      <- args[1]  # genome-wide log2FC bedgraph
out_dir        <- args[2]
shared_file    <- args[3]  # shared peaks + log2FC
wt_only_file   <- args[4]  # WT-only peaks + log2FC
ko_only_file   <- args[5]  # KO-only peaks + log2FC

library(data.table)
library(ggplot2)

setwd(out_dir)

# ------------------------
# 📊 Genome-wide Density Plot
# ------------------------

data <- fread(file_path, header = FALSE)
colnames(data) <- c("chr", "start", "end", "log2FC")
print("Loaded genome-wide log2FC")
print(head(data))

png("log2FC_density_plot.png")
plot(density(data$log2FC, na.rm = TRUE),
     main = "Genome-wide log2 Fold Change",
     xlab = "log2(KO / WT)",
     ylab = "Density",
     col = "blue",
     lwd = 2)
abline(v = c(-1,1), col = "red", lty = 2)
dev.off()

# ------------------------
# 📊 SEACR Peak-Type Density Plot (shared, WT-only, KO-only)
# ------------------------
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)  # for multi‐panel layouts
})

# -- Helper: attempt to fread; if file is missing or empty, return an empty data.table 
safe_fread <- function(path) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path, "→ returning empty data.table()"))
    # Return an empty table with at least 4 columns: V1,V2,V3, log2FC
    return( data.table(V1 = character(0),
                       V2 = integer(0),
                       V3 = integer(0),
                       log2FC = numeric(0)) )
  }
  dt <- tryCatch(fread(path, header = FALSE), 
                 error = function(e) {
                   warning(paste("Error reading", path, "→ returning empty data.table()"))
                   data.table()
                 })
  if (nrow(dt) == 0L) {
    # No rows at all
    if (ncol(dt) >= 4) {
      # If it still has columns, coerce to at least V1,V2,V3,log2FC
      return( data.table(V1 = character(0),
                         V2 = integer(0),
                         V3 = integer(0),
                         log2FC = numeric(0)) )
    } else {
      # If zero columns (even weird), force these four
      return( data.table(V1 = character(0),
                         V2 = integer(0),
                         V3 = integer(0),
                         log2FC = numeric(0)) )
    }
  }
  # If we got some data, ensure there's a log2FC column at the end:
  if (ncol(dt) < 4) {
    # This shouldn’t normally happen if you expect at least 4 cols, but guard just in case:
    dt[, log2FC := NA_real_]
  } else {
    # Rename the last column to “log2FC”
    setnames(dt, old = names(dt)[ncol(dt)], new = "log2FC")
  }
  # If columns 1:3 exist but not named V1,V2,V3, rename them so we can refer uniformly below
  if (!all(c("V1", "V2", "V3") %in% names(dt))) {
    old_names <- names(dt)
    # assume first three columns are chr,start,end if they aren’t already V1,V2,V3
    setnames(dt, old = old_names[1:3], new = c("V1", "V2", "V3"))
    # everything else (except the final log2FC) can be left alone
  }
  return(dt)
}

# — Replace these with your actual paths —
shared_file  <- "path/to/shared_file.bed"
wt_only_file <- "path/to/wt_only_file.bed"
ko_only_file <- "path/to/ko_only_file.bed"

# → Read (or get empty dt) for each group:
shared_dt <- safe_fread(shared_file)
wt_only_dt <- safe_fread(wt_only_file)
ko_only_dt <- safe_fread(ko_only_file)

# If any of those dt’s has no ‘log2FC’ column, add it explicitly:
for (dt in list(shared_dt, wt_only_dt, ko_only_dt)) {
  if (!"log2FC" %in% names(dt)) {
    dt[, log2FC := numeric( .N )] 
  }
}

# Build combined table via rbindlist, filling in missing columns automatically:
combined <- rbindlist(list(
  shared   = shared_dt[, .(chr = V1, start = as.integer(V2), end = as.integer(V3), log2FC)],
  `WT-only`= wt_only_dt[, .(chr = V1, start = as.integer(V2), end = as.integer(V3), log2FC)],
  `KO-only`= ko_only_dt[, .(chr = V1, start = as.integer(V2), end = as.integer(V3), log2FC)]
), idcol = "group", fill = TRUE)

# If any of the three original tables was empty, that group will simply contribute 0 rows.
# Set 'group' as a factor with a consistent ordering:
combined[, group := factor(group, levels = c("shared", "WT-only", "KO-only"))]
combined[, group := fcase(
  group == "shared",   "Shared",
  group == "WT-only",  "WT-only",
  group == "KO-only",  "KO-only"
)]

# Save the full annotated table:
fwrite(
  combined,
  file = "combined_peaks_log2FC.bed",
  quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = TRUE
)
message("Wrote combined_peaks_log2FC.bed")

# Count total entries per group (this will skip empty groups transparently):
group_counts <- combined[, .N, by = group][order(-N)]
print(group_counts)

# If you prefer to reorder by ascending N, simply use `order(N)`.
# Re‐factor 'group' so that ggplot layers them in ascending or descending order:
combined[, group := factor(group, levels = group_counts$group)]

# ------------------------
# 🟦 Top panel: Density without grouping
p1 <- ggplot(combined, aes(x = log2FC)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(
    title = "log2FC Distribution (All SEACR Peaks)",
    x = "log2(KO / WT)", y = "Density"
  ) +
  theme_minimal()

# 🟨 Bottom panel: Grouped by SEACR peak type
p2 <- ggplot(combined, aes(x = log2FC, fill = group, color = group)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(
    title = "log2FC Distribution by SEACR Peak Type",
    x = "log2(KO / WT)", y = "Density"
  ) +
  theme_minimal()

# 🖼 Combine the two density plots vertically:
final_plot <- p1 / p2
ggsave("log2FC_two_panel_density.pdf", plot = final_plot, width = 6, height = 6)

# ------------------------
# 🧱 Top panel: Histogram without grouping
h1 <- ggplot(combined, aes(x = log2FC)) +
  geom_histogram(binwidth = 0.25, fill = "steelblue", alpha = 0.6, color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(
    title = "log2FC Histogram (All SEACR Peaks)",
    x = "log2(KO / WT)", y = "Count"
  ) +
  theme_minimal()

# 🧱 Second panel: Stacked histogram by group
h2 <- ggplot(combined, aes(x = log2FC, fill = group)) +
  geom_histogram(binwidth = 0.25, alpha = 0.4, position = "stack") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  labs(
    title = "log2FC Histogram by SEACR Peak Type (Stacked, Ordered)",
    x = "log2(KO / WT)", y = "Count"
  ) +
  theme_minimal()

# 🖼 Combine histogram panels:
hist_plot <- h1 / h2
ggsave("log2FC_two_panel_histogram.pdf", plot = hist_plot, width = 6, height = 6)

# 🖼 Optionally combine p2 and h2 for a single “distribution” output:
combo <- p2 / h2
ggsave("log2FC_distribution.pdf", plot = combo, width = 6, height = 6)
