# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

ecoli.count  <- args[1]
mouse.count <- args[2]
total.count <- args[3]
out_file    <- args[4]

# ecoli.count  <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/3.Normalization/ecoli_counts.tsv"
# mouse.count <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/3.Normalization/mouse_counts.tsv"
# out_file    <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Run2/3.Normalization/scale_factors.tsv"

print(paste("Ecoli count file:", ecoli.count))
print(paste("Mouse count file:", mouse.count))
print(paste("Output file:", out_file))
# --------------------------
total_n <- read.table(total.count, header = TRUE, sep = "\t")
names(total_n) [2] <- "Total"
total_n$Total <- total_n$Total * 2
ecoli_n <- read.table(ecoli.count, header = TRUE, sep = "\t")
names(ecoli_n)[2] <- "Ecoli"
mouse_n <- read.table(mouse.count, header = TRUE, sep = "\t")
names(mouse_n)[2] <- "Mouse"

data1 <- merge(ecoli_n, mouse_n, by = "Sample")
data <- merge(data1, total_n, by = "Sample")
# Add percent Ecoli (optional, for reference)
# data$percent <- data$Ecoli / (data$Ecoli + data$Mouse)
data$percent <- data$Ecoli / (data$Total)
# Calculate scale factor: scale_factor = max(Ecoli) / Ecoli
# max_ecoli <- max(data$Ecoli)
# max_percent <- max(data$percent)
# data$scaled <- max_percent / data$percent 
# Compute median spike-in count
# Filter out samples that contain "IgG" in their names
non_IgG_data <- data[!grepl("IgG", data$Sample, ignore.case = TRUE), ]
median_percent <- median(data$percent, na.rm = TRUE)
# Compute max percent only among non-IgG samples
max_percent <- max(non_IgG_data$percent, na.rm = TRUE)
# Scale factor = median(Ecoli) / Ecoli_i
data$scale_factor <- (max_percent * 100)  / (data$percent * 100) 
# data$scale_factor <- 1 / (data$percent * 100) 
# data$scale_factor <- 1 / (data$percent * 100) 
# Write output
write.table(data, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)