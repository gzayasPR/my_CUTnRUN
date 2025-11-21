# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

aligment_results  <- args[1]
out_file    <- args[2]

# aligment_results <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project5_CREB1_V2/2.Alignment/alignment_summary.tsv"
# out_file    <- "/blue/mateescu/gzayas97/CUTRUN_Feng_Project/results/Project5_CREB1_V2/3.Normalization/scale_factors.tsv"

print(paste("Alignment results:",aligment_results))
print(paste("Output file:", out_file))
# --------------------------
align_results <- read.table(aligment_results, header = TRUE, sep = "\t")
names(align_results) <- c("Prject","Sample","Raw_Reads","Total","Mouse","Mouse_Alignment_Rate","Clean_Reads","Ecoli","Prop_Ecoli")
col_keep <- c("Sample","Total","Mouse","Ecoli","Clean_Reads")
data <- align_results[col_keep]
# Add percent Ecoli (optional, for reference)
# data$percent <- data$Ecoli / (data$Ecoli + data$Mouse)
data$Total <- data$Total/2
# data$percent <- data$Ecoli / (data$Total)
data$percent <- data$Ecoli / (data$Clean_Reads + data$Ecoli)
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
# max_percent <- max(data$percent, na.rm = TRUE)
# Scale factor = median(Ecoli) / Ecoli_i
# data$scale_factor <- (max_percent * 100)  / (data$percent * 100) 
# data$scale_factor <- 1 / (data$percent * 100) 

max_ecoli_reads <- max(non_IgG_data$Ecoli, na.rm = TRUE)
data$scale_factor <- max_ecoli_reads / (data$Ecoli) 
# Write output
write.table(data, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)