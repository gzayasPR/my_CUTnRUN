#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
aligment_results <- args[1]
out_file <- args[2]

message("Alignment results: ", aligment_results)
message("Output file: ", out_file)

x <- read.table(aligment_results, header = TRUE, sep = "\t", check.names = FALSE)
names(x) <- c("Project","Sample","Raw_Reads","TrimmedReads","MouseAligned","EcoliAligned",
              "TotalAligned","MouseClean","EcoliClean","TotalClean","MouseAlignRate","PortionOfEcoli")

# keep relevant
d <- x[, c("Sample","TrimmedReads","MouseClean","EcoliClean","TotalClean")]
# flag IgG / controls
d$is_IgG <- grepl("IgG", d$Sample, ignore.case = TRUE)

# sanity
d$MouseClean <- as.numeric(d$MouseClean)
d$EcoliClean <- as.numeric(d$EcoliClean)
if (anyNA(d$MouseClean) || anyNA(d$EcoliClean)) stop("MouseClean/EcoliClean has NA.")
if (any(d$MouseClean <= 0)) stop("MouseClean has zero or negative values; cannot compute sRPM.")

# spike-in per million mouse reads (stable)
d$spike_sRPM <- (d$EcoliClean / d$MouseClean) * 1e6

# choose robust reference among non-IgG
ref_pool <- d$spike_sRPM[!d$is_IgG & is.finite(d$spike_sRPM)]
if (length(ref_pool) < 2) stop("Not enough non-IgG samples to define reference.")
ref_value <- median(ref_pool, na.rm = TRUE)

# raw scale factor: bring all samples to the median spike-in level
d$scale_raw <- ref_value / d$spike_sRPM

# clamp extremes to avoid over-scaling from outliers (tweak if desired)
sf_min <- 1; sf_max <- 150
d$scale_clamped <- pmin(sf_max, pmax(sf_min, d$scale_raw))
d$clamped <- d$scale_clamped != d$scale_raw

# helpful extras
d$percent_spikein <- d$EcoliClean / (d$EcoliClean + d$MouseClean)
d$reference_spike_sRPM <- ref_value
d$notes <- ifelse(d$clamped, "clamped", "ok")

# arrange columns for your downstream script
# out <- d[, c("Sample","MouseClean","EcoliClean","TotalClean",
#              "spike_sRPM","percent_spikein","scale_raw","scale_clamped","notes")]
out <- d[, c("Sample","MouseClean","EcoliClean","TotalClean",
             "spike_sRPM","percent_spikein", "scale_clamped")]

write.table(out, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote: ", out_file)
