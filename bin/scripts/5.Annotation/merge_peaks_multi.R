#!/usr/bin/env Rscript

## ------------------------------------------------------------
## merge_multiinter_peaks.R
##
## Merge overlapping SEACR peaks across multiple samples,
## similar to bedtools multiinter, but:
##  - use sample names like WT1, WT2, KO1, KO2
##  - output per-sample presence
##  - add WT_any, KO_any, Group (Shared / WT-only / KO-only)
##
## Usage:
##   Rscript merge_multiinter_peaks.R samples.csv multiinter_raw.tsv
##
## samples.csv must have at least:
##   sample_name,bed_path
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(methods)
})

## ------------------------------------------------------------
## Helper: read sample table (CSV)
##   - col1: sample_name (WT1, WT2, KO1, KO2, etc.)
##   - col2: bed_path
## ------------------------------------------------------------
read_sample_table <- function(csv_path) {
  tab <- read.csv(csv_path, stringsAsFactors = FALSE, header = TRUE)
  if (ncol(tab) < 2) {
    stop("Sample CSV must have at least two columns: sample_name, bed_path")
  }
  colnames(tab)[1:2] <- c("sample_name", "bed_path")
  if (any(is.na(tab$sample_name) | tab$sample_name == "")) {
    stop("Found empty sample_name in first column of CSV.")
  }
  if (any(is.na(tab$bed_path) | tab$bed_path == "")) {
    stop("Found empty bed_path in second column of CSV.")
  }
  tab
}

## ------------------------------------------------------------
## Helper: read SEACR BED-like files into GRanges
##
## Expected format:
##   V1: chr
##   V2: start (0-based)
##   V3: end   (1-based, half-open)
##   V4+: extra cols (signal, max, summit string, etc.)
##
## We convert to 1-based GRanges internally:
##   GRanges start = BED_start + 1
##   GRanges end   = BED_end
## ------------------------------------------------------------
read_bed_as_granges_from_table <- function(sample_tab) {
  sample_names <- sample_tab$sample_name
  bed_files    <- sample_tab$bed_path
  
  gr_list <- vector("list", length(bed_files))
  names(gr_list) <- sample_names
  
  for (i in seq_along(bed_files)) {
    bed_path <- bed_files[i]
    if (!file.exists(bed_path)) {
      stop("BED file does not exist: ", bed_path)
    }
    
    bed_dt <- read.table(
      bed_path,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
    if (ncol(bed_dt) < 3) {
      stop("BED file must have at least 3 columns (chr, start, end): ", bed_path)
    }
    
    chr    <- bed_dt[[1]]
    start0 <- as.integer(bed_dt[[2]])  # 0-based start
    end    <- as.integer(bed_dt[[3]])  # BED end
    
    gr <- GRanges(
      seqnames = chr,
      ranges   = IRanges(start = start0 + 1L, end = end)
    )
    
    # Extra columns as metadata (signal, summit, etc.)
    if (ncol(bed_dt) > 3) {
      mcols(gr) <- bed_dt[, -(1:3), drop = FALSE]
    }
    mcols(gr)$sample_id <- sample_names[i]
    
    gr_list[[i]] <- gr
  }
  
  gr_list
}

## ------------------------------------------------------------
## Helper: compute presence/absence matrix
##
## presence_mat[union_region, sample] = 1 if sample has a peak
## overlapping that union region (with optional overlap thresholds).
## ------------------------------------------------------------
compute_presence_matrix <- function(union_gr,
                                    gr_list,
                                    min_overlap_frac = 0,
                                    reciprocal = FALSE,
                                    ignore_strand = TRUE) {
  n_union   <- length(union_gr)
  samples   <- names(gr_list)
  n_samples <- length(samples)
  
  presence_mat <- matrix(0L, nrow = n_union, ncol = n_samples)
  colnames(presence_mat) <- samples
  
  for (i in seq_along(gr_list)) {
    gr <- gr_list[[i]]
    if (length(gr) == 0L) next
    
    hits <- findOverlaps(union_gr, gr, ignore.strand = ignore_strand)
    if (length(hits) == 0L) next
    
    q_idx <- queryHits(hits)
    s_idx <- subjectHits(hits)
    
    if (min_overlap_frac > 0) {
      ov_gr    <- pintersect(union_gr[q_idx], gr[s_idx], drop.nohit.ranges = TRUE)
      ov_width <- width(ov_gr)
      frac_union  <- ov_width / width(union_gr[q_idx])
      frac_sample <- ov_width / width(gr[s_idx])
      
      if (reciprocal) {
        keep <- frac_union >= min_overlap_frac & frac_sample >= min_overlap_frac
      } else {
        keep <- frac_union >= min_overlap_frac
      }
      
      if (!any(keep)) next
      q_idx <- q_idx[keep]
    }
    
    presence_mat[unique(q_idx), i] <- 1L
  }
  
  presence_mat
}

## ------------------------------------------------------------
## Core: merge peaks across samples and build multiinter-like table
## ------------------------------------------------------------
merge_peaks_multi_from_table <- function(sample_tab,
                                         maxgap = 0L,
                                         ignore_strand = TRUE,
                                         min_overlap_frac = 0,
                                         reciprocal = FALSE) {
  ## 1. BEDs -> GRanges
  gr_list <- read_bed_as_granges_from_table(sample_tab)
  
  # sanity check
  is_gr <- vapply(gr_list, function(x) is(x, "GenomicRanges"), logical(1))
  if (!all(is_gr)) {
    stop("Internal error: some entries in gr_list are not GRanges: ",
         paste(names(gr_list)[!is_gr], collapse = ", "))
  }
  
  ## 2. Union of all peaks and reduction
  all_gr <- do.call(c, unname(gr_list))
  if (!is(all_gr, "GenomicRanges")) {
    all_gr <- unlist(as(gr_list, "GRangesList"))
  }
  
  union_gr <- GenomicRanges::reduce(
    all_gr,
    min.gapwidth  = as.integer(maxgap),
    ignore.strand = ignore_strand
  )
  
  ## 3. Presence/absence matrix
  presence_mat <- compute_presence_matrix(
    union_gr        = union_gr,
    gr_list         = gr_list,
    min_overlap_frac = min_overlap_frac,
    reciprocal      = reciprocal,
    ignore_strand   = ignore_strand
  )
  
  ## 4. n_files + which_files (using sample names, e.g., WT1, WT2, KO1)
  n_files <- rowSums(presence_mat)
  
  which_files <- apply(presence_mat, 1, function(x) {
    if (!any(x)) return(NA_character_)
    paste(colnames(presence_mat)[x == 1L], collapse = ",")
  })
  
  ## 5. WT_any, KO_any, Group (Shared / WT-only / KO-only)
  col_names <- colnames(presence_mat)
  wt_idx <- grep("^WT", col_names)
  ko_idx <- grep("^KO", col_names)
  
  WT_any <- if (length(wt_idx) > 0) {
    as.integer(rowSums(presence_mat[, wt_idx, drop = FALSE]) > 0)
  } else {
    integer(nrow(presence_mat))
  }
  
  KO_any <- if (length(ko_idx) > 0) {
    as.integer(rowSums(presence_mat[, ko_idx, drop = FALSE]) > 0)
  } else {
    integer(nrow(presence_mat))
  }
  
  Group <- rep("NA", nrow(presence_mat))
  Group[WT_any == 1 & KO_any == 1] <- "Shared"
  Group[WT_any == 1 & KO_any == 0] <- "WT-only"
  Group[WT_any == 0 & KO_any == 1] <- "KO-only"
  
  ## 6. Build output data.frame (back to 0-based BED start)
  out_df <- data.frame(
    chr         = as.character(GenomicRanges::seqnames(union_gr)),
    start       = GenomicRanges::start(union_gr) - 1L,  # 0-based
    end         = GenomicRanges::end(union_gr),
    n_files     = n_files,
    which_files = which_files,
    check.names = FALSE
  )
  
  presence_df <- as.data.frame(presence_mat, check.names = FALSE)
  
  out_df <- cbind(out_df,
                  presence_df,
                  WT_any = WT_any,
                  KO_any = KO_any,
                  Group  = Group)
  
  out_df
}

## ------------------------------------------------------------
## Main CLI
## ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage:\n",
      "  Rscript merge_multiinter_peaks.R samples.csv multiinter_raw.tsv\n\n",
      "samples.csv must have at least:\n",
      "  sample_name,bed_path\n",
      sep = "")
  quit(status = 1)
}

sample_csv <- args[1]
out_file   <- args[2]

## ---- Tunable defaults (you can edit these) ----
maxgap           <- 0L      # e.g. 50 to merge peaks within 50 bp
min_overlap_frac <- 0       # e.g. 0.2 => >=20% of union region
reciprocal       <- FALSE   # only matters if min_overlap_frac > 0
ignore_strand    <- TRUE    # no strand info in SEACR peaks

sample_tab <- read_sample_table(sample_csv)

merged_peaks <- merge_peaks_multi_from_table(
  sample_tab       = sample_tab,
  maxgap           = maxgap,
  ignore_strand    = ignore_strand,
  min_overlap_frac = min_overlap_frac,
  reciprocal       = reciprocal
)

write.table(
  merged_peaks,
  file      = out_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
