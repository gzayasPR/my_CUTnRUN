#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# CUT&RUN spike-in scaling (E. coli normalized by Mouse reads)
# - Computes spike_sRPM = (EcoliClean / MouseClean) * 1e6
# - Detects outliers among non-IgG using MAD or IQR (configurable)
# - Chooses a reference target from the non-IgG pool (minus outliers) using a
#   user-specified statistic: max/min/median/mean or a percentile (e.g., p90)
# - Scales all samples to that target; clamps extreme scale factors
# - Writes metadata columns for full provenance
#
# USAGE:
#   Rscript scale_spikein.R <alignment_results.tsv> <out.tsv> [REF_SPEC]
#
#   REF_SPEC (optional; also via env REF_SPEC) can be:
#     - "max" (default), "min", "median", "mean"
#     - Percentiles: "p90", "90%", "q0.9", or "0.9" (interpreted as 90th pct)
#
# ENV OVERRIDES (optional):
#   OUTLIER_METHOD = "MAD" | "IQR"        (default "MAD")
#   MAD_K          = 3.5                  (distance multiplier)
#   IQR_K          = 1.5
#   SF_MIN         = 1                    (min clamp for scale factor)
#   SF_MAX         = 150                  (max clamp)
#   MIN_REF_N      = 2                    (min non-IgG for reference)
#   REF_SPEC       = same as 3rd arg
#
# NOTES:
# - IgG samples are *never* used for the reference pool.
# - Outliers are detected only among non-IgG samples.
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript scale_spikein.R <alignment_results.tsv> <out.tsv> [REF_SPEC]")
}
aligment_results <- args[1]
out_file <- args[2]
# NEW: single spec for reference statistic (string)
# REF_SPEC <- if (length(args) >= 3) args[3] else "max"
REF_SPEC <- "max"
# -------------------- options (env overrides) ---------------------------------
OUTLIER_METHOD <- toupper(Sys.getenv("OUTLIER_METHOD", unset = "MAD")) # MAD or IQR
MAD_K          <- as.numeric(Sys.getenv("MAD_K", unset = "3.5"))
IQR_K          <- as.numeric(Sys.getenv("IQR_K", unset = "1.5"))
SF_MIN         <- as.numeric(Sys.getenv("SF_MIN", unset = "1"))
SF_MAX         <- as.numeric(Sys.getenv("SF_MAX", unset = "150"))
MIN_REF_N      <- as.integer(Sys.getenv("MIN_REF_N", unset = "2"))


message("Alignment results: ", aligment_results)
message("Output file: ", out_file)
message("Outlier method: ", OUTLIER_METHOD, "  |  REF_SPEC: ", REF_SPEC)

# -------------------- helpers -------------------------------------------------

# Detect outliers in a numeric vector using MAD or IQR (robust; NA/Inf safe)
detect_outliers <- function(v, method = "MAD", mad_k = 3.5, iqr_k = 1.5) {
  v <- as.numeric(v)
  v_ok <- is.finite(v)
  out <- rep(FALSE, length(v))
  if (sum(v_ok) < 2) return(out)
  vv <- v[v_ok]

  if (toupper(method) == "MAD") {
    med <- stats::median(vv, na.rm = TRUE)
    m   <- stats::mad(vv, center = med, na.rm = TRUE)  # 1.4826-consistent
    thresh <- mad_k * ifelse(m > 0, m, 0)
    out[v_ok] <- abs(vv - med) > thresh
  } else if (toupper(method) == "IQR") {
    qs <- stats::quantile(vv, probs = c(0.25, 0.75), na.rm = TRUE, type = 7)
    i  <- qs[2] - qs[1]
    if (i > 0) {
      lo <- qs[1] - iqr_k * i
      hi <- qs[2] + iqr_k * i
      out[v_ok] <- vv < lo | vv > hi
    }
  } else {
    stop("Unknown OUTLIER_METHOD: ", method)
  }
  out
}

# NEW: parse a REF_SPEC string into a method and (optional) quantile in [0,1]
# Accepted (case-insensitive):
#   "max","min","median","mean"
#   "p90","90%","q0.9","0.9"  -> 90th percentile
parse_ref_spec <- function(spec) {                                        # NEW
  s <- tolower(trimws(spec))
  # direct names
  if (s %in% c("max","min","median","mean")) {
    return(list(kind = s, q = NA_real_))
  }
  # patterns for percentiles
  s_nopct <- sub("%$", "", s)                # strip trailing %
  if (grepl("^p\\d{1,2}(\\.\\d+)?$", s)) {   # p90, p99.5
    num <- as.numeric(sub("^p", "", s))
    return(list(kind = "percentile", q = num / 100))
  }
  if (grepl("^q\\d*\\.?\\d+$", s)) {         # q0.9
    num <- as.numeric(sub("^q", "", s))
    return(list(kind = "percentile", q = num))
  }
  if (grepl("^\\d*\\.?\\d+$", s_nopct)) {    # "90", "90.0", "0.9"
    num <- as.numeric(s_nopct)
    # Interpret 0–1 as quantile; >1 as 0–100 percentage
    q <- if (num <= 1) num else num / 100
    return(list(kind = "percentile", q = q))
  }
  stop("REF_SPEC '", spec, "' not understood. Use max|min|median|mean or a percentile (e.g., p90, 90%, q0.9, 0.9).")
}

# NEW: compute the reference value given candidates and the parsed spec
compute_ref_value <- function(vals, spec_parsed) {                        # NEW
  vals <- vals[is.finite(vals)]
  if (!length(vals)) stop("No finite values in reference pool.")
  k <- spec_parsed$kind
  if (k == "max")    return(list(value = max(vals),    method = "max",    q = NA_real_))
  if (k == "min")    return(list(value = min(vals),    method = "min",    q = NA_real_))
  if (k == "median") return(list(value = median(vals), method = "median", q = NA_real_))
  if (k == "mean")   return(list(value = mean(vals),   method = "mean",   q = NA_real_))
  if (k == "percentile") {
    q <- spec_parsed$q
    if (!is.finite(q) || q < 0 || q > 1) stop("Percentile/quantile q must be in [0,1]. Got: ", q)
    val <- as.numeric(stats::quantile(vals, probs = q, na.rm = TRUE, type = 7))
    return(list(value = val, method = "percentile", q = q))
  }
  stop("Unhandled REF_SPEC kind: ", k)
}

# -------------------- read & basic checks ------------------------------------
x <- read.table(aligment_results, header = TRUE, sep = "\t", check.names = FALSE)
# Expecting a fixed column order from upstream; rename to stable names
names(x) <- c("Project","Sample","Raw_Reads","TrimmedReads","MouseAligned","EcoliAligned",
              "TotalAligned","MouseClean","EcoliClean","TotalClean","MouseAlignRate","PortionOfEcoli")

d <- x[, c("Sample","TrimmedReads","MouseClean","EcoliClean","TotalClean")]
d$is_IgG <- grepl("IgG", d$Sample, ignore.case = TRUE)

d$MouseClean <- as.numeric(d$MouseClean)
d$EcoliClean <- as.numeric(d$EcoliClean)

if (anyNA(d$MouseClean) || anyNA(d$EcoliClean)) {
  stop("MouseClean/EcoliClean has NA; cannot compute spike-in normalization.")
}
if (any(d$MouseClean <= 0)) {
  stop("MouseClean has zero or negative values; cannot compute sRPM.")
}

# -------------------- compute spike-in metric --------------------------------
# Spike-in per million mouse reads (stable to total depth)
d$spike_sRPM <- (d$EcoliClean / d$MouseClean) * 1e6

# -------------------- outlier detection among non-IgG only -------------------
nonIgG_idx <- which(!d$is_IgG & is.finite(d$spike_sRPM))
if (length(nonIgG_idx) < MIN_REF_N) {
  stop("Not enough non-IgG samples to define reference (need >= ", MIN_REF_N, ").")
}

is_out_nonIgG <- rep(FALSE, nrow(d))
if (length(nonIgG_idx) >= 2) {
  is_out_nonIgG[nonIgG_idx] <- detect_outliers(
    d$spike_sRPM[nonIgG_idx],
    method = OUTLIER_METHOD,
    mad_k  = MAD_K,
    iqr_k  = IQR_K
  )
}

d$is_outlier <- is_out_nonIgG  # only meaningful for non-IgG; IgG remain FALSE
n_out   <- sum(d$is_outlier & !d$is_IgG)
n_nonIgG <- sum(!d$is_IgG)

if (n_out > 0) {
  message("Outlier detection (", OUTLIER_METHOD, "): excluded ", n_out,
          " / ", n_nonIgG, " non-IgG sample(s) from reference pool.")
} else {
  message("Outlier detection (", OUTLIER_METHOD, "): no non-IgG outliers excluded.")
}

# -------------------- reference pool (non-IgG, exclude outliers) -------------
ref_candidates <- d$spike_sRPM[!d$is_IgG & !d$is_outlier & is.finite(d$spike_sRPM)]
if (length(ref_candidates) < MIN_REF_N) {
  warning("Too few non-IgG after outlier exclusion (", length(ref_candidates),
          "). Falling back to all non-IgG for ref_value.")
  ref_candidates <- d$spike_sRPM[!d$is_IgG & is.finite(d$spike_sRPM)]
}
if (length(ref_candidates) < MIN_REF_N) {
  stop("Not enough non-IgG samples to define reference even after fallback (need >= ",
       MIN_REF_N, ").")
}

# Summaries (for logs/metadata)
ref_pool_mean   <- mean(ref_candidates)
ref_pool_median <- median(ref_candidates)
ref_pool_max    <- max(ref_candidates)
ref_pool_min    <- min(ref_candidates)
ref_pool_n      <- length(ref_candidates)
ref_excluded    <- n_out

message(sprintf("Reference pool (non-IgG, post-filter): n=%d | mean=%.3f | median=%.3f | min=%.3f | max=%.3f",
                ref_pool_n, ref_pool_mean, ref_pool_median, ref_pool_min, ref_pool_max))

# -------------------- compute reference value per REF_SPEC -------------------
spec_parsed <- parse_ref_spec(REF_SPEC)                                      # NEW
ref_info    <- compute_ref_value(ref_candidates, spec_parsed)                # NEW
ref_value   <- ref_info$value
message(
  "ref_value (", ref_info$method,
  if (is.finite(ref_info$q)) paste0(", q=", ref_info$q) else "",
  ") = ", signif(ref_value, 6)
)

# -------------------- scaling + clamping -------------------------------------
d$scale_raw     <- ref_value / d$spike_sRPM
d$scale_clamped <- pmin(SF_MAX, pmax(SF_MIN, d$scale_raw))
d$clamped       <- d$scale_clamped != d$scale_raw

# -------------------- extras & metadata --------------------------------------
denom <- d$EcoliClean + d$MouseClean
d$percent_spikein        <- ifelse(denom > 0, d$EcoliClean / denom, NA_real_)
d$reference_spike_sRPM   <- ref_value
d$outlier_method         <- OUTLIER_METHOD

# Provenance / reference spec
d$ref_stat_method        <- ifelse(isTRUE(ref_info$method == "percentile"), "percentile", ref_info$method)  # NEW
d$ref_quantile           <- ifelse(is.finite(ref_info$q), ref_info$q, NA_real_)                              # NEW

# Reference pool summaries (helpful downstream)
d$ref_pool_mean          <- ref_pool_mean
d$ref_pool_median        <- ref_pool_median
d$ref_pool_min           <- ref_pool_min
d$ref_pool_max           <- ref_pool_max
d$ref_pool_n             <- ref_pool_n
d$ref_excluded_n         <- ref_excluded

d$notes <- ifelse(d$is_outlier & !d$is_IgG, "excluded_from_ref",
           ifelse(d$clamped, "clamped", "ok"))

# -------------------- output (original + new metadata) -----------------------
out <- d[, c("Sample","MouseClean","EcoliClean","TotalClean",
             "spike_sRPM","percent_spikein",
             "scale_clamped","is_outlier","notes",
             # metadata:
             "ref_stat_method","ref_quantile","reference_spike_sRPM",
             "ref_pool_mean","ref_pool_median","ref_pool_min","ref_pool_max",
             "ref_pool_n","ref_excluded_n","outlier_method")]

write.table(out, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote: ", out_file)
