# ==============================================================================
# Step 4B.2: TAC102 Flanking tBLASTn HSP & Boundary Characterization
# Project: Genomic & Transcriptomic Limits on kDNA Retention in T. equiperdum
# File: raw_read_analysis/scripts/04_tac102_flanking_hsps.R
# ==============================================================================

cat("======================================================\n")
cat(" Step 4B.2: TAC102 Flanking HSP & Boundary Analysis   \n")
cat("======================================================\n\n")

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("CONFIG")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

log_info("Starting Step 4B.2: TAC102 flanking HSP & boundary characterization.", LOG_FILE)

# ------------------------------------------------------------------------------
# 1. Load Inputs
# ------------------------------------------------------------------------------
per_read_path <- file.path(PATHS$reports, "raw_read_recovery_per_read.csv")
if (!file.exists(per_read_path)) {
  log_error(sprintf("Missing Step 3 output: %s", per_read_path), LOG_FILE)
  stop("Run 03_raw_read_recovery.R before Step 4B.2.")
}

per_read_df <- read.csv(per_read_path, stringsAsFactors = FALSE)
tac <- per_read_df[per_read_df$gene == "TAC102", , drop = FALSE]

tac$qstart_union <- as.character(tac$qstart_union)
tac$qend_union   <- as.character(tac$qend_union)

# Gap and flank definition
gap_start   <- 649L
gap_end     <- 732L
n_flank_min <- 550L
n_flank_max <- 648L
c_flank_min <- 733L
c_flank_max <- 830L

# ------------------------------------------------------------------------------
# 2. Extract All HSP Intervals & Per-Read Support
# ------------------------------------------------------------------------------
all_starts <- integer(0)
all_ends   <- integer(0)

n_flank_reads    <- character(0)
c_flank_reads    <- character(0)
both_flank_pairs <- character(0)

if (nrow(tac) > 0) {
  read_pairs <- split(tac, tac$read_id)
  
  for (r_id in names(read_pairs)) {
    sub_df <- read_pairs[[r_id]]
    has_n <- FALSE
    has_c <- FALSE
    
    for (i in seq_len(nrow(sub_df))) {
      starts <- as.integer(base::strsplit(sub_df$qstart_union[i], ";", fixed = TRUE)[[1]])
      ends   <- as.integer(base::strsplit(sub_df$qend_union[i], ";", fixed = TRUE)[[1]])
      
      all_starts <- c(all_starts, starts)
      all_ends   <- c(all_ends, ends)
      
      if (any(starts <= n_flank_max & ends >= n_flank_min)) has_n <- TRUE
      if (any(starts <= c_flank_max & ends >= c_flank_min)) has_c <- TRUE
    }
    
    if (has_n) n_flank_reads <- c(n_flank_reads, r_id)
    if (has_c) c_flank_reads <- c(c_flank_reads, r_id)
    if (has_n && has_c) both_flank_pairs <- c(both_flank_pairs, r_id)
  }
}

# ------------------------------------------------------------------------------
# 2B. Check for HSPs that directly enter the zero-support interval
# ------------------------------------------------------------------------------
gap_overlapping_hsps <- data.frame(
  start_aa = integer(0),
  end_aa   = integer(0)
)

if (length(all_starts) > 0) {
  overlaps <- all_starts <= gap_end & all_ends >= gap_start

  if (any(overlaps)) {
    gap_overlapping_hsps <- data.frame(
      start_aa = all_starts[overlaps],
      end_aa   = all_ends[overlaps]
    )
  }
}

# ------------------------------------------------------------------------------
# 3. Calculate Exact Boundary Extremes
# ------------------------------------------------------------------------------
# Maximum qend for HSPs starting at or before the N-flank boundary (<= 648)
n_side_hsps <- all_ends[all_starts <= n_flank_max]
max_n_qend  <- if (length(n_side_hsps) > 0) max(n_side_hsps) else NA_integer_

# Minimum qstart for HSPs ending at or after the C-flank boundary (>= 733)
c_side_hsps  <- all_starts[all_ends >= c_flank_min]
min_c_qstart <- if (length(c_side_hsps) > 0) min(c_side_hsps) else NA_integer_

# ------------------------------------------------------------------------------
# 4. Generate Report
# ------------------------------------------------------------------------------
report <- c(
  "==============================================",
  "TAC102 FLANKING HSP & BOUNDARY REPORT",
  "==============================================",
  "",
  sprintf("Date: %s", Sys.Date()),
  sprintf("Gene: TAC102 (reference: TREU927, 951 aa)"),
  sprintf("Gap interval: residues %d-%d aa (84 aa)", gap_start, gap_end),
  "",
  "-- Direct HSP Overlap With Zero-Support Interval --",
  sprintf("HSPs overlapping residues %d-%d: %d", gap_start, gap_end, nrow(gap_overlapping_hsps)),
  if (nrow(gap_overlapping_hsps) > 0) {
    c(
      "WARNING: Direct HSP overlap detected:",
      sprintf("  %d-%d aa", gap_overlapping_hsps$start_aa, gap_overlapping_hsps$end_aa)
    )
  } else {
    "No tBLASTn HSP directly overlaps the zero-support interval."
  },
  "",
  "-- Boundary Precision Metrics --",
  sprintf("N-flank maximum qend (HSPs starting <= %d): %s aa", 
          n_flank_max, ifelse(is.na(max_n_qend), "None", as.character(max_n_qend))),
  sprintf("C-flank minimum qstart (HSPs ending >= %d): %s aa", 
          c_flank_min, ifelse(is.na(min_c_qstart), "None", as.character(min_c_qstart))),
  sprintf("Boundary precision gap: %d-%d aa (%s aa missing)",
          ifelse(is.na(max_n_qend), gap_start, max_n_qend + 1L),
          ifelse(is.na(min_c_qstart), gap_end, min_c_qstart - 1L),
          ifelse(is.na(max_n_qend) || is.na(min_c_qstart), "NA",
                 as.character(min_c_qstart - max_n_qend - 1L))),
  "",
  "-- Flanking Read-Support Counts --",
  sprintf("N-flank zone (%d-%d aa) unique read pairs: %d", n_flank_min, n_flank_max, length(n_flank_reads)),
  sprintf("C-flank zone (%d-%d aa) unique read pairs: %d", c_flank_min, c_flank_max, length(c_flank_reads)),
  sprintf("Read pairs with support on both flanks:    %d", length(both_flank_pairs)),
  "",
  "NOTE ON TERMINOLOGY:",
  "'Read pairs with support on both flanks' indicates that the read-pair ID",
  "contains HSP hits in both flanking windows. Physical gap bridging requires",
  "direct nucleotide mapping and insert-size validation in Step 4B.3.",
  "=============================================="
)

report_out <- file.path(PATHS$reports, "tac102_flanking_hsps_report.txt")
writeLines(report, report_out)
log_info(sprintf("Saved flanking HSP report: %s", basename(report_out)), LOG_FILE)

# ------------------------------------------------------------------------------
# 5. Exit
# ------------------------------------------------------------------------------
log_info("Step 4B.2 TAC102 flanking HSP characterization completed.", LOG_FILE)
file.copy(LOG_FILE, file.path(PATHS$logs, "latest_execution.log"), overwrite = TRUE)

cat("\n[SUCCESS] Step 4B.2 completed.\n")
cat(sprintf("  Report: %s\n", report_out))
cat(sprintf("  Date: %s\n", Sys.Date()))
cat(sprintf("  Gene: TAC102 (reference: TREU927, 951 aa)\n"))