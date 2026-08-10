# ==============================================================================
# Step 4B.3a: TAC102 41-Pair Nucleotide & Query Geometry Characterization
# Project: Genomic & Transcriptomic Limits on kDNA Retention in T. equiperdum
# File: raw_read_analysis/scripts/04_tac102_nucleotide_geometry.R
# ==============================================================================

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("CONFIG")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

log_info("Starting Step 4B.3a: Extracting geometry for 41 both-flank pairs.", LOG_FILE)

# 1. Load data & extract TAC102 records
per_read_df <- read.csv(file.path(PATHS$reports, "raw_read_recovery_per_read.csv"), stringsAsFactors = FALSE)
tac <- per_read_df[per_read_df$gene == "TAC102", , drop = FALSE]

# Force character type on all coordinate union columns
tac$qstart_union <- as.character(tac$qstart_union)
tac$qend_union   <- as.character(tac$qend_union)
tac$sstart_union <- as.character(tac$sstart_union)
tac$send_union   <- as.character(tac$send_union)

n_flank_min <- 550L; n_flank_max <- 648L
c_flank_min <- 733L; c_flank_max <- 830L

# 2. Identify pair IDs with hits in both flanking windows
pair_summary <- do.call(rbind, lapply(split(tac, tac$pair_id), function(sub_df) {
  has_n <- FALSE
  has_c <- FALSE
  for (i in seq_len(nrow(sub_df))) {
    starts <- as.integer(base::strsplit(as.character(sub_df$qstart_union[i]), ";", fixed = TRUE)[[1]])
    ends   <- as.integer(base::strsplit(as.character(sub_df$qend_union[i]), ";", fixed = TRUE)[[1]])
    if (any(starts <= n_flank_max & ends >= n_flank_min)) has_n <- TRUE
    if (any(starts <= c_flank_max & ends >= c_flank_min)) has_c <- TRUE
  }
  data.frame(pair_id = sub_df$pair_id[1], has_N = has_n, has_C = has_c, stringsAsFactors = FALSE)
}))

both_flank_ids <- pair_summary$pair_id[pair_summary$has_N & pair_summary$has_C]
flank_records  <- tac[tac$pair_id %in% both_flank_ids, ]

# 3. Extract exact coordinate bounds and calculate gap span
nucleotide_pair_summary <- do.call(rbind, lapply(split(flank_records, flank_records$pair_id), function(df) {
  r1 <- df[df$mate == "R1", , drop = FALSE]
  r2 <- df[df$mate == "R2", , drop = FALSE]

  r1_qs <- if (nrow(r1) > 0) as.integer(base::strsplit(as.character(r1$qstart_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_
  r1_qe <- if (nrow(r1) > 0) as.integer(base::strsplit(as.character(r1$qend_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_
  r1_ss <- if (nrow(r1) > 0) as.integer(base::strsplit(as.character(r1$sstart_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_
  r1_se <- if (nrow(r1) > 0) as.integer(base::strsplit(as.character(r1$send_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_

  r2_qs <- if (nrow(r2) > 0) as.integer(base::strsplit(as.character(r2$qstart_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_
  r2_qe <- if (nrow(r2) > 0) as.integer(base::strsplit(as.character(r2$qend_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_
  r2_ss <- if (nrow(r2) > 0) as.integer(base::strsplit(as.character(r2$sstart_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_
  r2_se <- if (nrow(r2) > 0) as.integer(base::strsplit(as.character(r2$send_union[1]), ";", fixed = TRUE)[[1]][1]) else NA_integer_

  # Query gap distance (aa) between R2_qend (N-side) and R1_qstart (C-side)
  q_gap_aa <- if (!is.na(r1_qs) && !is.na(r2_qe)) r1_qs - r2_qe - 1L else NA_integer_

  data.frame(
    pair_id      = df$pair_id[1],
    r2_qstart    = r2_qs,
    r2_qend      = r2_qe,
    r2_sstart    = r2_ss,
    r2_send      = r2_se,
    r1_qstart    = r1_qs,
    r1_qend      = r1_qe,
    r1_sstart    = r1_ss,
    r1_send      = r1_se,
    query_gap_aa = q_gap_aa,
    stringsAsFactors = FALSE
  )
}))

# 4. Save CSV report
csv_out <- file.path(PATHS$reports, "tac102_41pairs_geometry.csv")
write.csv(nucleotide_pair_summary, csv_out, row.names = FALSE)
log_info(sprintf("Saved 41-pair geometry report to %s", csv_out), LOG_FILE)

# 5. Output summary
cat("\n======================================================\n")
cat(sprintf(" TAC102 Geometry Extraction Complete (%d pairs)\n", nrow(nucleotide_pair_summary)))
cat("======================================================\n\n")
print(head(nucleotide_pair_summary, 10))