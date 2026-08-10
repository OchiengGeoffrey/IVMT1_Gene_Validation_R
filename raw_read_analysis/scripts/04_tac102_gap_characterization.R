# ==============================================================================
# Step 4A: TAC102 Raw-Read Support Gap Characterization
# Project: Genomic & Transcriptomic Limits on kDNA Retention in T. equiperdum
#
# Objective:
#   Identify the exact TAC102 protein residues with zero tBLASTn support
#   in the IVM-t1 raw reads and summarize the resulting continuous gaps.
#
# Important:
#   "Support depth" refers to protein-homology support from tBLASTn.
#   It is NOT nucleotide sequencing depth.
#
# Input:
#   raw_read_recovery_per_read.csv
#
# Output:
#   tac102_support_depth_by_residue.csv
#   tac102_uncovered_intervals.csv
#   tac102_gap_report.txt
# ==============================================================================

cat("======================================================\n")
cat(" Step 4A: TAC102 Gap Characterization                 \n")
cat("======================================================\n\n")

# ------------------------------------------------------------------------------
# 1. Environment
# ------------------------------------------------------------------------------

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("CONFIG")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

log_info("Starting Step 4A: TAC102 gap characterization.", LOG_FILE)

# ------------------------------------------------------------------------------
# 2. Input
# ------------------------------------------------------------------------------

PER_READ <- file.path(
  PATHS$reports,
  "raw_read_recovery_per_read.csv"
)

if (!file.exists(PER_READ)) {
  log_error(
    sprintf("Missing Step 3 output: %s", PER_READ),
    LOG_FILE
  )
  stop("Run 03_raw_read_recovery.R before Step 4A.")
}

per_read_df <- read.csv(
  PER_READ,
  stringsAsFactors = FALSE
)

gene <- "TAC102"
protein_len <- 951L

tac <- per_read_df[
  per_read_df$gene == gene,
  ,
  drop = FALSE
]

# Force character type to prevent numeric type inference by read.csv
tac$qstart_union <- as.character(tac$qstart_union)
tac$qend_union   <- as.character(tac$qend_union)

log_info(
  sprintf(
    "Loaded %d TAC102 read-level evidence rows.",
    nrow(tac)
  ),
  LOG_FILE
)

# ------------------------------------------------------------------------------
# 3. Validate required columns
# ------------------------------------------------------------------------------

required_cols <- c(
  "gene",
  "read_id",
  "mate",
  "qstart_union",
  "qend_union"
)

missing_cols <- setdiff(required_cols, names(tac))

if (length(missing_cols) > 0) {
  stop(
    sprintf(
      "Missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------------------------
# 4. Calculate protein-position support depth
# ------------------------------------------------------------------------------

support_depth <- integer(protein_len)

if (nrow(tac) > 0) {

  for (i in seq_len(nrow(tac))) {

    starts <- as.integer(
      base::strsplit(
        as.character(tac$qstart_union[i]),
        ";",
        fixed = TRUE
      )[[1]]
    )

    ends <- as.integer(
      base::strsplit(
        as.character(tac$qend_union[i]),
        ";",
        fixed = TRUE
      )[[1]]
    )

    if (length(starts) != length(ends)) {
      stop(
        sprintf(
          "Coordinate mismatch in TAC102 row %d.",
          i
        )
      )
    }

    for (k in seq_along(starts)) {

      s <- max(1L, starts[k])
      e <- min(protein_len, ends[k])

      if (s <= e) {
        support_depth[s:e] <- support_depth[s:e] + 1L
      }
    }
  }
}

# ------------------------------------------------------------------------------
# 5. Per-residue table
# ------------------------------------------------------------------------------

residue_df <- data.frame(
  gene = gene,
  protein_position = seq_len(protein_len),
  support_depth = support_depth,
  supported = support_depth > 0,
  stringsAsFactors = FALSE
)

residue_out <- file.path(
  PATHS$reports,
  "tac102_support_depth_by_residue.csv"
)

write.csv(
  residue_df,
  residue_out,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Identify uncovered residues
# ------------------------------------------------------------------------------

uncovered <- which(support_depth == 0)

coverage_pct <- mean(support_depth > 0) * 100
uncovered_pct <- mean(support_depth == 0) * 100

# ------------------------------------------------------------------------------
# 7. Convert uncovered residues into continuous intervals
# ------------------------------------------------------------------------------

if (length(uncovered) > 0) {

  starts <- uncovered[
    c(TRUE, diff(uncovered) != 1)
  ]

  ends <- uncovered[
    c(diff(uncovered) != 1, TRUE)
  ]

  gap_df <- data.frame(
    gene = gene,
    start_aa = starts,
    end_aa = ends,
    length_aa = ends - starts + 1L,
    stringsAsFactors = FALSE
  )

} else {

  gap_df <- data.frame(
    gene = character(),
    start_aa = integer(),
    end_aa = integer(),
    length_aa = integer(),
    stringsAsFactors = FALSE
  )
}

gap_out <- file.path(
  PATHS$reports,
  "tac102_uncovered_intervals.csv"
)

write.csv(
  gap_df,
  gap_out,
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Gap report
# ------------------------------------------------------------------------------

report <- c(
  "==============================================",
  "TAC102 RAW-READ SUPPORT GAP REPORT",
  "==============================================",
  "",
  sprintf("Date: %s", Sys.Date()),
  sprintf("Gene: %s", gene),
  sprintf("Protein length: %d aa", protein_len),
  "",
  "Evidence basis:",
  "Protein-vs-translated-read tBLASTn evidence from Step 3.",
  "Support depth represents the number of distinct read-level",
  "protein-homology observations covering each reference residue.",
  "",
  sprintf("Supported residues: %d / %d (%.2f%%)",
          sum(support_depth > 0),
          protein_len,
          coverage_pct),
  sprintf("Uncovered residues: %d / %d (%.2f%%)",
          sum(support_depth == 0),
          protein_len,
          uncovered_pct),
  "",
  "Continuous zero-support intervals:"
)

if (nrow(gap_df) > 0) {

  for (i in seq_len(nrow(gap_df))) {

    report <- c(
      report,
      sprintf(
        "  %d-%d aa (%d aa)",
        gap_df$start_aa[i],
        gap_df$end_aa[i],
        gap_df$length_aa[i]
      )
    )
  }

} else {

  report <- c(
    report,
    "  NONE"
  )
}

report <- c(
  report,
  "",
  "IMPORTANT:",
  "Zero tBLASTn support does not by itself establish gene absence.",
  "Possible explanations include sequence divergence, low complexity,",
  "repetitive sequence, read/assembly limitations, or true absence.",
  "",
  "The uncovered intervals should therefore be compared against:",
  "1. the Phase 1B assembly-based TAC102 recovery, and",
  "2. the underlying raw-read sequence evidence.",
  "",
  "=============================================="
)

report_out <- file.path(
  PATHS$reports,
  "tac102_gap_report.txt"
)

writeLines(report, report_out)

# ------------------------------------------------------------------------------
# 9. Logging
# ------------------------------------------------------------------------------

log_info(
  sprintf(
    "TAC102 support coverage: %.2f%%; uncovered: %.2f%%.",
    coverage_pct,
    uncovered_pct
  ),
  LOG_FILE
)

log_info(
  sprintf(
    "Identified %d continuous zero-support interval(s).",
    nrow(gap_df)
  ),
  LOG_FILE
)

log_info(
  "Step 4A TAC102 gap characterization completed.",
  LOG_FILE
)

file.copy(
  LOG_FILE,
  file.path(PATHS$logs, "latest_execution.log"),
  overwrite = TRUE
)

cat("\n[SUCCESS] Step 4A completed.\n")
cat(sprintf(
  "  Residue table: %s\n",
  residue_out
))
cat(sprintf(
  "  Gap intervals: %s\n",
  gap_out
))
cat(sprintf(
  "  Gap report: %s\n",
  report_out
))