# ==============================================================================
# Step 4B.3b PRE-FLIGHT:
# TAC102 Raw HSP Structure and Identifier Validation
#
# Project: Genomic & Transcriptomic Limits on kDNA Retention in T. equiperdum
#
# PURPOSE
# -------
# Before performing codon-anchored extension into TAC102 residues 649–732,
# validate:
#
#   1. The structure/columns of the raw HSP table
#   2. The identifier relationship between pair_id/read_id and sseqid
#   3. The raw HSP coordinates for representative both-flank pairs
#   4. Whether the original HSPs appear suitable for +3/-3 codon extrapolation
#   5. Whether any representative HSP contains evidence of an alignment gap
#
# IMPORTANT
# ---------
# This script DOES NOT classify A/B/C/D.
# It DOES NOT extract gap-facing sequence.
# It DOES NOT alter the frozen protocol.
#
# It is a pre-flight validation step before Step 4B.3b.
# ==============================================================================


cat("==============================================================\n")
cat(" Step 4B.3b PRE-FLIGHT: TAC102 Raw HSP Validation            \n")
cat("==============================================================\n\n")


# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("CONFIG")) {
  source(here::here(
    "raw_read_analysis",
    "scripts",
    "00_setup.R"
  ))
}

log_info(
  "Starting TAC102 Step 4B.3b pre-flight HSP validation.",
  LOG_FILE
)


# ------------------------------------------------------------------------------
# 1. Load per-read summary and raw HSP table
# ------------------------------------------------------------------------------

per_read_file <- file.path(
  PATHS$reports,
  "raw_read_recovery_per_read.csv"
)

all_hsps_file <- file.path(
  PATHS$reports,
  "raw_read_recovery_all_hsps.csv"
)

per_read_df <- read.csv(
  per_read_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all_hsps <- read.csv(
  all_hsps_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

tac_pr <- per_read_df[
  per_read_df$gene == "TAC102",
  ,
  drop = FALSE
]

tac_hsps <- all_hsps[
  all_hsps$gene == "TAC102",
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# 2. Report available columns
# ------------------------------------------------------------------------------

cat("\n--------------------------------------------------------------\n")
cat("1. TAC102 RAW HSP COLUMN STRUCTURE\n")
cat("--------------------------------------------------------------\n\n")

cat("Number of TAC102 raw HSP rows:", nrow(tac_hsps), "\n\n")

cat("TAC102 raw HSP columns:\n")
print(names(tac_hsps))

cat("\n\nPer-read TAC102 columns:\n")
print(names(tac_pr))


# ------------------------------------------------------------------------------
# 3. Check expected raw HSP columns
# ------------------------------------------------------------------------------

expected_hsp_columns <- c(
  "gene",
  "mate",
  "sseqid",
  "qstart",
  "qend",
  "sstart",
  "send",
  "pident",
  "evalue",
  "bitscore"
)

missing_hsp_columns <- setdiff(
  expected_hsp_columns,
  names(tac_hsps)
)

cat("\n--------------------------------------------------------------\n")
cat("2. EXPECTED HSP COLUMN CHECK\n")
cat("--------------------------------------------------------------\n\n")

if (length(missing_hsp_columns) == 0) {

  cat("PASS: All expected core HSP columns are present.\n")

} else {

  cat("WARNING: Missing expected columns:\n")
  print(missing_hsp_columns)
}


# ------------------------------------------------------------------------------
# 4. Inspect data types of coordinate columns
# ------------------------------------------------------------------------------

cat("\n--------------------------------------------------------------\n")
cat("3. HSP COORDINATE DATA TYPES\n")
cat("--------------------------------------------------------------\n\n")

coordinate_columns <- intersect(
  c("qstart", "qend", "sstart", "send"),
  names(tac_hsps)
)

for (col in coordinate_columns) {

  cat(
    sprintf(
      "%-10s : %s\n",
      col,
      paste(class(tac_hsps[[col]]), collapse = ", ")
    )
  )
}


# ------------------------------------------------------------------------------
# 5. Check identifier conventions
# ------------------------------------------------------------------------------

cat("\n--------------------------------------------------------------\n")
cat("4. IDENTIFIER CONVENTION CHECK\n")
cat("--------------------------------------------------------------\n\n")

pair_ids <- unique(tac_pr$pair_id)
read_ids <- unique(tac_pr$read_id)

cat("Unique TAC102 pair IDs :", length(pair_ids), "\n")
cat("Unique TAC102 read IDs :", length(read_ids), "\n")
cat("Unique TAC102 sseqids  :", length(unique(tac_hsps$sseqid)), "\n\n")

cat("Example pair_id values:\n")
print(head(pair_ids, 10))

cat("\nExample read_id values:\n")
print(head(read_ids, 10))

cat("\nExample sseqid values:\n")
print(head(unique(tac_hsps$sseqid), 10))


# ------------------------------------------------------------------------------
# 6. Explicit identifier matching test
# ------------------------------------------------------------------------------

pair_in_sseqid <- pair_ids %in% tac_hsps$sseqid
read_in_sseqid <- read_ids %in% tac_hsps$sseqid

cat("\nPair IDs found directly in sseqid:\n")
cat(
  sprintf(
    "%d / %d\n",
    sum(pair_in_sseqid),
    length(pair_in_sseqid)
  )
)

cat("\nRead IDs found directly in sseqid:\n")
cat(
  sprintf(
    "%d / %d\n",
    sum(read_in_sseqid),
    length(read_in_sseqid)
  )
)


# ------------------------------------------------------------------------------
# 7. Representative pairs
#
# Includes:
#   10334997       clean 84-aa case
#   10897477       clean 84-aa case
#   21299          mate-orientation flip
#   26367709       mate-orientation flip
# ------------------------------------------------------------------------------

representative_ids <- c(
  "SRR7910035.10334997",
  "SRR7910035.10897477",
  "SRR7910035.21299",
  "SRR7910035.26367709"
)

representative_ids <- representative_ids[
  representative_ids %in% pair_ids
]

cat("\n--------------------------------------------------------------\n")
cat("5. REPRESENTATIVE PAIRS\n")
cat("--------------------------------------------------------------\n\n")

print(representative_ids)


# ------------------------------------------------------------------------------
# 8. Display ALL raw HSP records for representatives
# ------------------------------------------------------------------------------

for (pid in representative_ids) {

  cat("\n==============================================================\n")
  cat("PAIR:", pid, "\n")
  cat("==============================================================\n\n")

  sub <- tac_hsps[
    tac_hsps$sseqid == pid,
    ,
    drop = FALSE
  ]

  if (nrow(sub) == 0) {

    cat("NO RAW HSP RECORD FOUND USING DIRECT sseqid MATCH.\n")

    # Try read_id as fallback diagnostic only
    fallback <- tac_hsps[
      tac_hsps$sseqid %in% unique(
        tac_pr$read_id[tac_pr$pair_id == pid]
      ),
      ,
      drop = FALSE
    ]

    if (nrow(fallback) > 0) {

      cat("\nFallback read_id match FOUND:\n")
      print(fallback)

    } else {

      cat("No fallback read_id match found either.\n")
    }

  } else {

    cat("Raw HSP records found:", nrow(sub), "\n\n")

    print(sub)
  }
}


# ------------------------------------------------------------------------------
# 9. Check orientation of representative HSPs
# ------------------------------------------------------------------------------

cat("\n--------------------------------------------------------------\n")
cat("6. REPRESENTATIVE HSP STRAND ORIENTATION\n")
cat("--------------------------------------------------------------\n\n")

orientation_summary <- do.call(
  rbind,
  lapply(
    representative_ids,
    function(pid) {

      sub <- tac_hsps[
        tac_hsps$sseqid == pid,
        ,
        drop = FALSE
      ]

      if (nrow(sub) == 0) {
        return(
          data.frame(
            pair_id = pid,
            n_hsps = 0L,
            forward_hsps = NA_integer_,
            reverse_hsps = NA_integer_,
            stringsAsFactors = FALSE
          )
        )
      }

      forward <- sub$sstart < sub$send
      reverse <- sub$sstart > sub$send

      data.frame(
        pair_id = pid,
        n_hsps = nrow(sub),
        forward_hsps = sum(forward, na.rm = TRUE),
        reverse_hsps = sum(reverse, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )
)

print(orientation_summary)


# ------------------------------------------------------------------------------
# 10. Check whether qstart/qend and sstart/send imply a simple codon mapping
#
# For a gap-free protein-to-nucleotide alignment:
#
#       subject_nt_span = protein_span * 3
#
# where:
#
#       protein_span = abs(qend - qstart) + 1
#       nucleotide_span = abs(send - sstart) + 1
#
# If this equality holds, the HSP is consistent with a continuous
# codon-by-codon relationship at the coordinate level.
#
# NOTE:
# This does NOT prove absence of gaps in the actual BLAST alignment.
# It is a diagnostic consistency check only.
# ------------------------------------------------------------------------------

cat("\n--------------------------------------------------------------\n")
cat("7. CODON-SPAN CONSISTENCY CHECK\n")
cat("--------------------------------------------------------------\n\n")

if (
  all(
    c("qstart", "qend", "sstart", "send") %in%
      names(tac_hsps)
  )
) {

  tac_hsps$q_span_aa <-
    abs(tac_hsps$qend - tac_hsps$qstart) + 1L

  tac_hsps$s_span_nt <-
    abs(tac_hsps$send - tac_hsps$sstart) + 1L

  tac_hsps$expected_nt_span <-
    tac_hsps$q_span_aa * 3L

  tac_hsps$span_matches_3x <-
    tac_hsps$s_span_nt == tac_hsps$expected_nt_span

  cat(
    "Total TAC102 HSPs:",
    nrow(tac_hsps),
    "\n"
  )

  cat(
    "HSPs with nucleotide span exactly 3x protein span:",
    sum(tac_hsps$span_matches_3x, na.rm = TRUE),
    "\n"
  )

  cat(
    "HSPs failing the 3x span check:",
    sum(!tac_hsps$span_matches_3x, na.rm = TRUE),
    "\n"
  )

  cat("\nDistribution:\n")
  print(
    table(
      factor(
        tac_hsps$span_matches_3x,
        levels = c(TRUE, FALSE)
      )
    )
  )

} else {

  cat(
    "Cannot perform span check because required coordinate columns are missing.\n"
  )
}


# ------------------------------------------------------------------------------
# 11. Show HSPs that fail the simple 3x codon-span relationship
# ------------------------------------------------------------------------------

if ("span_matches_3x" %in% names(tac_hsps)) {

  failed_span <- tac_hsps[
    !tac_hsps$span_matches_3x,
    ,
    drop = FALSE
  ]

  cat("\n--------------------------------------------------------------\n")
  cat("8. HSPs FAILING 3x CODON-SPAN CHECK\n")
  cat("--------------------------------------------------------------\n\n")

  cat(
    "Number failing:",
    nrow(failed_span),
    "\n\n"
  )

  if (nrow(failed_span) > 0) {

    display_cols <- intersect(
      c(
        "gene",
        "mate",
        "sseqid",
        "qstart",
        "qend",
        "q_span_aa",
        "sstart",
        "send",
        "s_span_nt",
        "expected_nt_span",
        "pident",
        "evalue",
        "bitscore"
      ),
      names(failed_span)
    )

    print(
      failed_span[, display_cols, drop = FALSE]
    )
  }
}


# ------------------------------------------------------------------------------
# 12. Check for suspicious qstart/qend ordering
# ------------------------------------------------------------------------------

cat("\n--------------------------------------------------------------\n")
cat("9. QUERY COORDINATE SANITY CHECK\n")
cat("--------------------------------------------------------------\n\n")

bad_q_order <- tac_hsps[
  tac_hsps$qstart > tac_hsps$qend,
  ,
  drop = FALSE
]

cat(
  "HSPs with qstart > qend:",
  nrow(bad_q_order),
  "\n"
)

if (nrow(bad_q_order) > 0) {
  print(bad_q_order)
}


# ------------------------------------------------------------------------------
# 13. Save diagnostic table
# ------------------------------------------------------------------------------

diagnostic_file <- file.path(
  PATHS$reports,
  "tac102_hsp_preflight_diagnostics.csv"
)

if (
  all(
    c(
      "q_span_aa",
      "s_span_nt",
      "expected_nt_span",
      "span_matches_3x"
    ) %in%
    names(tac_hsps)
  )
) {

  write.csv(
    tac_hsps,
    diagnostic_file,
    row.names = FALSE
  )

  cat(
    "\nDiagnostic HSP table saved to:\n",
    diagnostic_file,
    "\n"
  )
}


# ------------------------------------------------------------------------------
# 14. Final interpretation
# ------------------------------------------------------------------------------

cat("\n==============================================================\n")
cat(" PRE-FLIGHT COMPLETE\n")
cat("==============================================================\n\n")

cat(
  "IMPORTANT:\n",
  "This script does NOT perform Step 4B.3b classification.\n",
  "It only validates the raw HSP coordinate structure needed for\n",
  "codon-anchored extension.\n\n"
)

cat(
  "Next decision depends on:\n",
  "  1. Whether sseqid directly matches pair/read IDs\n",
  "  2. Whether the representative HSPs are structurally consistent\n",
  "  3. Whether the 3x codon-span diagnostic is satisfied\n",
  "  4. Whether the raw HSP representation contains enough information\n",
  "     to guarantee frame-preserving extension\n"
)

log_info(
  "TAC102 Step 4B.3b pre-flight validation completed.",
  LOG_FILE
)