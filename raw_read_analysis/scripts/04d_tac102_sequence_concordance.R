# ==============================================================
# Step 4B.3d: TAC102 Pair-to-Pair Sequence Concordance
# ==============================================================
#
# Purpose:
#   Evaluate amino-acid sequence concordance among the 28
#   TAC102 paired reads that have two-sided evidence spanning
#   the complete TAC102 gap region (AA 649-732).
#
# Important methodological constraints:
#   1. The TAC102 nucleotide reference is unavailable.
#   2. Therefore, this is NOT nucleotide-level validation.
#   3. The analysis uses the original-HSP-frame,
#      codon-anchored translated candidate sequences generated
#      in Step 4B.3b.
#   4. aligned_aa is a COUNT, not an amino-acid sequence.
#      Therefore it must not be split into characters.
#   5. Positional amino-acid reconstruction uses:
#        translated_candidate
#        gap_position_start
#        gap_position_end
#   6. Only sequences whose translated length is exactly
#      compatible with their reported TAC102 coordinate span
#      are used for positional concordance.
#
# ==============================================================

library(dplyr)
library(readr)
library(tidyr)

# --------------------------------------------------------------
# 1. PATH DEFINITIONS
# --------------------------------------------------------------

csv_pairs <- paste0(
  "raw_read_analysis/reports/",
  "tac102_41pairs_pair_level_evidence.csv"
)

csv_evidence <- paste0(
  "raw_read_analysis/reports/",
  "tac102_41pairs_gap_evidence.csv"
)

csv_tails <- paste0(
  "raw_read_analysis/reports/",
  "tac102_41pairs_unaligned_tails.csv"
)

output_positional <- paste0(
  "raw_read_analysis/reports/",
  "tac102_gap_positional_concordance.csv"
)

output_sequence <- paste0(
  "raw_read_analysis/reports/",
  "tac102_gap_sequence_concordance.csv"
)

output_consensus <- paste0(
  "raw_read_analysis/reports/",
  "tac102_gap_consensus_peptide.fasta"
)

output_report <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3d_sequence_concordance_report.txt"
)

# --------------------------------------------------------------
# 2. HEADER
# --------------------------------------------------------------

cat("==============================================================\n")
cat(" Step 4B.3d: TAC102 Pair-to-Pair Sequence Concordance\n")
cat("==============================================================\n\n")

cat("Target TAC102 coordinate range: AA 649 - 732 (84 AA positions)\n\n")

# --------------------------------------------------------------
# 3. CHECK INPUT FILES
# --------------------------------------------------------------

required_files <- c(
  csv_pairs,
  csv_evidence,
  csv_tails
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required input file(s):",
      paste(missing_files, collapse = "\n")
    )
  )
}

# --------------------------------------------------------------
# 4. LOAD INPUT DATA
# --------------------------------------------------------------

pairs_df <- read_csv(
  csv_pairs,
  show_col_types = FALSE
)

ev_df <- read_csv(
  csv_evidence,
  show_col_types = FALSE
)

tails_df <- read_csv(
  csv_tails,
  show_col_types = FALSE
)

cat("Input rows:\n")
cat("Pair-level:", nrow(pairs_df), "\n")
cat("Gap-evidence:", nrow(ev_df), "\n")
cat("Unaligned tails:", nrow(tails_df), "\n\n")

# --------------------------------------------------------------
# 5. VALIDATE REQUIRED COLUMNS
# --------------------------------------------------------------

required_pair_cols <- c(
  "pair_id",
  "full_gap_649_732"
)

required_ev_cols <- c(
  "pair_id",
  "mate",
  "orientation_role",
  "gap_position_start",
  "gap_position_end",
  "evidence_class"
)

required_tail_cols <- c(
  "pair_id",
  "mate",
  "orientation_role",
  "hsp_strand",
  "gap_facing_nt_length",
  "gap_facing_sequence",
  "translated_candidate",
  "gap_position_start",
  "gap_position_end",
  "aligned_aa",
  "pident_aa",
  "evidence_class"
)

missing_pair_cols <- setdiff(
  required_pair_cols,
  names(pairs_df)
)

missing_ev_cols <- setdiff(
  required_ev_cols,
  names(ev_df)
)

missing_tail_cols <- setdiff(
  required_tail_cols,
  names(tails_df)
)

if (length(missing_pair_cols) > 0) {
  stop(
    paste(
      "Missing required pair-level columns:",
      paste(missing_pair_cols, collapse = ", ")
    )
  )
}

if (length(missing_ev_cols) > 0) {
  stop(
    paste(
      "Missing required gap-evidence columns:",
      paste(missing_ev_cols, collapse = ", ")
    )
  )
}

if (length(missing_tail_cols) > 0) {
  stop(
    paste(
      "Missing required tail columns:",
      paste(missing_tail_cols, collapse = ", ")
    )
  )
}

# --------------------------------------------------------------
# 6. IDENTIFY THE 28 QUALIFYING PAIRS
# --------------------------------------------------------------

qualifying_pairs <- pairs_df %>%
  filter(full_gap_649_732 == TRUE) %>%
  pull(pair_id)

cat("Qualifying full-gap pairs:", length(qualifying_pairs), "\n")
cat("Expected:", 28, "\n\n")

if (length(qualifying_pairs) != 28) {
  warning(
    paste(
      "Expected 28 qualifying pairs, but found",
      length(qualifying_pairs)
    )
  )
}

# --------------------------------------------------------------
# 7. EXTRACT THE 28 PAIR-LEVEL RECORDS
# --------------------------------------------------------------

qual_pairs_df <- pairs_df %>%
  filter(pair_id %in% qualifying_pairs)

# --------------------------------------------------------------
# 8. EXTRACT THE CORRESPONDING TAIL SEQUENCES
# --------------------------------------------------------------

qual_tails <- tails_df %>%
  filter(pair_id %in% qualifying_pairs)

cat("Tail rows belonging to qualifying pairs:",
    nrow(qual_tails), "\n")

cat(
  "Expected approximately:",
  length(qualifying_pairs) * 2,
  "\n\n"
)

# --------------------------------------------------------------
# 9. KEEP ONLY SEQUENCES WITH A VALID GAP COORDINATE RANGE
# --------------------------------------------------------------

qual_tails <- qual_tails %>%
  mutate(
    valid_gap_coordinates =
      !is.na(gap_position_start) &
      !is.na(gap_position_end) &
      gap_position_start >= 649 &
      gap_position_end <= 732 &
      gap_position_start <= gap_position_end,

    translated_candidate =
      ifelse(
        is.na(translated_candidate),
        "",
        translated_candidate
      ),

    translated_candidate =
      toupper(translated_candidate),

    translated_length =
      nchar(translated_candidate),

    expected_coordinate_length =
      ifelse(
        valid_gap_coordinates,
        gap_position_end - gap_position_start + 1,
        NA_real_
      ),

    sequence_coordinate_length_match =
      valid_gap_coordinates &
      translated_length == expected_coordinate_length
  )

cat("=== SEQUENCE / COORDINATE VALIDATION ===\n")

cat(
  "Rows with valid TAC102 gap coordinates:",
  sum(qual_tails$valid_gap_coordinates),
  "\n"
)

cat(
  "Rows with sequence length matching coordinate span:",
  sum(qual_tails$sequence_coordinate_length_match),
  "\n\n"
)

# --------------------------------------------------------------
# 10. SHOW ANY LENGTH / COORDINATE DISCORDANCE
# --------------------------------------------------------------

length_discordance <- qual_tails %>%
  filter(
    valid_gap_coordinates &
      !sequence_coordinate_length_match
  ) %>%
  select(
    pair_id,
    mate,
    orientation_role,
    gap_position_start,
    gap_position_end,
    expected_coordinate_length,
    translated_length,
    aligned_aa,
    pident_aa,
    evidence_class
  )

if (nrow(length_discordance) > 0) {

  cat("WARNING: sequence/coordinate length discordance detected:\n")
  print(length_discordance, n = Inf)
  cat("\n")

} else {

  cat(
    "All usable translated candidates have lengths ",
    "consistent with their reported TAC102 coordinate spans.\n\n",
    sep = ""
  )
}

# --------------------------------------------------------------
# 11. KEEP POSITIONALLY USABLE SEQUENCES
# --------------------------------------------------------------

usable_tails <- qual_tails %>%
  filter(sequence_coordinate_length_match)

cat(
  "Positionally usable translated sequences:",
  nrow(usable_tails),
  "\n\n"
)

# --------------------------------------------------------------
# 12. EXPAND EACH TRANSLATED SEQUENCE INTO TAC102 POSITIONS
# --------------------------------------------------------------

positional_aa <- usable_tails %>%
  rowwise() %>%
  do({

    row_data <- .

    aa_vec <- strsplit(
      row_data$translated_candidate,
      "",
      fixed = TRUE
    )[[1]]

    pos_seq <- seq(
      from = row_data$gap_position_start,
      to   = row_data$gap_position_end
    )

    if (length(aa_vec) != length(pos_seq)) {

      tibble()

    } else {

      tibble(
        pair_id          = row_data$pair_id,
        mate             = row_data$mate,
        orientation_role = row_data$orientation_role,
        evidence_class   = row_data$evidence_class,
        pident_aa        = row_data$pident_aa,
        tac102_pos       = pos_seq,
        amino_acid       = aa_vec
      )
    }

  }) %>%
  ungroup()

cat(
  "Expanded positional amino-acid observations:",
  nrow(positional_aa),
  "\n\n"
)

# --------------------------------------------------------------
# 13. RESTRICT TO THE COMPLETE 649-732 REGION
# --------------------------------------------------------------

positional_aa <- positional_aa %>%
  filter(
    tac102_pos >= 649,
    tac102_pos <= 732
  )

# --------------------------------------------------------------
# 14. POSITIONAL CONCORDANCE
# --------------------------------------------------------------

position_summary <- positional_aa %>%
  group_by(tac102_pos) %>%
  summarise(

    coverage_depth = n(),

    unique_aa_count =
      n_distinct(amino_acid),

    dominant_aa =
      names(
        which.max(
          table(amino_acid)
        )
      ),

    dominant_freq =
      max(
        table(amino_acid)
      ),

    concordance_pct =
      round(
        dominant_freq /
          coverage_depth *
          100,
        1
      ),

    all_observed_aa =
      paste(
        sort(
          unique(amino_acid)
        ),
        collapse = "/"
      ),

    .groups = "drop"
  ) %>%
  arrange(tac102_pos)

# --------------------------------------------------------------
# 15. COMPLETE POSITION COVERAGE CHECK
# --------------------------------------------------------------

expected_positions <- tibble(
  tac102_pos = 649:732
)

position_summary_complete <- expected_positions %>%
  left_join(
    position_summary,
    by = "tac102_pos"
  )

# --------------------------------------------------------------
# 16. OVERALL CONCORDANCE METRICS
# --------------------------------------------------------------

total_positions <- 84

positions_covered <- sum(
  !is.na(position_summary_complete$coverage_depth)
)

positions_missing <- total_positions - positions_covered

if (positions_covered > 0) {

  mean_depth <- round(
    mean(
      position_summary_complete$coverage_depth,
      na.rm = TRUE
    ),
    2
  )

  high_concordance_positions <- sum(
    position_summary_complete$concordance_pct >= 90,
    na.rm = TRUE
  )

  complete_concordance_positions <- sum(
    position_summary_complete$concordance_pct == 100,
    na.rm = TRUE
  )

} else {

  mean_depth <- NA_real_
  high_concordance_positions <- 0
  complete_concordance_positions <- 0
}

cat("==============================================================\n")
cat(" CONCORDANCE METRICS\n")
cat("==============================================================\n\n")

cat(
  "Target positions:",
  total_positions,
  "\n"
)

cat(
  "Positions covered:",
  positions_covered,
  "\n"
)

cat(
  "Positions missing:",
  positions_missing,
  "\n"
)

cat(
  "Mean sequence depth per covered AA position:",
  mean_depth,
  "\n"
)

cat(
  "Positions with >=90% AA concordance:",
  high_concordance_positions,
  "\n"
)

cat(
  "Positions with 100% AA concordance:",
  complete_concordance_positions,
  "\n\n"
)

# --------------------------------------------------------------
# 17. RECONSTRUCT CONSENSUS PEPTIDE
# --------------------------------------------------------------

consensus_peptide <- paste(
  position_summary_complete$dominant_aa[
    !is.na(position_summary_complete$dominant_aa)
  ],
  collapse = ""
)

cat("==============================================================\n")
cat(" RECONSTRUCTED TAC102 GAP CONSENSUS PEPTIDE\n")
cat("==============================================================\n\n")

cat(
  "Coordinates: 649-732\n"
)

cat(
  "Consensus length:",
  nchar(consensus_peptide),
  "aa\n\n"
)

cat(
  consensus_peptide,
  "\n\n"
)

# --------------------------------------------------------------
# 18. PAIR-LEVEL SEQUENCE CONCORDANCE
# --------------------------------------------------------------
#
# For each qualifying pair, compare the N- and C-flank
# observations wherever both have the same TAC102 coordinate.
#
# This is useful because Step 4B.3c established that 28 pairs
# have both flanks reaching the gap and therefore span the
# complete 649-732 region at the pair level.

pair_sequence_concordance <- positional_aa %>%
  group_by(
    pair_id,
    tac102_pos
  ) %>%
  summarise(
    n_observations = n(),
    observed_aa = paste(
      sort(unique(amino_acid)),
      collapse = "/"
    ),
    concordant =
      n_distinct(amino_acid) == 1,
    .groups = "drop"
  ) %>%
  group_by(pair_id) %>%
  summarise(

    positions_observed = n(),

    concordant_positions =
      sum(concordant),

    discordant_positions =
      sum(!concordant),

    pair_concordance_pct =
      round(
        concordant_positions /
          positions_observed *
          100,
        1
      ),

    .groups = "drop"
  ) %>%
  arrange(desc(pair_concordance_pct))

# --------------------------------------------------------------
# 19. PRINT PAIR CONCORDANCE
# --------------------------------------------------------------

cat("==============================================================\n")
cat(" PAIR-LEVEL SEQUENCE CONCORDANCE\n")
cat("==============================================================\n\n")

print(
  pair_sequence_concordance,
  n = Inf
)

cat("\n")

# --------------------------------------------------------------
# 20. SAVE POSITIONAL REPORT
# --------------------------------------------------------------

write_csv(
  position_summary_complete,
  output_positional
)

cat(
  "Saved positional report:\n",
  output_positional,
  "\n"
)

# --------------------------------------------------------------
# 21. SAVE PAIR-LEVEL SEQUENCE CONCORDANCE
# --------------------------------------------------------------

write_csv(
  pair_sequence_concordance,
  output_sequence
)

cat(
  "Saved pair-level concordance report:\n",
  output_sequence,
  "\n"
)

# --------------------------------------------------------------
# 22. SAVE CONSENSUS FASTA
# --------------------------------------------------------------

consensus_lines <- c(
  ">TAC102_AA649_732_raw_read_consensus",
  consensus_peptide
)

writeLines(
  consensus_lines,
  output_consensus
)

cat(
  "Saved consensus peptide:\n",
  output_consensus,
  "\n"
)

# --------------------------------------------------------------
# 23. WRITE TEXT REPORT
# --------------------------------------------------------------

report_lines <- c(
  "==============================================================",
  "TAC102 PAIR-TO-PAIR SEQUENCE CONCORDANCE",
  "Step 4B.3d",
  "==============================================================",
  "",
  paste("Target region: TAC102 AA", 649, "-", 732),
  paste("Target length:", total_positions, "aa"),
  "",
  "COHORT",
  paste("Qualifying full-gap pairs:", length(qualifying_pairs)),
  paste("Positionally usable translated sequences:", nrow(usable_tails)),
  "",
  "POSITIONAL COVERAGE",
  paste("Positions covered:", positions_covered),
  paste("Positions missing:", positions_missing),
  paste("Mean depth:", mean_depth),
  paste(
    "Positions >=90% concordant:",
    high_concordance_positions
  ),
  paste(
    "Positions 100% concordant:",
    complete_concordance_positions
  ),
  "",
  "CONSENSUS",
  paste("Consensus length:", nchar(consensus_peptide)),
  consensus_peptide,
  "",
  "METHODOLOGICAL LIMITATION",
  "TAC102 nucleotide reference is unavailable.",
  "Therefore this analysis is not nucleotide-level identity validation.",
  "Sequence reconstruction uses the original-HSP-frame,",
  "codon-anchored translated candidate sequences generated in",
  "Step 4B.3b.",
  "",
  "==============================================================",
  "END Step 4B.3d",
  "=============================================================="
)

writeLines(
  report_lines,
  output_report
)

cat(
  "Saved text report:\n",
  output_report,
  "\n\n"
)

# --------------------------------------------------------------
# 24. FINAL STATUS
# --------------------------------------------------------------

cat("==============================================================\n")
cat("[COMPLETE] Step 4B.3d sequence concordance analysis finished.\n")
cat("==============================================================\n")