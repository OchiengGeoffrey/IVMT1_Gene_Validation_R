library(dplyr)
library(readr)
library(tidyr)

cat("==============================================================\n")
cat(" Step 4B.3f: Pair-Normalized TAC102 Sequence Concordance\n")
cat("==============================================================\n\n")


# ==============================================================
# 0. PATHS AND TARGET
# ==============================================================

csv_pairs <- "raw_read_analysis/reports/tac102_41pairs_pair_level_evidence.csv"
csv_tails <- "raw_read_analysis/reports/tac102_41pairs_unaligned_tails.csv"

out_positional <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3f_pair_normalized_position_concordance.csv"
)

out_conflicts <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3f_pair_conflicts.csv"
)

out_summary <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3f_pair_normalized_summary.txt"
)

target_start  <- 649
target_end    <- 732
target_length <- target_end - target_start + 1


cat(
  "Target TAC102 coordinate range: AA",
  target_start, "-", target_end,
  "(", target_length, "AA positions)\n\n"
)


# ==============================================================
# 1. LOAD DATA
# ==============================================================

if (!file.exists(csv_pairs)) {
  stop("Missing file: ", csv_pairs)
}

if (!file.exists(csv_tails)) {
  stop("Missing file: ", csv_tails)
}

pairs_df <- read_csv(
  csv_pairs,
  show_col_types = FALSE
)

tails_df <- read_csv(
  csv_tails,
  show_col_types = FALSE
)

cat("Input rows:\n")
cat("Pair-level:", nrow(pairs_df), "\n")
cat("Unaligned tails:", nrow(tails_df), "\n\n")


# ==============================================================
# 2. VALIDATE REQUIRED COLUMNS
# ==============================================================

required_pair_cols <- c(
  "pair_id",
  "full_gap_649_732"
)

required_tail_cols <- c(
  "pair_id",
  "orientation_role",
  "translated_candidate",
  "gap_position_start",
  "gap_position_end"
)

missing_pair_cols <- setdiff(
  required_pair_cols,
  names(pairs_df)
)

missing_tail_cols <- setdiff(
  required_tail_cols,
  names(tails_df)
)

if (length(missing_pair_cols) > 0) {
  stop(
    "Missing pair-level columns: ",
    paste(missing_pair_cols, collapse = ", ")
  )
}

if (length(missing_tail_cols) > 0) {
  stop(
    "Missing tail columns: ",
    paste(missing_tail_cols, collapse = ", ")
  )


}


# ==============================================================
# 3. VALIDATE ORIENTATION CODING
# ==============================================================

cat("=== ORIENTATION ROLE VALIDATION ===\n")

print(
  table(
    tails_df$orientation_role,
    useNA = "ifany"
  )
)

valid_roles <- c("N", "C")

observed_roles <- unique(
  na.omit(tails_df$orientation_role)
)

unexpected_roles <- setdiff(
  observed_roles,
  valid_roles
)

if (length(unexpected_roles) > 0) {

  stop(
    "Unexpected orientation_role values detected: ",
    paste(unexpected_roles, collapse = ", "),
    "\nExpected values are exactly: N and C."
  )

}

cat(
  "Orientation coding validated: N and C.\n\n"
)


# ==============================================================
# 4. IDENTIFY QUALIFYING COMPLETE-GAP PAIRS
# ==============================================================

qualifying_pairs <- pairs_df %>%
  filter(full_gap_649_732 == TRUE) %>%
  pull(pair_id) %>%
  unique()

cat(
  "Qualifying complete-gap pairs:",
  length(qualifying_pairs),
  "\n"
)

cat("Expected: 28\n\n")

if (length(qualifying_pairs) != 28) {
  warning(
    "Expected 28 qualifying pairs but found ",
    length(qualifying_pairs)
  )
}


# ==============================================================
# 5. RESTRICT TAILS TO QUALIFYING PAIRS
# ==============================================================

qual_tails <- tails_df %>%
  filter(pair_id %in% qualifying_pairs) %>%
  filter(
    !is.na(translated_candidate),
    translated_candidate != "",
    !is.na(gap_position_start),
    !is.na(gap_position_end)
  ) %>%
  mutate(
    candidate_length = nchar(translated_candidate),
    coordinate_length =
      gap_position_end - gap_position_start + 1
  )

cat(
  "Tail rows belonging to qualifying pairs:",
  nrow(qual_tails),
  "\n"
)

cat(
  "Expected approximately:",
  length(qualifying_pairs) * 2,
  "\n\n"
)


# ==============================================================
# 6. SEQUENCE / COORDINATE VALIDATION
# ==============================================================

invalid_rows <- qual_tails %>%
  filter(candidate_length != coordinate_length)

cat("=== SEQUENCE / COORDINATE VALIDATION ===\n")

cat(
  "Rows with valid coordinates:",
  nrow(qual_tails),
  "\n"
)

cat(
  "Rows with matching sequence/coordinate lengths:",
  sum(
    qual_tails$candidate_length ==
      qual_tails$coordinate_length
  ),
  "\n"
)

if (nrow(invalid_rows) > 0) {

  print(
    invalid_rows %>%
      select(
        pair_id,
        orientation_role,
        gap_position_start,
        gap_position_end,
        candidate_length,
        coordinate_length
      )
  )

}

qual_tails <- qual_tails %>%
  filter(candidate_length == coordinate_length)

cat(
  "Positionally usable translated sequences:",
  nrow(qual_tails),
  "\n\n"
)


# ==============================================================
# 7. EXPAND SEQUENCES TO TAC102 POSITIONS
# ==============================================================

positional_list <- lapply(
  seq_len(nrow(qual_tails)),
  function(i) {

    row <- qual_tails[i, ]

    aa <- strsplit(
      row$translated_candidate,
      "",
      fixed = TRUE
    )[[1]]

    pos <- seq(
      row$gap_position_start,
      row$gap_position_end
    )

    if (length(aa) != length(pos)) {
      return(NULL)
    }

    tibble(
      pair_id = row$pair_id,
      orientation_role = row$orientation_role,
      tac102_pos = pos,
      amino_acid = aa
    )
  }
)

positional_aa <- bind_rows(positional_list) %>%
  filter(
    tac102_pos >= target_start,
    tac102_pos <= target_end
  )

cat(
  "Expanded positional observations:",
  nrow(positional_aa),
  "\n\n"
)


# ==============================================================
# 8. COLLAPSE TO ONE CALL PER PAIR / POSITION
# ==============================================================

pair_pos_resolved <- positional_aa %>%

  group_by(
    pair_id,
    tac102_pos
  ) %>%

  summarise(

    has_N =
      any(orientation_role == "N"),

    has_C =
      any(orientation_role == "C"),

    aa_N = {
      x <- amino_acid[
        orientation_role == "N"
      ]

      if (length(x) == 0) {
        NA_character_
      } else {
        x[1]
      }
    },

    aa_C = {
      x <- amino_acid[
        orientation_role == "C"
      ]

      if (length(x) == 0) {
        NA_character_
      } else {
        x[1]
      }
    },

    .groups = "drop"
  ) %>%

  mutate(

    status = case_when(

      has_N & !has_C ~
        "N_only",

      !has_N & has_C ~
        "C_only",

      has_N & has_C &
        !is.na(aa_N) &
        !is.na(aa_C) &
        aa_N == aa_C ~
        "NC_concordant",

      has_N & has_C &
        !is.na(aa_N) &
        !is.na(aa_C) &
        aa_N != aa_C ~
        "NC_conflict",

      TRUE ~
        "Unknown"
    )
  )


cat("=== PAIR/POSITION STATUS ===\n")

print(
  table(
    pair_pos_resolved$status,
    useNA = "ifany"
  )
)

cat("\n")


# ==============================================================
# 9. HARD VALIDATION: NO UNKNOWN CALLS
# ==============================================================

unknown_count <- sum(
  pair_pos_resolved$status == "Unknown"
)

if (unknown_count > 0) {

  stop(
    "ANALYTICAL ERROR: ",
    unknown_count,
    " pair-position observations were classified as Unknown."
  )

}


# ==============================================================
# 10. INTERNAL N/C CONFLICTS
# ==============================================================

pair_conflicts <- pair_pos_resolved %>%
  filter(status == "NC_conflict") %>%
  select(
    pair_id,
    tac102_pos,
    aa_N,
    aa_C
  ) %>%
  arrange(
    pair_id,
    tac102_pos
  )

cat(
  "Internal N/C overlap conflicts:",
  nrow(pair_conflicts),
  "\n\n"
)


# ==============================================================
# 11. GENERATE ONE PAIR-LEVEL CALL PER POSITION
# ==============================================================

pair_calls_for_consensus <- bind_rows(

  pair_pos_resolved %>%

    filter(
      status %in%
        c(
          "N_only",
          "C_only",
          "NC_concordant"
        )
    ) %>%

    mutate(

      aa = case_when(
        status == "N_only" ~ aa_N,
        status == "C_only" ~ aa_C,
        status == "NC_concordant" ~ aa_N,
        TRUE ~ NA_character_
      ),

      weight = 1.0

    ) %>%

    select(
      pair_id,
      tac102_pos,
      aa,
      weight
    ),

  pair_pos_resolved %>%

    filter(status == "NC_conflict") %>%

    transmute(
      pair_id,
      tac102_pos,
      aa = aa_N,
      weight = 0.5
    ),

  pair_pos_resolved %>%

    filter(status == "NC_conflict") %>%

    transmute(
      pair_id,
      tac102_pos,
      aa = aa_C,
      weight = 0.5
    )
) %>%

  filter(
    !is.na(aa),
    aa != ""
  )


# ==============================================================
# 12. POSITIONAL PAIR DEPTH
# ==============================================================

position_summary <- pair_pos_resolved %>%

  group_by(tac102_pos) %>%

  summarise(

    pair_depth =
      n_distinct(pair_id),

    n_only_count =
      sum(status == "N_only"),

    c_only_count =
      sum(status == "C_only"),

    nc_concordant_count =
      sum(status == "NC_concordant"),

    nc_conflict_count =
      sum(status == "NC_conflict"),

    .groups = "drop"
  )


# ==============================================================
# 13. WEIGHTED CONSENSUS
# ==============================================================

position_consensus <- pair_calls_for_consensus %>%

  group_by(
    tac102_pos,
    aa
  ) %>%

  summarise(
    weighted_count = sum(weight),
    .groups = "drop"
  ) %>%

  group_by(tac102_pos) %>%

  summarise(

    unique_aa_count =
      n_distinct(aa),

    dominant_aa =
      aa[
        which.max(weighted_count)
      ],

    dominant_pair_count =
      max(weighted_count),

    total_weight =
      sum(weighted_count),

    all_observed_aa =
      paste(
        sort(unique(aa)),
        collapse = "/"
      ),

    .groups = "drop"
  )


# ==============================================================
# 14. FINAL POSITIONAL TABLE
# ==============================================================

pair_norm_summary <- position_summary %>%

  left_join(
    position_consensus,
    by = "tac102_pos"
  ) %>%

  mutate(

    pair_concordance_pct =
      round(
        dominant_pair_count /
          pair_depth *
          100,
        1
      ),

    discordant_pair_count =
      pair_depth -
      dominant_pair_count,

    discordance_pct =
      round(
        discordant_pair_count /
          pair_depth *
          100,
        1
      )
  ) %>%

  arrange(tac102_pos)


# ==============================================================
# 15. FINAL VALIDATION
# ==============================================================

if (
  any(
    is.na(pair_norm_summary$dominant_aa)
  )
) {

  stop(
    "ANALYTICAL ERROR: one or more target positions lack ",
    "a dominant amino-acid call."
  )

}

max_pair_depth <- max(
  pair_norm_summary$pair_depth,
  na.rm = TRUE
)

if (
  max_pair_depth >
    length(qualifying_pairs)
) {

  stop(
    "ANALYTICAL ERROR: pair depth exceeds number of ",
    "independent qualifying pairs."
  )

}

if (
  nrow(pair_norm_summary) != target_length
) {

  stop(
    "Expected ",
    target_length,
    " positions but found ",
    nrow(pair_norm_summary)
  )

}


# ==============================================================
# 16. CONSENSUS PEPTIDE
# ==============================================================

consensus_peptide <- paste(
  pair_norm_summary$dominant_aa,
  collapse = ""
)


cat("=== PAIR-NORMALIZED CONSENSUS PEPTIDE (649-732) ===\n")
cat(consensus_peptide, "\n\n")


# ==============================================================
# 17. GLOBAL METRICS
# ==============================================================

mean_pair_depth <- round(
  mean(pair_norm_summary$pair_depth),
  2
)

mean_pair_concordance <- round(
  mean(pair_norm_summary$pair_concordance_pct),
  2
)

positions_ge90 <- sum(
  pair_norm_summary$pair_concordance_pct >= 90
)

positions_ge95 <- sum(
  pair_norm_summary$pair_concordance_pct >= 95
)

positions_100 <- sum(
  pair_norm_summary$pair_concordance_pct == 100
)

positions_lt75 <- sum(
  pair_norm_summary$pair_concordance_pct < 75
)


cat("=== PAIR-NORMALIZED GLOBAL METRICS ===\n")

cat(
  "Target positions:",
  target_length,
  "\n"
)

cat(
  "Mean pair depth per position (max 28):",
  mean_pair_depth,
  "\n"
)

cat(
  "Maximum pair depth:",
  max_pair_depth,
  "\n"
)

cat(
  "Mean pair-normalized concordance:",
  mean_pair_concordance,
  "%\n"
)

cat(
  "Positions >=90% concordance:",
  positions_ge90,
  "\n"
)

cat(
  "Positions >=95% concordance:",
  positions_ge95,
  "\n"
)

cat(
  "Positions with 100% concordance:",
  positions_100,
  "\n"
)

cat(
  "Positions <75% concordance:",
  positions_lt75,
  "\n"
)

cat(
  "Internal N/C pair conflicts:",
  nrow(pair_conflicts),
  "\n\n"
)


# ==============================================================
# 18. PAIR-LEVEL SUMMARY
# ==============================================================

pair_level_summary <- pair_pos_resolved %>%

  group_by(pair_id) %>%

  summarise(

    positions_observed =
      n_distinct(tac102_pos),

    N_only_positions =
      sum(status == "N_only"),

    C_only_positions =
      sum(status == "C_only"),

    NC_concordant_positions =
      sum(status == "NC_concordant"),

    NC_conflict_positions =
      sum(status == "NC_conflict"),

    .groups = "drop"
  ) %>%

  mutate(

    usable_positions =
      N_only_positions +
      C_only_positions +
      NC_concordant_positions +
      NC_conflict_positions,

    pair_internal_concordance_pct =
      round(
        (
          N_only_positions +
          C_only_positions +
          NC_concordant_positions
        ) /
        usable_positions *
        100,
        1
      )
  )


# ==============================================================
# 19. SAVE OUTPUTS
# ==============================================================

write_csv(
  pair_norm_summary,
  out_positional
)

write_csv(
  pair_conflicts,
  out_conflicts
)


summary_text <- c(

  "Step 4B.3f: Pair-Normalized TAC102 Sequence Concordance",

  "==============================================================",

  "",

  paste(
    "Target range:",
    target_start,
    "-",
    target_end,
    "(",
    target_length,
    "AA)"
  ),

  paste(
    "Qualifying complete-gap pairs:",
    length(qualifying_pairs)
  ),

  paste(
    "Maximum theoretical pair depth:",
    length(qualifying_pairs)
  ),

  paste(
    "Observed maximum pair depth:",
    max_pair_depth
  ),

  paste(
    "Mean pair depth per position:",
    mean_pair_depth
  ),

  paste(
    "Mean pair-normalized concordance:",
    mean_pair_concordance,
    "%"
  ),

  paste(
    "Positions >=90% concordance:",
    positions_ge90
  ),

  paste(
    "Positions >=95% concordance:",
    positions_ge95
  ),

  paste(
    "Positions with 100% concordance:",
    positions_100
  ),

  paste(
    "Positions <75% concordance:",
    positions_lt75
  ),

  paste(
    "Internal N/C pair conflicts:",
    nrow(pair_conflicts)
  ),

  "",

  "Pair-Normalized Consensus Peptide (649-732):",

  consensus_peptide,

  "",

  "Analytical unit:",

  "Complete paired-end read pair.",

  "",

  "N and C flank observations from the same pair are not treated",

  "as independent observations. When both cover the same position",

  "they are collapsed to one pair-level call.",

  "Exact N/C disagreements are retained as internal pair conflicts.",

  "",

  "Outputs:",

  out_positional,

  out_conflicts,

  out_summary
)

writeLines(
  summary_text,
  out_summary
)


cat("Saved reports:\n")
cat(" ", out_positional, "\n")
cat(" ", out_conflicts, "\n")
cat(" ", out_summary, "\n")

cat(
  "\n[COMPLETE] Step 4B.3f pair-normalized concordance analysis finished.\n"
)