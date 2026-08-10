library(dplyr)
library(readr)
library(tidyr)

cat("==============================================================\n")
cat(" Step 4B.3e: TAC102 Sequence Discordance Audit\n")
cat("==============================================================\n\n")

# ==============================================================
# PATHS
# ==============================================================

csv_pairs <- "raw_read_analysis/reports/tac102_41pairs_pair_level_evidence.csv"
csv_ev    <- "raw_read_analysis/reports/tac102_41pairs_gap_evidence.csv"
csv_tails <- "raw_read_analysis/reports/tac102_41pairs_unaligned_tails.csv"

out_positional <- "raw_read_analysis/reports/tac102_4B3e_position_discordance.csv"
out_pairwise   <- "raw_read_analysis/reports/tac102_4B3e_pair_discordance.csv"
out_summary    <- "raw_read_analysis/reports/tac102_4B3e_discordance_summary.txt"

# ==============================================================
# PARAMETERS
# ==============================================================

target_start <- 649
target_end   <- 732
target_length <- target_end - target_start + 1

cat("Target TAC102 coordinate range: AA",
    target_start, "-", target_end,
    "(", target_length, "AA positions)\n\n")

# ==============================================================
# LOAD DATA
# ==============================================================

pairs_df <- read_csv(csv_pairs, show_col_types = FALSE)
ev_df    <- read_csv(csv_ev, show_col_types = FALSE)
tails_df <- read_csv(csv_tails, show_col_types = FALSE)

cat("Input rows:\n")
cat("Pair-level:", nrow(pairs_df), "\n")
cat("Gap-evidence:", nrow(ev_df), "\n")
cat("Unaligned tails:", nrow(tails_df), "\n\n")

# ==============================================================
# REQUIRED COLUMN CHECK
# ==============================================================

required_pairs <- c(
  "pair_id",
  "full_gap_649_732"
)

required_tails <- c(
  "pair_id",
  "orientation_role",
  "gap_position_start",
  "gap_position_end",
  "translated_candidate"
)

missing_pairs <- setdiff(required_pairs, names(pairs_df))
missing_tails <- setdiff(required_tails, names(tails_df))

if (length(missing_pairs) > 0) {
  stop(
    "Missing pair-level columns: ",
    paste(missing_pairs, collapse = ", ")
  )
}

if (length(missing_tails) > 0) {
  stop(
    "Missing tail columns: ",
    paste(missing_tails, collapse = ", ")
  )
}

# ==============================================================
# QUALIFYING PAIRS
# ==============================================================

qualifying_pairs <- pairs_df %>%
  filter(full_gap_649_732 == TRUE) %>%
  pull(pair_id)

cat("Qualifying complete-gap pairs:", length(qualifying_pairs), "\n\n")

# ==============================================================
# SELECT USABLE TRANSLATED SEQUENCES
# ==============================================================

qual_tails <- tails_df %>%
  filter(pair_id %in% qualifying_pairs) %>%
  filter(!is.na(translated_candidate)) %>%
  filter(translated_candidate != "") %>%
  filter(!is.na(gap_position_start)) %>%
  filter(!is.na(gap_position_end)) %>%
  mutate(
    candidate_length = nchar(translated_candidate),
    coordinate_length =
      gap_position_end - gap_position_start + 1
  )

cat("Usable translated sequences:", nrow(qual_tails), "\n")

# ==============================================================
# COORDINATE / SEQUENCE VALIDATION
# ==============================================================

qual_tails <- qual_tails %>%
  mutate(
    length_match = candidate_length == coordinate_length
  )

cat(
  "Sequences with matching sequence/coordinate lengths:",
  sum(qual_tails$length_match),
  "\n"
)

if (any(!qual_tails$length_match)) {

  cat("\nWARNING: Some sequences have inconsistent lengths.\n")

  print(
    qual_tails %>%
      filter(!length_match) %>%
      select(
        pair_id,
        orientation_role,
        gap_position_start,
        gap_position_end,
        candidate_length,
        coordinate_length
      ),
    n = Inf
  )
}

# Retain only mathematically valid sequences
qual_tails <- qual_tails %>%
  filter(length_match)

# ==============================================================
# EXPAND EACH SEQUENCE INTO POSITIONAL OBSERVATIONS
# ==============================================================

positional_aa <- bind_rows(
  lapply(seq_len(nrow(qual_tails)), function(i) {

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

    tibble(
      pair_id = row$pair_id,
      orientation_role = row$orientation_role,
      tac102_pos = pos,
      amino_acid = aa
    )
  })
)

positional_aa <- positional_aa %>%
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
# POSITIONAL DISCORDANCE
# ==============================================================

position_summary <- positional_aa %>%
  group_by(tac102_pos) %>%
  summarise(

    coverage_depth = n(),

    unique_aa_count =
      n_distinct(amino_acid),

    dominant_aa =
      names(which.max(table(amino_acid))),

    dominant_count =
      max(table(amino_acid)),

    dominant_frequency =
      round(
        dominant_count / coverage_depth,
        4
      ),

    concordance_pct =
      round(
        dominant_count / coverage_depth * 100,
        1
      ),

    all_observed_aa =
      paste(
        sort(unique(amino_acid)),
        collapse = "/"
      ),

    .groups = "drop"
  ) %>%
  mutate(

    discordant_count =
      coverage_depth - dominant_count,

    discordance_pct =
      round(
        discordant_count / coverage_depth * 100,
        1
      ),

    concordance_class = case_when(

      concordance_pct == 100 ~ "100%",

      concordance_pct >= 95 ~
        ">=95%",

      concordance_pct >= 90 ~
        "90-94.9%",

      concordance_pct >= 75 ~
        "75-89.9%",

      TRUE ~
        "<75%"
    )
  ) %>%
  arrange(tac102_pos)

# ==============================================================
# PRINT POSITIONAL SUMMARY
# ==============================================================

cat("=== POSITIONAL DISCORDANCE SUMMARY ===\n\n")

print(
  position_summary,
  n = Inf
)

# ==============================================================
# CLASS COUNTS
# ==============================================================

cat("\n=== CONCORDANCE CLASS COUNTS ===\n")

print(
  table(
    position_summary$concordance_class
  )
)

# ==============================================================
# HIGHLY DISCORDANT POSITIONS
# ==============================================================

cat("\n=== MOST DISCORDANT POSITIONS ===\n")

most_discordant <- position_summary %>%
  arrange(
    concordance_pct,
    desc(coverage_depth)
  ) %>%
  select(
    tac102_pos,
    coverage_depth,
    dominant_aa,
    dominant_frequency,
    concordance_pct,
    discordance_pct,
    unique_aa_count,
    all_observed_aa
  )

print(
  most_discordant,
  n = Inf
)

# ==============================================================
# POSITIONS WITH MULTIPLE AMINO ACIDS
# ==============================================================

variable_positions <- position_summary %>%
  filter(unique_aa_count > 1)

cat(
  "\nPositions with >1 observed amino acid:",
  nrow(variable_positions),
  "\n"
)

print(
  variable_positions,
  n = Inf
)

# ==============================================================
# POSITIONS WITH 100% CONCORDANCE
# ==============================================================

perfect_positions <- position_summary %>%
  filter(concordance_pct == 100)

cat(
  "\nPositions with 100% concordance:",
  nrow(perfect_positions),
  "\n"
)

print(
  perfect_positions,
  n = Inf
)

# ==============================================================
# BUILD CONSENSUS
# ==============================================================

consensus_peptide <- paste(
  position_summary$dominant_aa,
  collapse = ""
)

cat("\n=== CONSENSUS PEPTIDE ===\n")
cat(consensus_peptide, "\n\n")

# ==============================================================
# PAIR-LEVEL OBSERVATION MATRIX
# ==============================================================

pair_position_matrix <- positional_aa %>%
  select(
    pair_id,
    tac102_pos,
    amino_acid
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = tac102_pos,
    values_from = amino_acid,
    names_prefix = "AA_"
  )

# ==============================================================
# PAIR-LEVEL CONCORDANCE AGAINST COHORT CONSENSUS
# ==============================================================

consensus_lookup <- position_summary %>%
  select(
    tac102_pos,
    dominant_aa
  )

pair_discordance <- positional_aa %>%
  left_join(
    consensus_lookup,
    by = "tac102_pos"
  ) %>%
  mutate(
    agrees_with_consensus =
      amino_acid == dominant_aa
  ) %>%
  group_by(pair_id) %>%
  summarise(

    positions_observed =
      n(),

    concordant_positions =
      sum(agrees_with_consensus),

    discordant_positions =
      sum(!agrees_with_consensus),

    pair_concordance_pct =
      round(
        concordant_positions /
          positions_observed * 100,
        1
      ),

    .groups = "drop"
  ) %>%
  arrange(pair_concordance_pct)

cat("\n=== PAIR-LEVEL DISCORDANCE ===\n\n")

print(
  pair_discordance,
  n = Inf
)

# ==============================================================
# DISCORDANCE BY N/C FLANK
# ==============================================================

flank_discordance <- positional_aa %>%
  left_join(
    consensus_lookup,
    by = "tac102_pos"
  ) %>%
  mutate(
    agrees_with_consensus =
      amino_acid == dominant_aa
  ) %>%
  group_by(orientation_role) %>%
  summarise(

    observations = n(),

    concordant =
      sum(agrees_with_consensus),

    discordant =
      sum(!agrees_with_consensus),

    concordance_pct =
      round(
        concordant / observations * 100,
        1
      ),

    .groups = "drop"
  )

cat("\n=== N/C FLANK CONCORDANCE ===\n\n")

print(flank_discordance)

# ==============================================================
# GLOBAL METRICS
# ==============================================================

total_positions <- nrow(position_summary)

mean_depth <- round(
  mean(position_summary$coverage_depth),
  2
)

positions_100 <- sum(
  position_summary$concordance_pct == 100
)

positions_90 <- sum(
  position_summary$concordance_pct >= 90
)

positions_variable <- sum(
  position_summary$unique_aa_count > 1
)

mean_concordance <- round(
  mean(position_summary$concordance_pct),
  2
)

cat("\n=== GLOBAL DISCORDANCE METRICS ===\n\n")

cat(
  "Target positions:",
  target_length,
  "\n"
)

cat(
  "Positions observed:",
  total_positions,
  "\n"
)

cat(
  "Mean depth:",
  mean_depth,
  "\n"
)

cat(
  "Mean positional concordance:",
  mean_concordance,
  "%\n"
)

cat(
  "Positions >=90% concordance:",
  positions_90,
  "\n"
)

cat(
  "Positions with 100% concordance:",
  positions_100,
  "\n"
)

cat(
  "Variable positions:",
  positions_variable,
  "\n"
)

# ==============================================================
# SAVE REPORTS
# ==============================================================

write_csv(
  position_summary,
  out_positional
)

write_csv(
  pair_discordance,
  out_pairwise
)

writeLines(
  c(
    "Step 4B.3e: TAC102 Sequence Discordance Audit",
    "==============================================================",
    "",
    paste(
      "Target range:",
      target_start,
      "-",
      target_end
    ),
    paste(
      "Target length:",
      target_length,
      "AA"
    ),
    paste(
      "Qualifying pairs:",
      length(qualifying_pairs)
    ),
    paste(
      "Usable translated sequences:",
      nrow(qual_tails)
    ),
    paste(
      "Expanded positional observations:",
      nrow(positional_aa)
    ),
    "",
    paste(
      "Positions observed:",
      total_positions
    ),
    paste(
      "Mean depth:",
      mean_depth
    ),
    paste(
      "Mean positional concordance:",
      mean_concordance,
      "%"
    ),
    paste(
      "Positions >=90% concordance:",
      positions_90
    ),
    paste(
      "Positions with 100% concordance:",
      positions_100
    ),
    paste(
      "Variable positions:",
      positions_variable
    ),
    "",
    "Consensus peptide:",
    consensus_peptide
  ),
  out_summary
)

cat("\nSaved:\n")
cat(out_positional, "\n")
cat(out_pairwise, "\n")
cat(out_summary, "\n")

cat("\n[COMPLETE] Step 4B.3e discordance audit finished.\n")