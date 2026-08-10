library(dplyr)
library(readr)
library(tidyr)
library(stringr)

cat("==============================================================\n")
cat(" Step 4B.3g: TAC102 Internal Pair-Conflict Resolution Audit\n")
cat("==============================================================\n\n")

# ==============================================================
# 0. PATHS
# ==============================================================

csv_conflicts <- "raw_read_analysis/reports/tac102_4B3f_pair_conflicts.csv"
csv_tails     <- "raw_read_analysis/reports/tac102_41pairs_unaligned_tails.csv"

out_audit   <- "raw_read_analysis/reports/tac102_4B3g_conflict_detailed_audit.csv"
out_class   <- "raw_read_analysis/reports/tac102_4B3g_conflict_classification_summary.csv"
out_summary <- "raw_read_analysis/reports/tac102_4B3g_resolution_summary.txt"

target_start <- 649
target_end   <- 732

# ==============================================================
# 1. INPUT VALIDATION
# ==============================================================

if (!file.exists(csv_conflicts)) {
  stop(
    "Missing conflict file:\n",
    csv_conflicts,
    "\nRun Step 4B.3f first."
  )
}

if (!file.exists(csv_tails)) {
  stop(
    "Missing tail file:\n",
    csv_tails,
    "\nRun Step 4B.3b first."
  )
}

conflicts_df <- readr::read_csv(
  csv_conflicts,
  show_col_types = FALSE
)

tails_df <- readr::read_csv(
  csv_tails,
  show_col_types = FALSE
)

cat("Input rows:\n")
cat("  Internal conflicts:", nrow(conflicts_df), "\n")
cat("  Tail observations:", nrow(tails_df), "\n\n")

# ==============================================================
# 2. REQUIRED COLUMN VALIDATION
# ==============================================================

required_conflict <- c(
  "pair_id",
  "tac102_pos",
  "aa_N",
  "aa_C"
)

missing_conflict <- setdiff(
  required_conflict,
  names(conflicts_df)
)

if (length(missing_conflict) > 0) {
  stop(
    "Missing conflict columns: ",
    paste(missing_conflict, collapse = ", ")
  )
}

required_tail <- c(
  "pair_id",
  "orientation_role",
  "gap_position_start",
  "gap_position_end",
  "translated_candidate"
)

missing_tail <- setdiff(
  required_tail,
  names(tails_df)
)

if (length(missing_tail) > 0) {
  stop(
    "Missing tail columns: ",
    paste(missing_tail, collapse = ", ")
  )
}

# ==============================================================
# 3. IDENTIFY NUCLEOTIDE SEQUENCE COLUMN
# ==============================================================

nt_candidates <- c(
  "gap_facing_sequence",
  "unaligned_seq",
  "gap_facing_nt_sequence",
  "sequence"
)

available_nt <- nt_candidates[
  nt_candidates %in% names(tails_df)
]

if (length(available_nt) == 0) {
  stop(
    "No nucleotide sequence column found. Expected one of:\n",
    paste(nt_candidates, collapse = ", ")
  )
}

nt_column <- available_nt[1]

cat("Nucleotide sequence column:", nt_column, "\n\n")

# ==============================================================
# 4. STANDARDIZE ORIENTATION
# ==============================================================

tails_df <- tails_df %>%
  dplyr::mutate(
    orientation_role = dplyr::case_when(
      orientation_role %in% c("N", "N_flank") ~ "N",
      orientation_role %in% c("C", "C_flank") ~ "C",
      TRUE ~ as.character(orientation_role)
    )
  )

# ==============================================================
# 5. PREPARE USABLE TAILS
# ==============================================================

tails_use <- tails_df %>%
  dplyr::filter(
    orientation_role %in% c("N", "C"),
    !is.na(translated_candidate),
    translated_candidate != "",
    !is.na(gap_position_start),
    !is.na(gap_position_end)
  ) %>%
  dplyr::mutate(
    aa_length = nchar(translated_candidate),

    coordinate_length =
      gap_position_end - gap_position_start + 1,

    nt_sequence = .data[[nt_column]],

    nt_length = dplyr::if_else(
      is.na(nt_sequence),
      NA_integer_,
      nchar(nt_sequence)
    )
  ) %>%
  dplyr::filter(
    aa_length == coordinate_length
  )

cat("Usable tail observations:", nrow(tails_use), "\n")
cat(
  "N-flank observations:",
  sum(tails_use$orientation_role == "N"),
  "\n"
)
cat(
  "C-flank observations:",
  sum(tails_use$orientation_role == "C"),
  "\n\n"
)

# ==============================================================
# 6. BUILD N AND C TABLES
# ==============================================================

n_tail <- tails_use %>%
  dplyr::filter(orientation_role == "N") %>%
  dplyr::select(
    pair_id,
    start_N = gap_position_start,
    end_N   = gap_position_end,
    seq_N_aa = translated_candidate,
    seq_N_nt = nt_sequence
  )

c_tail <- tails_use %>%
  dplyr::filter(orientation_role == "C") %>%
  dplyr::select(
    pair_id,
    start_C = gap_position_start,
    end_C   = gap_position_end,
    seq_C_aa = translated_candidate,
    seq_C_nt = nt_sequence
  )

# ==============================================================
# 7. JOIN CONFLICTS WITH N/C EVIDENCE
# ==============================================================

conflict_joined <- conflicts_df %>%
  dplyr::left_join(n_tail, by = "pair_id") %>%
  dplyr::left_join(c_tail, by = "pair_id")

both_present <- sum(
  !is.na(conflict_joined$seq_N_aa) &
  !is.na(conflict_joined$seq_C_aa)
)

cat(
  "Conflicts with both N and C sequence records:",
  both_present,
  "\n\n"
)

# ==============================================================
# 8. CODON EXTRACTION
# ==============================================================

extract_codon <- function(
    nt_sequence,
    aa_start,
    aa_end,
    tac102_position
) {

  if (
    is.na(nt_sequence) ||
    is.na(aa_start) ||
    is.na(aa_end) ||
    is.na(tac102_position)
  ) {
    return(NA_character_)
  }

  relative_aa <-
    tac102_position - aa_start + 1

  if (
    relative_aa < 1 ||
    relative_aa > (aa_end - aa_start + 1)
  ) {
    return(NA_character_)
  }

  nt_start <-
    (relative_aa - 1) * 3 + 1

  nt_end <-
    nt_start + 2

  if (nt_end > nchar(nt_sequence)) {
    return(NA_character_)
  }

  substr(
    nt_sequence,
    nt_start,
    nt_end
  )
}

# ==============================================================
# 9. FRAME DIAGNOSTIC
# ==============================================================

frame_status <- function(nt_sequence) {

  if (is.na(nt_sequence)) {
    return("missing_nt")
  }

  if (nchar(nt_sequence) %% 3 == 0) {
    return("frame_consistent")
  }

  "frame_inconsistent"
}

# ==============================================================
# 10. CODON CLASSIFICATION HELPERS
# ==============================================================

is_lys_codon <- function(codon) {

  if (is.na(codon)) {
    return(FALSE)
  }

  codon %in% c("AAA", "AAG")
}

is_glu_codon <- function(codon) {

  if (is.na(codon)) {
    return(FALSE)
  }

  codon %in% c("GAA", "GAG")
}

is_low_complexity_codon <- function(codon) {

  if (is.na(codon)) {
    return(FALSE)
  }

  codon %in% c(
    "AAA",
    "AAG",
    "GAA",
    "GAG"
  )
}

# ==============================================================
# 11. CONFLICT RESOLUTION FUNCTION
# ==============================================================

resolve_one_conflict <- function(row) {

  pos <- row$tac102_pos

  codon_N <- extract_codon(
    row$seq_N_nt,
    row$start_N,
    row$end_N,
    pos
  )

  codon_C <- extract_codon(
    row$seq_C_nt,
    row$start_C,
    row$end_C,
    pos
  )

  frame_N <- frame_status(row$seq_N_nt)
  frame_C <- frame_status(row$seq_C_nt)

  # ------------------------------------------------------------
  # Determine overlap
  # ------------------------------------------------------------

  if (
    is.na(row$start_N) ||
    is.na(row$start_C) ||
    is.na(row$end_N) ||
    is.na(row$end_C)
  ) {

    overlap_start <- NA_integer_
    overlap_end   <- NA_integer_
    overlap_len   <- 0L

  } else {

    overlap_start <-
      max(row$start_N, row$start_C)

    overlap_end <-
      min(row$end_N, row$end_C)

    overlap_len <-
      max(
        0L,
        overlap_end - overlap_start + 1
      )
  }

  # ------------------------------------------------------------
  # Classification
  # ------------------------------------------------------------

  classification <- dplyr::case_when(

    is.na(row$aa_N) ||
      is.na(row$aa_C) ~
      "missing_aa_call",

    is.na(codon_N) ||
      is.na(codon_C) ~
      "unresolved_missing_nt",

    codon_N == codon_C &&
      row$aa_N != row$aa_C ~
      "nucleotide_identical_translation_discrepancy",

    frame_N == "frame_inconsistent" ||
      frame_C == "frame_inconsistent" ~
      "frame_inconsistent_tail",

    codon_N != codon_C &&
      row$aa_N == row$aa_C ~
      "synonymous_nucleotide_difference",

    is_low_complexity_codon(codon_N) &&
      is_low_complexity_codon(codon_C) ~
      "low_complexity_Lys_Glu_region",

    codon_N != codon_C &&
      row$aa_N != row$aa_C ~
      "nonsynonymous_nucleotide_difference",

    TRUE ~
      "unresolved"
  )

  tibble::tibble(
    pair_id = row$pair_id,
    tac102_pos = pos,

    aa_N = row$aa_N,
    aa_C = row$aa_C,

    codon_N = codon_N,
    codon_C = codon_C,

    frame_N = frame_N,
    frame_C = frame_C,

    nt_length_N =
      ifelse(
        is.na(row$seq_N_nt),
        NA_integer_,
        nchar(row$seq_N_nt)
      ),

    nt_length_C =
      ifelse(
        is.na(row$seq_C_nt),
        NA_integer_,
        nchar(row$seq_C_nt)
      ),

    overlap_start = overlap_start,
    overlap_end   = overlap_end,
    overlap_len_aa = overlap_len,

    N_is_Lys =
      is_lys_codon(codon_N),

    C_is_Lys =
      is_lys_codon(codon_C),

    N_is_Glu =
      is_glu_codon(codon_N),

    C_is_Glu =
      is_glu_codon(codon_C),

    conflict_class = classification
  )
}

# ==============================================================
# 12. RUN THE AUDIT
# ==============================================================

detailed_audit <- dplyr::bind_rows(
  lapply(
    seq_len(nrow(conflict_joined)),
    function(i) {

      resolve_one_conflict(
        conflict_joined[i, ]
      )

    }
  )
)

cat(
  "Detailed conflict records generated:",
  nrow(detailed_audit),
  "\n\n"
)

# ==============================================================
# 13. CLASSIFICATION SUMMARY
# ==============================================================

class_summary <- detailed_audit %>%
  dplyr::group_by(conflict_class) %>%
  dplyr::summarise(
    conflict_count = dplyr::n(),
    percentage =
      round(
        dplyr::n() /
          nrow(detailed_audit) *
          100,
        2
      ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(conflict_count)
  )

cat("=== CONFLICT CLASSIFICATION SUMMARY ===\n\n")

print(class_summary)

cat("\n")

# ==============================================================
# 14. FRAME / LOW-COMPLEXITY SUMMARY
# ==============================================================

diagnostic_summary <- detailed_audit %>%
  dplyr::summarise(

    total_conflicts = dplyr::n(),

    frame_inconsistent_N =
      sum(
        frame_N == "frame_inconsistent",
        na.rm = TRUE
      ),

    frame_inconsistent_C =
      sum(
        frame_C == "frame_inconsistent",
        na.rm = TRUE
      ),

    Lys_N =
      sum(
        N_is_Lys,
        na.rm = TRUE
      ),

    Lys_C =
      sum(
        C_is_Lys,
        na.rm = TRUE
      ),

    Glu_N =
      sum(
        N_is_Glu,
        na.rm = TRUE
      ),

    Glu_C =
      sum(
        C_is_Glu,
        na.rm = TRUE
      ),

    low_complexity_conflicts =
      sum(
        conflict_class ==
          "low_complexity_Lys_Glu_region"
      )
  )

cat("=== FRAME / LOW-COMPLEXITY SUMMARY ===\n\n")

print(diagnostic_summary)

cat("\n")

# ==============================================================
# 15. MOST COMMON AA CONFLICTS
# ==============================================================

aa_conflict_summary <- detailed_audit %>%
  dplyr::group_by(
    aa_N,
    aa_C
  ) %>%
  dplyr::summarise(
    conflict_count = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(conflict_count)
  )

cat("=== MOST COMMON AA CONFLICTS ===\n\n")

print(
  head(
    aa_conflict_summary,
    20
  )
)

cat("\n")

# ==============================================================
# 16. POSITIONAL SUMMARY
# ==============================================================

position_summary <- detailed_audit %>%
  dplyr::group_by(tac102_pos) %>%
  dplyr::summarise(

    conflicts = dplyr::n(),

    unique_N =
      dplyr::n_distinct(aa_N),

    unique_C =
      dplyr::n_distinct(aa_C),

    N_calls =
      paste(
        sort(unique(aa_N)),
        collapse = "/"
      ),

    C_calls =
      paste(
        sort(unique(aa_C)),
        collapse = "/"
      ),

    classes =
      paste(
        sort(unique(conflict_class)),
        collapse = ";"
      ),

    .groups = "drop"
  ) %>%
  dplyr::arrange(tac102_pos)

# ==============================================================
# 17. SAVE OUTPUTS
# ==============================================================

readr::write_csv(
  detailed_audit,
  out_audit
)

readr::write_csv(
  class_summary,
  out_class
)

# ==============================================================
# 18. TEXT SUMMARY
# ==============================================================

summary_text <- c(

  "Step 4B.3g: TAC102 Internal Pair-Conflict Resolution Audit",

  "==============================================================",

  paste(
    "Total internal N/C conflicts:",
    nrow(detailed_audit)
  ),

  paste(
    "Nucleotide sequence column:",
    nt_column
  ),

  "",

  "Conflict classification summary:",

  capture.output(
    print(class_summary)
  ),

  "",

  "Frame / low-complexity diagnostics:",

  capture.output(
    print(diagnostic_summary)
  ),

  "",

  "Most common amino-acid conflicts:",

  capture.output(
    print(
      head(
        aa_conflict_summary,
        20
      )
    )
  ),

  "",

  "Outputs:",

  out_audit,
  out_class,
  out_summary
)

writeLines(
  summary_text,
  out_summary
)

# ==============================================================
# 19. FINAL MESSAGE
# ==============================================================

cat("Saved:\n")
cat(" ", out_audit, "\n")
cat(" ", out_class, "\n")
cat(" ", out_summary, "\n\n")

cat(
  "[COMPLETE] Step 4B.3g conflict resolution audit finished.\n"
)