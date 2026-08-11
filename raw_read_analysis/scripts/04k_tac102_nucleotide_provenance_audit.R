# ==============================================================
# Step 4B.3k: TAC102 Targeted Nucleotide Provenance Audit
# ==============================================================
#
# Purpose:
#   Determine the provenance of the remaining same-length
#   amino-acid discordance cases identified in Step 4B.3j.
#
# This audit specifically asks:
#
#   1. Does the nucleotide codon extracted from the stored
#      gap-facing sequence actually encode the reported AA?
#
#   2. Can a +1 or +2 nucleotide local frame interpretation
#      recover the reported AA?
#
#   3. Does the reconstructed codon differ from a codon capable
#      of encoding the reported AA by only one nucleotide?
#
#   4. Does the evidence instead indicate a larger extraction/
#      provenance discrepancy?
#
# IMPORTANT:
#   An amino-acid mismatch alone is NOT treated as evidence of
#   a nucleotide mutation.
#
# ==============================================================


# ==============================================================
# 1. LOAD REQUIRED PACKAGES
# ==============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})


# ==============================================================
# 2. PATHS
# ==============================================================

csv_audit_detail <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3j_discordance_detailed_audit.csv"
)

csv_tails <- paste0(
  "raw_read_analysis/reports/",
  "tac102_41pairs_unaligned_tails.csv"
)

out_detail <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3k_provenance_detailed_audit.csv"
)

out_class <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3k_provenance_classification_summary.csv"
)

out_sum <- paste0(
  "raw_read_analysis/reports/",
  "tac102_4B3k_provenance_cause_summary.txt"
)


# ==============================================================
# 3. HEADER
# ==============================================================

cat("==============================================================\n")
cat(" Step 4B.3k: TAC102 Targeted Nucleotide Provenance Audit\n")
cat("==============================================================\n\n")


# ==============================================================
# 4. CHECK INPUT FILES
# ==============================================================

if (!file.exists(csv_audit_detail)) {
  stop(
    "Missing required 4B.3j audit file:\n",
    csv_audit_detail
  )
}

if (!file.exists(csv_tails)) {
  stop(
    "Missing required raw-tail file:\n",
    csv_tails
  )
}


# ==============================================================
# 5. LOAD INPUT DATA
# ==============================================================

audit_3j <- readr::read_csv(
  csv_audit_detail,
  show_col_types = FALSE
)

raw_tails <- readr::read_csv(
  csv_tails,
  show_col_types = FALSE
)


cat("4B.3j audit rows:", nrow(audit_3j), "\n")
cat("Raw tail rows:", nrow(raw_tails), "\n\n")


# ==============================================================
# 6. CHECK REQUIRED COLUMNS
# ==============================================================

required_audit_columns <- c(
  "pair_id",
  "mate",
  "orientation_role",
  "gap_position_start",
  "gap_position_end",
  "reported_aa",
  "reconstructed_aa"
)

required_tail_columns <- c(
  "pair_id",
  "mate",
  "orientation_role",
  "hsp_qstart",
  "hsp_qend",
  "hsp_sstart",
  "hsp_send",
  "gap_facing_nt_start",
  "gap_facing_nt_end",
  "gap_facing_sequence"
)

missing_audit <- setdiff(
  required_audit_columns,
  names(audit_3j)
)

missing_tail <- setdiff(
  required_tail_columns,
  names(raw_tails)
)

if (length(missing_audit) > 0) {
  stop(
    "Missing columns in 4B.3j audit:\n",
    paste(missing_audit, collapse = ", ")
  )
}

if (length(missing_tail) > 0) {
  stop(
    "Missing columns in raw tail dataset:\n",
    paste(missing_tail, collapse = ", ")
  )
}


# ==============================================================
# 7. SELECT THE 43 SAME-LENGTH DISCORDANCE CASES
# ==============================================================

same_length <- audit_3j %>%
  dplyr::filter(
    discordance_class ==
      "same_length_amino_acid_discordance"
  )

cat(
  "Same-length amino-acid discordance cases:",
  nrow(same_length),
  "\n"
)


# ==============================================================
# 8. JOIN BACK TO ORIGINAL TAIL INFORMATION
# ==============================================================

tail_coordinates <- raw_tails %>%
  dplyr::select(
    pair_id,
    mate,
    orientation_role,
    hsp_qstart,
    hsp_qend,
    hsp_sstart,
    hsp_send,
    gap_facing_nt_start,
    gap_facing_nt_end,
    gap_facing_sequence
  ) %>%
  dplyr::distinct(
    pair_id,
    mate,
    orientation_role,
    .keep_all = TRUE
  )

audit_df <- same_length %>%
  dplyr::inner_join(
    tail_coordinates,
    by = c(
      "pair_id",
      "mate",
      "orientation_role"
    )
  )

cat(
  "Successfully joined tail coordinates:",
  nrow(audit_df),
  "\n\n"
)


# ==============================================================
# 9. STANDARD GENETIC CODE
# ==============================================================

codon_map <- c(

  "TTT" = "F",
  "TTC" = "F",
  "TTA" = "L",
  "TTG" = "L",

  "TCT" = "S",
  "TCC" = "S",
  "TCA" = "S",
  "TCG" = "S",

  "TAT" = "Y",
  "TAC" = "Y",
  "TAA" = "*",
  "TAG" = "*",

  "TGT" = "C",
  "TGC" = "C",
  "TGG" = "W",

  "CTT" = "L",
  "CTC" = "L",
  "CTA" = "L",
  "CTG" = "L",

  "CCT" = "P",
  "CCC" = "P",
  "CCA" = "P",
  "CCG" = "P",

  "CAT" = "H",
  "CAC" = "H",
  "CAA" = "Q",
  "CAG" = "Q",

  "CGT" = "R",
  "CGC" = "R",
  "CGA" = "R",
  "CGG" = "R",

  "ATT" = "I",
  "ATC" = "I",
  "ATA" = "I",
  "ATG" = "M",

  "ACT" = "T",
  "ACC" = "T",
  "ACA" = "T",
  "ACG" = "T",

  "AAT" = "N",
  "AAC" = "N",
  "AAA" = "K",
  "AAG" = "K",

  "AGT" = "S",
  "AGC" = "S",
  "AGA" = "R",
  "AGG" = "R",

  "GTT" = "V",
  "GTC" = "V",
  "GTA" = "V",
  "GTG" = "V",

  "GCT" = "A",
  "GCC" = "A",
  "GCA" = "A",
  "GCG" = "A",

  "GAT" = "D",
  "GAC" = "D",
  "GAA" = "E",
  "GAG" = "E",

  "GGT" = "G",
  "GGC" = "G",
  "GGA" = "G",
  "GGG" = "G"
)


# ==============================================================
# 10. CODON TRANSLATION HELPER
# ==============================================================

translate_single_codon <- function(codon) {

  if (is.na(codon)) {
    return("X")
  }

  codon <- toupper(as.character(codon))

  if (nchar(codon) != 3) {
    return("X")
  }

  aa <- unname(codon_map[codon])

  if (length(aa) == 0 || is.na(aa)) {
    return("X")
  }

  aa
}


# ==============================================================
# 11. GET STANDARD CODONS FOR AN AMINO ACID
# ==============================================================

get_standard_codons <- function(aa) {

  if (is.na(aa) || length(aa) == 0) {
    return(character(0))
  }

  aa <- toupper(as.character(aa))

  names(codon_map)[
    unname(codon_map) == aa
  ]
}


# ==============================================================
# 12. MINIMUM CODON HAMMING DISTANCE
# ==============================================================

min_codon_hamming_dist <- function(
  codon,
  target_aa
) {

  if (
    is.na(codon) ||
    is.na(target_aa) ||
    nchar(codon) != 3
  ) {
    return(NA_integer_)
  }

  standard_codons <- get_standard_codons(
    target_aa
  )

  if (length(standard_codons) == 0) {
    return(NA_integer_)
  }

  codon_chars <- strsplit(
    toupper(codon),
    "",
    fixed = TRUE
  )[[1]]

  distances <- vapply(
    standard_codons,
    function(reference_codon) {

      reference_chars <- strsplit(
        reference_codon,
        "",
        fixed = TRUE
      )[[1]]  # FIXED: added missing closing bracket

      sum(
        codon_chars != reference_chars
      )
    },
    integer(1)
  )

  min(distances)
}


# ==============================================================
# 13. EXTRACT CODON FROM GAP-FACING SEQUENCE
# ==============================================================

extract_gap_codon <- function(
  gap_sequence,
  gap_facing_nt_start,
  codon_nt_start,
  aa_index
) {

  if (
    is.na(gap_sequence) ||
    is.na(gap_facing_nt_start) ||
    is.na(codon_nt_start) ||
    is.na(aa_index)
  ) {
    return(NA_character_)
  }

  gap_sequence <- toupper(
    as.character(gap_sequence)
  )

  local_start <- (
    codon_nt_start -
      gap_facing_nt_start +
      1L
  ) +
    ((aa_index - 1L) * 3L)

  local_end <- local_start + 2L

  if (local_start < 1L) {
    return(NA_character_)
  }

  if (local_end > nchar(gap_sequence)) {
    return(NA_character_)
  }

  substr(
    gap_sequence,
    local_start,
    local_end
  )
}


# ==============================================================
# 14. DETERMINE LOCAL FRAME ALTERNATIVES
# ==============================================================

extract_frame_codon <- function(
  gap_sequence,
  local_start,
  offset
) {

  if (
    is.na(gap_sequence) ||
    is.na(local_start)
  ) {
    return(NA_character_)
  }

  start <- local_start + offset
  end <- start + 2L

  if (start < 1L || end > nchar(gap_sequence)) {
    return(NA_character_)
  }

  substr(
    gap_sequence,
    start,
    end
  )
}


# ==============================================================
# 15. PROVENANCE AUDIT
# ==============================================================

provenance_rows <- vector(
  "list",
  nrow(audit_df)
)

for (i in seq_len(nrow(audit_df))) {

  row <- audit_df[i, ]

  reported_aa <- as.character(
    row$reported_aa[[1]]
  )

  reconstructed_aa <- as.character(
    row$reconstructed_aa[[1]]
  )

  gap_sequence <- as.character(
    row$gap_facing_sequence[[1]]
  )

  gap_position_start <- suppressWarnings(
    as.integer(row$gap_position_start[[1]])
  )

  gap_position_end <- suppressWarnings(
    as.integer(row$gap_position_end[[1]])
  )

  gf_start <- suppressWarnings(
    as.integer(row$gap_facing_nt_start[[1]])
  )

  gf_end <- suppressWarnings(
    as.integer(row$gap_facing_nt_end[[1]])
  )

  hsp_sstart <- suppressWarnings(
    as.integer(row$hsp_sstart[[1]])
  )

  hsp_send <- suppressWarnings(
    as.integer(row$hsp_send[[1]])
  )

  hsp_qstart <- suppressWarnings(
    as.integer(row$hsp_qstart[[1]])
  )


  # ------------------------------------------------------------
  # Validate AA strings
  # ------------------------------------------------------------

  if (
    is.na(reported_aa) ||
    is.na(reconstructed_aa) ||
    reported_aa == "" ||
    reconstructed_aa == ""
  ) {
    next
  }

  reported_chars <- strsplit(
    reported_aa,
    "",
    fixed = TRUE
  )[[1]]

  reconstructed_chars <- strsplit(
    reconstructed_aa,
    "",
    fixed = TRUE
  )[[1]]


  # ------------------------------------------------------------
  # Identify all AA mismatches
  # ------------------------------------------------------------

  n_compare <- min(
    length(reported_chars),
    length(reconstructed_chars)
  )

  mismatch_idx <- which(
    reported_chars[seq_len(n_compare)] !=
      reconstructed_chars[seq_len(n_compare)]
  )

  if (length(mismatch_idx) == 0) {
    next
  }


  # ------------------------------------------------------------
  # First discordant position
  # ------------------------------------------------------------

  first_idx <- mismatch_idx[1]

  first_tac102_pos <- if (
    !is.na(gap_position_start)
  ) {
    gap_position_start +
      first_idx -
      1L
  } else {
    NA_integer_
  }

  reported_first_aa <-
    reported_chars[first_idx]

  reconstructed_first_aa <-
    reconstructed_chars[first_idx]


  # ------------------------------------------------------------
  # Determine codon-anchored nucleotide start
  #
  # This follows the same coordinate relationship established
  # during the 4B.3i window check.
  # ------------------------------------------------------------

  forward <- if (
    !is.na(hsp_sstart) &&
    !is.na(hsp_send)
  ) {
    hsp_sstart < hsp_send
  } else {
    TRUE
  }


  # The codon corresponding to TAC102 position
  # gap_position_start is anchored at the first nucleotide
  # of the codon-anchored extension.
  #
  # We derive the nucleotide coordinate of the first codon
  # from the original HSP coordinates.

  if (
    !is.na(hsp_qstart) &&
    !is.na(hsp_sstart) &&
    !is.na(gap_position_start)
  ) {

    if (forward) {

      codon_nt_start <- (
        hsp_sstart +
          (
            gap_position_start -
              hsp_qstart
          ) * 3L
      )

    } else {

      codon_nt_high <- (
        hsp_sstart -
          (
            gap_position_start -
              hsp_qstart
          ) * 3L
      )

      codon_nt_start <-
        codon_nt_high - 2L
    }

  } else {

    codon_nt_start <- NA_integer_
  }


  # ------------------------------------------------------------
  # Extract codon corresponding to first mismatch
  # ------------------------------------------------------------

  reconstructed_codon <- extract_gap_codon(
    gap_sequence = gap_sequence,
    gap_facing_nt_start = gf_start,
    codon_nt_start = codon_nt_start,
    aa_index = first_idx
  )

  reconstructed_codon_aa <-
    translate_single_codon(
      reconstructed_codon
    )


  # ------------------------------------------------------------
  # Calculate local sequence position
  # ------------------------------------------------------------

  local_first_codon_start <- if (
    !is.na(gf_start) &&
    !is.na(codon_nt_start)
  ) {

    (
      codon_nt_start -
        gf_start +
        1L
    ) +
      (
        first_idx -
          1L
      ) * 3L

  } else {

    NA_integer_
  }


  # ------------------------------------------------------------
  # Local +1 / +2 frame alternatives
  # ------------------------------------------------------------

  codon_frame0 <- extract_frame_codon(
    gap_sequence,
    local_first_codon_start,
    0L
  )

  codon_frame1 <- extract_frame_codon(
    gap_sequence,
    local_first_codon_start,
    1L
  )

  codon_frame2 <- extract_frame_codon(
    gap_sequence,
    local_first_codon_start,
    2L
  )

  aa_frame0 <- translate_single_codon(
    codon_frame0
  )

  aa_frame1 <- translate_single_codon(
    codon_frame1
  )

  aa_frame2 <- translate_single_codon(
    codon_frame2
  )


  # ------------------------------------------------------------
  # Does a local boundary shift recover reported AA?
  # ------------------------------------------------------------

  frame_shift_restores_reported <-
    aa_frame1 == reported_first_aa ||
    aa_frame2 == reported_first_aa


  # ------------------------------------------------------------
  # Hamming distance
  #
  # IMPORTANT:
  # This is only a diagnostic measure.
  # It does NOT establish that the biological nucleotide
  # sequence actually contains a mutation.
  # ------------------------------------------------------------

  hamming_distance <-
    min_codon_hamming_dist(
      reconstructed_codon,
      reported_first_aa
    )


  # ------------------------------------------------------------
  # Codon-level interpretation
  # ------------------------------------------------------------

  codon_directly_encodes_reported <-
    reconstructed_codon_aa ==
    reported_first_aa


  codon_directly_encodes_reconstructed <-
    reconstructed_codon_aa ==
    reconstructed_first_aa


  # ------------------------------------------------------------
  # Determine whether nucleotide sequence is available
  # ------------------------------------------------------------

  nucleotide_available <- !is.na(
    reconstructed_codon
  ) &&
    nchar(reconstructed_codon) == 3


  # ------------------------------------------------------------
  # Classification
  # ------------------------------------------------------------

  provenance_class <- dplyr::case_when(

    !nucleotide_available ~
      "E_insufficient_nucleotide_evidence",

    codon_directly_encodes_reported ~
      "A_translation_bookkeeping_mismatch",

    frame_shift_restores_reported ~
      "C_codon_boundary_frame_interpretation",

    hamming_distance == 1L ~
      "D_single_nucleotide_difference_candidate",

    hamming_distance >= 2L ~
      "B_extraction_provenance_nucleotide_mismatch",

    TRUE ~
      "E_ambiguous_evidence"
  )


  # ------------------------------------------------------------
  # Store result
  # ------------------------------------------------------------

  provenance_rows[[i]] <- tibble::tibble(

    pair_id =
      as.character(row$pair_id[[1]]),

    mate =
      as.character(row$mate[[1]]),

    orientation_role =
      as.character(row$orientation_role[[1]]),

    hsp_qstart =
      hsp_qstart,

    hsp_sstart =
      hsp_sstart,

    hsp_send =
      hsp_send,

    strand =
      ifelse(forward, "+", "-"),

    gap_position_start =
      gap_position_start,

    gap_position_end =
      gap_position_end,

    gap_facing_nt_start =
      gf_start,

    gap_facing_nt_end =
      gf_end,

    codon_nt_start =
      codon_nt_start,

    total_aa_length =
      length(reported_chars),

    total_mismatch_positions =
      length(mismatch_idx),

    first_discordant_aa_idx =
      first_idx,

    first_discordant_tac102_pos =
      first_tac102_pos,

    reported_aa =
      reported_first_aa,

    reconstructed_aa =
      reconstructed_first_aa,

    reconstructed_codon =
      reconstructed_codon,

    reconstructed_codon_translation =
      reconstructed_codon_aa,

    codon_directly_encodes_reported =
      codon_directly_encodes_reported,

    codon_directly_encodes_reconstructed =
      codon_directly_encodes_reconstructed,

    frame0_codon =
      codon_frame0,

    frame0_aa =
      aa_frame0,

    frame_plus1_codon =
      codon_frame1,

    frame_plus1_aa =
      aa_frame1,

    frame_plus2_codon =
      codon_frame2,

    frame_plus2_aa =
      aa_frame2,

    frame_shift_restores_reported =
      frame_shift_restores_reported,

    hamming_dist_to_reported =
      hamming_distance,

    nucleotide_evidence_available =
      nucleotide_available,

    provenance_classification =
      provenance_class
  )
}


# ==============================================================
# 16. COMBINE RESULTS
# ==============================================================

provenance_df <- dplyr::bind_rows(
  provenance_rows
)


if (nrow(provenance_df) == 0) {
  stop(
    "No provenance records were generated."
  )
}


# ==============================================================
# 17. CLASSIFICATION SUMMARY
# ==============================================================
#
# Explicitly use dplyr::count().
#
# This avoids the previous failure where matrixStats::count()
# was dispatched instead of dplyr::count().
# ==============================================================

class_summary <- provenance_df %>%
  dplyr::count(
    provenance_classification,
    name = "case_count"
  ) %>%
  dplyr::mutate(
    percentage = round(
      100 *
        case_count /
        sum(case_count),
      2
    )
  ) %>%
  dplyr::arrange(
    dplyr::desc(case_count)
  )


# ==============================================================
# 18. FIRST DISCORDANT POSITION DISTRIBUTION
# ==============================================================

pos_distribution <- provenance_df %>%
  dplyr::count(
    first_discordant_tac102_pos,
    orientation_role,
    reported_aa,
    reconstructed_aa,
    name = "case_count"
  ) %>%
  dplyr::arrange(
    dplyr::desc(case_count)
  )


# ==============================================================
# 19. CODON DISCORDANCE SUMMARY
# ==============================================================

codon_summary <- provenance_df %>%
  dplyr::count(
    reported_aa,
    reconstructed_aa,
    reconstructed_codon,
    reconstructed_codon_translation,
    hamming_dist_to_reported,
    provenance_classification,
    name = "case_count"
  ) %>%
  dplyr::arrange(
    dplyr::desc(case_count)
  )


# ==============================================================
# 20. FRAME-SHIFT SUMMARY
# ==============================================================

frame_summary <- provenance_df %>%
  dplyr::count(
    frame_shift_restores_reported,
    name = "case_count"
  ) %>%
  dplyr::mutate(
    percentage = round(
      100 *
        case_count /
        sum(case_count),
      2
    )
  )


# ==============================================================
# 21. PRINT MAIN RESULTS
# ==============================================================

cat("\n")
cat("==============================================================\n")
cat(" POSITION-BY-POSITION PROVENANCE AUDIT SUMMARY\n")
cat("==============================================================\n\n")

print(
  class_summary,
  n = Inf
)


cat("\n")
cat("==============================================================\n")
cat(" FIRST DISCORDANT POSITION DISTRIBUTION\n")
cat("==============================================================\n\n")

print(
  pos_distribution,
  n = 30
)


cat("\n")
cat("==============================================================\n")
cat(" CODON DISCORDANCE BY AMINO-ACID PAIR\n")
cat("==============================================================\n\n")

print(
  codon_summary,
  n = 30
)


cat("\n")
cat("==============================================================\n")
cat(" LOCAL FRAME-SHIFT SUMMARY\n")
cat("==============================================================\n\n")

print(
  frame_summary,
  n = Inf
)


# ==============================================================
# 22. ADDITIONAL DIAGNOSTIC COUNTS
# ==============================================================

cat("\n")
cat("==============================================================\n")
cat(" ADDITIONAL DIAGNOSTICS\n")
cat("==============================================================\n\n")

cat(
  "Audited discordant tail records:",
  nrow(provenance_df),
  "\n"
)

cat(
  "Records with usable reconstructed codon:",
  sum(
    provenance_df$nucleotide_evidence_available,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Records where reconstructed codon directly encodes reported AA:",
  sum(
    provenance_df$codon_directly_encodes_reported,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Records where +1/+2 local shift restores reported AA:",
  sum(
    provenance_df$frame_shift_restores_reported,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Records with minimum codon Hamming distance = 1:",
  sum(
    provenance_df$hamming_dist_to_reported == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Records with minimum codon Hamming distance >= 2:",
  sum(
    provenance_df$hamming_dist_to_reported >= 2,
    na.rm = TRUE
  ),
  "\n"
)


# ==============================================================
# 23. SAVE DETAILED OUTPUT
# ==============================================================

readr::write_csv(
  provenance_df,
  out_detail
)


# ==============================================================
# 24. SAVE CLASSIFICATION SUMMARY
# ==============================================================

readr::write_csv(
  class_summary,
  out_class
)


# ==============================================================
# 25. WRITE TEXT SUMMARY
# ==============================================================

summary_text <- c(

  "Step 4B.3k: TAC102 Targeted Nucleotide Provenance Audit",

  "==============================================================",

  "",

  paste(
    "Total audited same-length discordant cases:",
    nrow(provenance_df)
  ),

  "",

  "Classification Breakdown:",

  capture.output(
    print(
      class_summary,
      n = Inf
    )
  ),

  "",

  "Frame-Shift Summary:",

  capture.output(
    print(
      frame_summary,
      n = Inf
    )
  ),

  "",

  "Interpretation:",

  paste(
    "This audit evaluates the nucleotide provenance of the",
    "same-length amino-acid discordances identified in Step",
    "4B.3j."
  ),

  paste(
    "The reconstructed codon is derived from the stored",
    "gap-facing nucleotide sequence using the codon-anchored",
    "coordinate system established during Step 4B.3i."
  ),

  paste(
    "Amino-acid discordance alone is not interpreted as evidence",
    "of nucleotide variation."
  ),

  paste(
    "A Hamming distance of one identifies only a candidate",
    "single-nucleotide explanation and does not independently",
    "establish a biological SNP."
  ),

  "",

  "Outputs:",

  out_detail,

  out_class,

  out_sum
)


writeLines(
  summary_text,
  out_sum
)


# ==============================================================
# 26. FINAL STATUS
# ==============================================================

cat("\n")
cat("==============================================================\n")
cat(" OUTPUTS\n")
cat("==============================================================\n\n")

cat(
  out_detail,
  "\n"
)

cat(
  out_class,
  "\n"
)

cat(
  out_sum,
  "\n"
)

cat("\n")
cat(
  "[COMPLETE] Step 4B.3k nucleotide provenance audit finished.\n"
)