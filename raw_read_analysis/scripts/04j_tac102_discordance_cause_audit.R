library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)

cat("==============================================================\n")
cat(" Step 4B.3j: TAC102 Direct Cause-of-Discordance Audit\n")
cat("==============================================================\n\n")


# ==============================================================
# 1. PATHS
# ==============================================================

csv_tails <- "raw_read_analysis/reports/tac102_41pairs_unaligned_tails.csv"

out_detail <- "raw_read_analysis/reports/tac102_4B3j_discordance_detailed_audit.csv"
out_class  <- "raw_read_analysis/reports/tac102_4B3j_discordance_classification_summary.csv"
out_sum    <- "raw_read_analysis/reports/tac102_4B3j_discordance_cause_summary.txt"


# ==============================================================
# 2. LOAD DATA
# ==============================================================

if (!file.exists(csv_tails)) {
  stop("Missing input file: ", csv_tails)
}

tails_df <- read.csv(
  csv_tails,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Input rows:", nrow(tails_df), "\n")
cat("Input columns:", ncol(tails_df), "\n\n")


# ==============================================================
# 3. CODON TRANSLATION MAP
# ==============================================================

codon_map <- c(
  "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
  "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
  "TAT"="Y", "TAC"="Y",
  "TAA"="*", "TAG"="*", "TGA"="*",
  "TGT"="C", "TGC"="C", "TGG"="W",

  "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
  "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
  "CAT"="H", "CAC"="H",
  "CAA"="Q", "CAG"="Q",
  "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",

  "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
  "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
  "AAT"="N", "AAC"="N",
  "AAA"="K", "AAG"="K",
  "AGT"="S", "AGC"="S",
  "AGA"="R", "AGG"="R",

  "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
  "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
  "GAT"="D", "GAC"="D",
  "GAA"="E", "GAG"="E",
  "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
)


translate_codons <- function(nt) {

  if (is.na(nt) || nchar(nt) < 3) {
    return("")
  }

  nt <- toupper(nt)

  if (nchar(nt) %% 3 != 0) {
    return("")
  }

  codons <- substring(
    nt,
    seq(1, nchar(nt), by = 3),
    seq(3, nchar(nt), by = 3)
  )

  aa <- unname(codon_map[codons])

  aa[is.na(aa)] <- "X"

  paste(aa, collapse = "")
}


# ==============================================================
# 4. HELPER: EXTRACT EXPECTED CODON WINDOW
# ==============================================================

extract_expected_window <- function(row) {

  gf_start <- suppressWarnings(as.integer(row[["gap_facing_nt_start"]]))
  gf_end   <- suppressWarnings(as.integer(row[["gap_facing_nt_end"]]))

  p_start <- suppressWarnings(as.integer(row[["gap_position_start"]]))
  p_end   <- suppressWarnings(as.integer(row[["gap_position_end"]]))

  seq_gf <- row[["gap_facing_sequence"]]

  if (
    is.na(gf_start) ||
    is.na(gf_end) ||
    is.na(p_start) ||
    is.na(p_end) ||
    is.na(seq_gf) ||
    seq_gf == ""
  ) {
    return(tibble::tibble(
      expected_nt = NA_character_,
      expected_aa = NA_character_,
      expected_nt_start = NA_integer_,
      expected_nt_end = NA_integer_,
      expected_length_bp = NA_integer_,
      extraction_status = "missing_coordinates_or_sequence"
    ))
  }

  # ------------------------------------------------------------
  # Reconstruct codon-anchored window exactly as in 04i
  # ------------------------------------------------------------

  qstart_raw <- suppressWarnings(as.integer(row[["hsp_qstart"]]))
  sstart_raw <- suppressWarnings(as.integer(row[["hsp_sstart"]]))
  send_raw   <- suppressWarnings(as.integer(row[["hsp_send"]]))

  if (
    is.na(qstart_raw) ||
    is.na(sstart_raw) ||
    is.na(send_raw)
  ) {
    return(tibble::tibble(
      expected_nt = NA_character_,
      expected_aa = NA_character_,
      expected_nt_start = NA_integer_,
      expected_nt_end = NA_integer_,
      expected_length_bp = NA_integer_,
      extraction_status = "missing_HSP_coordinates"
    ))
  }

  forward <- sstart_raw < send_raw

  if (forward) {

    nt_start <- sstart_raw +
      (p_start - qstart_raw) * 3L

    nt_end <- sstart_raw +
      (p_end - qstart_raw) * 3L + 2L

  } else {

    nt_end <- sstart_raw -
      (p_start - qstart_raw) * 3L

    nt_start <- sstart_raw -
      (p_end - qstart_raw) * 3L - 2L
  }

  # ------------------------------------------------------------
  # Convert read-relative coordinates to positions inside
  # gap_facing_sequence
  # ------------------------------------------------------------

  local_start <- nt_start - gf_start + 1L
  local_end   <- nt_end   - gf_start + 1L

  seq_len <- nchar(seq_gf)

  if (
    local_start < 1 ||
    local_end > seq_len ||
    local_start > local_end
  ) {
    return(tibble::tibble(
      expected_nt = NA_character_,
      expected_aa = NA_character_,
      expected_nt_start = nt_start,
      expected_nt_end = nt_end,
      expected_length_bp = max(
        0L,
        nt_end - nt_start + 1L
      ),
      extraction_status = "codon_window_outside_gap_facing_sequence"
    ))
  }

  expected_nt <- substr(
    seq_gf,
    local_start,
    local_end
  )

  if (nchar(expected_nt) %% 3 != 0) {
    return(tibble::tibble(
      expected_nt = expected_nt,
      expected_aa = NA_character_,
      expected_nt_start = nt_start,
      expected_nt_end = nt_end,
      expected_length_bp = nchar(expected_nt),
      extraction_status = "extracted_window_not_multiple_of_3"
    ))
  }

  expected_aa <- translate_codons(expected_nt)

  tibble::tibble(
    expected_nt = expected_nt,
    expected_aa = expected_aa,
    expected_nt_start = nt_start,
    expected_nt_end = nt_end,
    expected_length_bp = nchar(expected_nt),
    extraction_status = "success"
  )
}


# ==============================================================
# 5. SELECT USABLE TAILS
# ==============================================================

usable <- tails_df %>%
  dplyr::filter(
    !is.na(gap_facing_sequence),
    gap_facing_sequence != "",
    !is.na(gap_position_start),
    !is.na(gap_position_end),
    !is.na(translated_candidate),
    translated_candidate != ""
  )

cat("Usable translated tail observations:", nrow(usable), "\n\n")


# ==============================================================
# 6. INDEPENDENTLY RECONSTRUCT EACH TAIL
# ==============================================================

cat("--- Independently reconstructing nucleotide-to-protein evidence ---\n\n")

audit_rows <- lapply(seq_len(nrow(usable)), function(i) {

  row <- usable[i, ]

  reconstructed <- extract_expected_window(row)

  reported_aa <- toupper(row$translated_candidate)

  expected_aa <- reconstructed$expected_aa

  aa_identical <- (
    !is.na(expected_aa) &&
    expected_aa != "" &&
    expected_aa == reported_aa
  )

  aa_same_length <- (
    !is.na(expected_aa) &&
    nchar(expected_aa) == nchar(reported_aa)
  )

  nt_identical_to_expected <- FALSE

  if (
    reconstructed$extraction_status == "success" &&
    !is.na(reconstructed$expected_nt)
  ) {
    nt_identical_to_expected <- TRUE
  }

  mismatch_count <- NA_integer_

  if (aa_same_length) {
    mismatch_count <- sum(
      strsplit(expected_aa, "", fixed = TRUE)[[1]] !=
        strsplit(reported_aa, "", fixed = TRUE)[[1]]
    )
  }

  classification <- dplyr::case_when(
    reconstructed$extraction_status != "success" ~
      reconstructed$extraction_status,

    aa_identical ~
      "direct_translation_reproduces_reported_candidate",

    !aa_same_length ~
      "translation_length_discordance",

    aa_same_length && mismatch_count > 0 ~
      "same_length_amino_acid_discordance",

    TRUE ~
      "unclassified"
  )

  tibble::tibble(
    pair_id = row$pair_id,
    mate = row$mate,
    orientation_role = row$orientation_role,

    gap_position_start = row$gap_position_start,
    gap_position_end = row$gap_position_end,

    reported_aa = reported_aa,
    reconstructed_aa = expected_aa,

    reported_aa_length = nchar(reported_aa),
    reconstructed_aa_length =
      ifelse(is.na(expected_aa), NA_integer_, nchar(expected_aa)),

    aa_mismatch_count = mismatch_count,

    codon_nt_start = reconstructed$expected_nt_start,
    codon_nt_end   = reconstructed$expected_nt_end,

    expected_nt = reconstructed$expected_nt,

    extraction_status = reconstructed$extraction_status,

    direct_nt_translation_match = aa_identical,

    discordance_class = classification
  )
})

audit_df <- dplyr::bind_rows(audit_rows)


# ==============================================================
# 7. SUMMARY OF DIRECT RECONSTRUCTION
# ==============================================================

cat("=== DIRECT RECONSTRUCTION SUMMARY ===\n\n")

class_summary <- audit_df %>%
  dplyr::count(discordance_class, name = "case_count") %>%
  dplyr::mutate(
    percentage = round(
      100 * case_count / nrow(audit_df),
      2
    )
  ) %>%
  dplyr::arrange(desc(case_count))

print(class_summary)

cat("\n")


# ==============================================================
# 8. REPORT EXACT MATCHES
# ==============================================================

n_exact <- sum(
  audit_df$discordance_class ==
    "direct_translation_reproduces_reported_candidate",
  na.rm = TRUE
)

n_length <- sum(
  audit_df$discordance_class ==
    "translation_length_discordance",
  na.rm = TRUE
)

n_aa_diff <- sum(
  audit_df$discordance_class ==
    "same_length_amino_acid_discordance",
  na.rm = TRUE
)

cat("Exact direct translations:", n_exact, "\n")
cat("Translation-length discordances:", n_length, "\n")
cat("Same-length AA discordances:", n_aa_diff, "\n\n")


# ==============================================================
# 9. CHARACTERIZE SAME-LENGTH AA DISCORDANCE
# ==============================================================

same_length_df <- audit_df %>%
  dplyr::filter(
    discordance_class ==
      "same_length_amino_acid_discordance"
  )

cat("=== SAME-LENGTH AA DISCORDANCE ===\n\n")

if (nrow(same_length_df) > 0) {

  mismatch_summary <- same_length_df %>%
    dplyr::select(
      pair_id,
      mate,
      orientation_role,
      gap_position_start,
      gap_position_end,
      reported_aa,
      reconstructed_aa,
      aa_mismatch_count
    )

  print(mismatch_summary, row.names = FALSE)

} else {

  cat("No same-length amino-acid discordances detected.\n")
}

cat("\n")


# ==============================================================
# 10. POSITION-LEVEL DISCORDANCE ANALYSIS
# ==============================================================

cat("=== POSITION-LEVEL AMINO-ACID DISCORDANCE ===\n\n")

position_rows <- list()

if (nrow(same_length_df) > 0) {

  for (i in seq_len(nrow(same_length_df))) {

    row <- same_length_df[i, ]

    reported_chars <- strsplit(row$reported_aa, "", fixed = TRUE)[[1]]
    reconstructed_chars <- strsplit(row$reconstructed_aa, "", fixed = TRUE)[[1]]

    mismatch_positions <- which(reported_chars != reconstructed_chars)

    if (length(mismatch_positions) > 0) {

      for (j in mismatch_positions) {

        tac102_abs_pos <- as.integer(row$gap_position_start) + j - 1L

        position_rows[[length(position_rows) + 1L]] <- tibble::tibble(
          pair_id = row$pair_id,
          mate = row$mate,
          orientation_role = row$orientation_role,
          local_aa_index = j,
          tac102_pos = tac102_abs_pos,
          reported_aa = reported_chars[j],
          reconstructed_aa = reconstructed_chars[j]
        )
      }
    }
  }
}

position_df <- dplyr::bind_rows(position_rows)

if (nrow(position_df) > 0) {

  position_discordance_summary <- position_df %>%
    dplyr::count(
      reported_aa,
      reconstructed_aa,
      name = "discordance_count"
    ) %>%
    dplyr::arrange(desc(discordance_count))

  print(position_discordance_summary)

  cat("\n=== TOP AFFECTED TAC102 POSITIONS ===\n\n")

  top_positions <- position_df %>%
    dplyr::count(
      tac102_pos,
      orientation_role,
      reported_aa,
      reconstructed_aa,
      name = "position_count"
    ) %>%
    dplyr::arrange(desc(position_count))

  print(top_positions, n = 30)

} else {

  cat("No positional amino-acid discordances detected.\n")
}

cat("\n")


# ==============================================================
# 11. LOOK FOR LOW-COMPLEXITY SIGNATURES
# ==============================================================

cat("=== LOW-COMPLEXITY SIGNATURE CHECK ===\n\n")

if (nrow(same_length_df) > 0) {

  low_complexity_df <- same_length_df %>%
    dplyr::mutate(
      reported_K_fraction = str_count(reported_aa, "K") / nchar(reported_aa),
      reconstructed_K_fraction = str_count(reconstructed_aa, "K") / nchar(reconstructed_aa),
      reported_E_fraction = str_count(reported_aa, "E") / nchar(reported_aa),
      reconstructed_E_fraction = str_count(reconstructed_aa, "E") / nchar(reconstructed_aa)
    ) %>%
    dplyr::select(
      pair_id,
      orientation_role,
      reported_K_fraction,
      reconstructed_K_fraction,
      reported_E_fraction,
      reconstructed_E_fraction
    )

  print(low_complexity_df, row.names = FALSE)

} else {

  cat("No same-length discordances available for low-complexity analysis.\n")
}

cat("\n")


# ==============================================================
# 12. SAVE OUTPUTS
# ==============================================================

readr::write_csv(
  audit_df,
  out_detail
)

readr::write_csv(
  class_summary,
  out_class
)


# ==============================================================
# 13. SUMMARY REPORT
# ==============================================================

summary_text <- c(

  "Step 4B.3j: TAC102 Direct Cause-of-Discordance Audit",
  "==============================================================",
  "",

  paste(
    "Total usable translated tail observations:",
    nrow(audit_df)
  ),

  paste(
    "Direct translations reproducing reported candidate:",
    n_exact
  ),

  paste(
    "Translation-length discordances:",
    n_length
  ),

  paste(
    "Same-length amino-acid discordances:",
    n_aa_diff
  ),

  "",

  "Discordance classification:",

  capture.output(
    print(class_summary)
  ),

  "",

  "Interpretation:",
  "",

  paste(
    "This audit independently reconstructs the nucleotide window",
    "corresponding to each reported gap-position range from the",
    "stored gap-facing sequence."
  ),

  paste(
    "The independently translated nucleotide sequence is then",
    "compared directly with translated_candidate."
  ),

  paste(
    "This distinguishes coordinate/extraction problems from",
    "true nucleotide-supported amino-acid discordance."
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
# 14. FINAL MESSAGE
# ==============================================================

cat("Saved reports:\n")
cat("  ", out_detail, "\n")
cat("  ", out_class, "\n")
cat("  ", out_sum, "\n\n")

cat("[COMPLETE] Step 4B.3j discordance-cause audit finished.\n")