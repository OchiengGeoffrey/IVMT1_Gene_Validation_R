library(dplyr)
library(readr)
library(tidyr)
library(stringr)

cat("==============================================================\n")
cat(" Step 4B.3h: Direct Nucleotide-to-Protein Frame & Consistency Audit\n")
cat("==============================================================\n\n")

# ==============================================================
# PATH DEFINITIONS
# ==============================================================

csv_conflicts <- "raw_read_analysis/reports/tac102_4B3f_pair_conflicts.csv"
csv_tails     <- "raw_read_analysis/reports/tac102_41pairs_unaligned_tails.csv"
csv_04g_audit <- "raw_read_analysis/reports/tac102_4B3g_conflict_detailed_audit.csv"

out_55_audit <- "raw_read_analysis/reports/tac102_4B3h_55_anomaly_audit.csv"
out_4_snps   <- "raw_read_analysis/reports/tac102_4B3h_4_nonsynonymous_candidates.csv"
out_summary  <- "raw_read_analysis/reports/tac102_4B3h_frame_audit_summary.txt"


# ==============================================================
# STANDARD CODON TRANSLATION MAP
# ==============================================================

codon_map <- c(
  "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
  "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
  "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*",
  "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",

  "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
  "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
  "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
  "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",

  "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
  "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
  "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
  "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",

  "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
  "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
  "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
  "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
)


# ==============================================================
# TRANSLATION FUNCTIONS
# ==============================================================

translate_nt <- function(nt_seq, frame = 1) {

  if (is.na(nt_seq) || nchar(nt_seq) < 3) {
    return("")
  }

  nt_seq <- toupper(nt_seq)

  seq_sub <- substr(
    nt_seq,
    frame,
    nchar(nt_seq)
  )

  usable_length <- floor(nchar(seq_sub) / 3) * 3

  if (usable_length < 3) {
    return("")
  }

  seq_sub <- substr(
    seq_sub,
    1,
    usable_length
  )

  codon_starts <- seq(
    1,
    usable_length - 2,
    by = 3
  )

  codon_ends <- seq(
    3,
    usable_length,
    by = 3
  )

  codons <- substring(
    seq_sub,
    codon_starts,
    codon_ends
  )

  aas <- unname(
    codon_map[codons]
  )

  aas[is.na(aas)] <- "X"

  paste(
    aas,
    collapse = ""
  )
}


rev_comp <- function(nt_seq) {

  if (is.na(nt_seq)) {
    return(NA_character_)
  }

  nt_seq <- toupper(nt_seq)

  comp <- chartr(
    "ATCG",
    "TAGC",
    nt_seq
  )

  chars <- strsplit(
    comp,
    "",
    fixed = TRUE
  )[[1]]

  paste(
    rev(chars),
    collapse = ""
  )
}


# ==============================================================
# SAFE FIRST-ROW EXTRACTION
# ==============================================================
#
# Important:
# We deliberately use base-R indexing rather than dplyr::slice().
# In this environment BiocGenerics/IRanges can mask slice().
# ==============================================================

get_tail <- function(df, pair_id_value, orientation_value) {

  hits <- df[
    df$pair_id == pair_id_value &
      df$orientation_role == orientation_value,
    ,
    drop = FALSE
  ]

  if (nrow(hits) == 0) {
    return(NULL)
  }

  hits[1, , drop = FALSE]
}


# ==============================================================
# LOAD DATA
# ==============================================================

if (!file.exists(csv_conflicts)) {
  stop("Missing input file: ", csv_conflicts)
}

if (!file.exists(csv_tails)) {
  stop("Missing input file: ", csv_tails)
}

if (!file.exists(csv_04g_audit)) {
  stop("Missing input file: ", csv_04g_audit)
}

conflicts_df <- read_csv(
  csv_conflicts,
  show_col_types = FALSE
)

tails_df <- read_csv(
  csv_tails,
  show_col_types = FALSE
)

audit_04g_df <- read_csv(
  csv_04g_audit,
  show_col_types = FALSE
)


# ==============================================================
# INPUT SUMMARY
# ==============================================================

cat("Input rows:\n")
cat("Internal conflicts:", nrow(conflicts_df), "\n")
cat("Tail observations:", nrow(tails_df), "\n")
cat("04B.3g audit rows:", nrow(audit_04g_df), "\n\n")


# ==============================================================
# IDENTIFY NUCLEOTIDE COLUMN
# ==============================================================

nt_candidates <- c(
  "gap_facing_sequence",
  "unaligned_seq"
)

nt_column <- nt_candidates[
  nt_candidates %in% names(tails_df)
][1]

if (is.na(nt_column)) {
  stop(
    "Could not identify nucleotide sequence column. ",
    "Expected one of: ",
    paste(nt_candidates, collapse = ", ")
  )
}

cat(
  "Nucleotide sequence column:",
  nt_column,
  "\n\n"
)


# ==============================================================
# PREPARE USABLE TAIL DATA
# ==============================================================

tails_usable <- tails_df %>%
  filter(
    !is.na(.data[[nt_column]]),
    .data[[nt_column]] != "",
    !is.na(translated_candidate),
    translated_candidate != "",
    !is.na(gap_position_start),
    !is.na(gap_position_end)
  )

cat(
  "Usable tail observations:",
  nrow(tails_usable),
  "\n"
)

cat(
  "N-flank observations:",
  sum(tails_usable$orientation_role == "N"),
  "\n"
)

cat(
  "C-flank observations:",
  sum(tails_usable$orientation_role == "C"),
  "\n\n"
)


# ==============================================================
# DIRECT FRAME AUDIT OF ALL USABLE TAIL SEQUENCES
# ==============================================================

cat("--- Direct frame audit of all usable tail sequences ---\n\n")

frame_audit <- bind_rows(
  lapply(seq_len(nrow(tails_usable)), function(i) {

    row <- tails_usable[i, ]

    nt <- as.character(
      row[[nt_column]]
    )

    reported_aa <- as.character(
      row$translated_candidate
    )

    forward_frames <- sapply(
      1:3,
      function(f) translate_nt(nt, f)
    )

    reverse_frames <- sapply(
      1:3,
      function(f) {
        translate_nt(
          rev_comp(nt),
          f
        )
      }
    )

    forward_match <- which(
      forward_frames == reported_aa
    )

    reverse_match <- which(
      reverse_frames == reported_aa
    )

    if (length(forward_match) > 0) {

      matched_frame <- forward_match[1]

      result_class <- ifelse(
        matched_frame == 1,
        "frame_1_match",
        "non_frame_1_forward"
      )

    } else if (length(reverse_match) > 0) {

      matched_frame <- reverse_match[1]

      result_class <- "reverse_complement_required"

    } else {

      matched_frame <- NA_integer_

      result_class <- "translation_not_reproduced"
    }

    tibble(
      pair_id = as.character(row$pair_id),
      orientation_role = as.character(row$orientation_role),
      gap_position_start = as.numeric(row$gap_position_start),
      gap_position_end = as.numeric(row$gap_position_end),
      reported_translation = reported_aa,
      nucleotide_length = nchar(nt),
      forward_frame_1 = forward_frames[1],
      forward_frame_2 = forward_frames[2],
      forward_frame_3 = forward_frames[3],
      reverse_frame_1 = reverse_frames[1],
      reverse_frame_2 = reverse_frames[2],
      reverse_frame_3 = reverse_frames[3],
      matched_frame = matched_frame,
      direct_frame_class = result_class
    )
  })
)

frame_1_matches <- sum(
  frame_audit$direct_frame_class == "frame_1_match",
  na.rm = TRUE
)

non_frame_1 <- sum(
  frame_audit$direct_frame_class == "non_frame_1_forward",
  na.rm = TRUE
)

reverse_required <- sum(
  frame_audit$direct_frame_class == "reverse_complement_required",
  na.rm = TRUE
)

not_reproduced <- sum(
  frame_audit$direct_frame_class == "translation_not_reproduced",
  na.rm = TRUE
)

cat(
  "Tail sequences matching reported translation:",
  frame_1_matches,
  "\n"
)

cat(
  "Tail sequences requiring non-frame-1 forward translation:",
  non_frame_1,
  "\n"
)

cat(
  "Tail sequences requiring reverse complement:",
  reverse_required,
  "\n"
)

cat(
  "Tail sequences where reported translation cannot be reproduced:",
  not_reproduced,
  "\n\n"
)


# ==============================================================
# PART A
# AUDIT THE 55 IDENTICAL-NUCLEOTIDE / DIFFERENT-AA ANOMALIES
# ==============================================================

cat("--- Part A: Auditing 55 identical-nucleotide / different-AA anomalies ---\n\n")

anomalies_55 <- audit_04g_df %>%
  filter(
    conflict_class ==
      "nucleotide_identical_translation_discrepancy"
  )

cat(
  "Found",
  nrow(anomalies_55),
  "cases.\n\n"
)


# ==============================================================
# AUDIT EACH 55 ANOMALY
# ==============================================================

audit_55_results <- bind_rows(

  lapply(
    seq_len(nrow(anomalies_55)),
    function(i) {

      row <- anomalies_55[i, ]

      pair_id_value <- as.character(
        row$pair_id
      )

      tail_N <- get_tail(
        tails_df,
        pair_id_value,
        "N"
      )

      tail_C <- get_tail(
        tails_df,
        pair_id_value,
        "C"
      )

      if (is.null(tail_N) || is.null(tail_C)) {

        return(
          tibble(
            pair_id = pair_id_value,
            tac102_pos = row$tac102_pos,
            aa_N = row$aa_N,
            aa_C = row$aa_C,
            extracted_codon = row$codon_N,
            frame_N_matched = "Missing",
            frame_C_matched = "Missing",
            diagnostic_cause =
              "Missing N/C tail record"
          )
        )
      }

      nt_N <- as.character(
        tail_N[[nt_column]][1]
      )

      nt_C <- as.character(
        tail_C[[nt_column]][1]
      )

      reported_N <- as.character(
        tail_N$translated_candidate[1]
      )

      reported_C <- as.character(
        tail_C$translated_candidate[1]
      )

      f_N <- sapply(
        1:3,
        function(f) translate_nt(nt_N, f)
      )

      f_C <- sapply(
        1:3,
        function(f) translate_nt(nt_C, f)
      )

      rc_N <- sapply(
        1:3,
        function(f) {
          translate_nt(
            rev_comp(nt_N),
            f
          )
        }
      )

      rc_C <- sapply(
        1:3,
        function(f) {
          translate_nt(
            rev_comp(nt_C),
            f
          )
        }
      )

      match_frame_N <- which(
        f_N == reported_N
      )

      match_frame_C <- which(
        f_C == reported_C
      )

      rc_match_N <- which(
        rc_N == reported_N
      )

      rc_match_C <- which(
        rc_C == reported_C
      )

      if (
        length(match_frame_N) > 0 &&
        match_frame_N[1] == 1 &&
        length(match_frame_C) > 0 &&
        match_frame_C[1] == 1
      ) {

        diagnostic_cause <-
          "Both reported translations reproduce in frame 1; conflict is a translation bookkeeping discrepancy"

      } else if (
        length(match_frame_N) > 0 &&
        match_frame_N[1] != 1
      ) {

        diagnostic_cause <-
          paste0(
            "N-tail frame offset (+",
            match_frame_N[1] - 1,
            " bp)"
          )

      } else if (
        length(match_frame_C) > 0 &&
        match_frame_C[1] != 1
      ) {

        diagnostic_cause <-
          paste0(
            "C-tail frame offset (+",
            match_frame_C[1] - 1,
            " bp)"
          )

      } else if (
        length(rc_match_N) > 0 ||
        length(rc_match_C) > 0
      ) {

        diagnostic_cause <-
          "Strand orientation/reverse-complement required"

      } else {

        diagnostic_cause <-
          "Identical extracted nucleotide triplets but reported translations cannot be independently reproduced"
      }

      tibble(
        pair_id = pair_id_value,
        tac102_pos = row$tac102_pos,
        aa_N = row$aa_N,
        aa_C = row$aa_C,
        extracted_codon = row$codon_N,
        frame_N_matched = ifelse(
          length(match_frame_N) > 0,
          paste0("+", match_frame_N[1]),
          "RC/None"
        ),
        frame_C_matched = ifelse(
          length(match_frame_C) > 0,
          paste0("+", match_frame_C[1]),
          "RC/None"
        ),
        diagnostic_cause = diagnostic_cause
      )
    }
  )
)


# ==============================================================
# 55-ANOMALY SUMMARY
# ==============================================================

cat("=== 55-ANOMALY DIAGNOSTIC SUMMARY ===\n\n")

if (nrow(audit_55_results) > 0) {

  anomaly_summary <- audit_55_results %>%
    group_by(diagnostic_cause) %>%
    summarise(
      case_count = n(),
      percentage = round(
        100 * case_count / nrow(audit_55_results),
        2
      ),
      .groups = "drop"
    ) %>%
    arrange(
      desc(case_count)
    )

  print(anomaly_summary)

} else {

  anomaly_summary <- tibble(
    diagnostic_cause = character(),
    case_count = integer(),
    percentage = numeric()
  )

  cat("No 55-anomaly cases were present.\n")
}

cat("\n")


# ==============================================================
# PART B
# INSPECT CANDIDATE NONSYNONYMOUS NUCLEOTIDE DIFFERENCES
# ==============================================================

cat("--- Part B: Inspecting candidate nonsynonymous nucleotide differences ---\n\n")

snps_4 <- audit_04g_df %>%
  filter(
    conflict_class ==
      "nucleotide_different"
  )

cat(
  "Candidate nonsynonymous nucleotide differences:",
  nrow(snps_4),
  "\n\n"
)


# ==============================================================
# IMPORTANT:
# Always create an empty SNP tibble WITH THE EXPECTED COLUMNS.
# This prevents select() from failing when there are zero candidates.
# ==============================================================

snp_template <- tibble(
  pair_id = character(),
  tac102_pos = integer(),
  aa_N = character(),
  aa_C = character(),
  codon_N = character(),
  codon_C = character(),
  context_NT_N = character(),
  context_NT_C = character(),
  is_valid_snp = logical()
)


if (nrow(snps_4) == 0) {

  snp_results <- snp_template

} else {

  snp_results <- bind_rows(

    lapply(
      seq_len(nrow(snps_4)),
      function(i) {

        row <- snps_4[i, ]

        pair_id_value <- as.character(
          row$pair_id
        )

        tail_N <- get_tail(
          tails_df,
          pair_id_value,
          "N"
        )

        tail_C <- get_tail(
          tails_df,
          pair_id_value,
          "C"
        )

        if (is.null(tail_N) || is.null(tail_C)) {

          return(
            tibble(
              pair_id = pair_id_value,
              tac102_pos = row$tac102_pos,
              aa_N = row$aa_N,
              aa_C = row$aa_C,
              codon_N = row$codon_N,
              codon_C = row$codon_C,
              context_NT_N = NA_character_,
              context_NT_C = NA_character_,
              is_valid_snp = FALSE
            )
          )
        }

        pos <- as.numeric(
          row$tac102_pos
        )

        start_N <- as.numeric(
          tail_N$gap_position_start[1]
        )

        start_C <- as.numeric(
          tail_C$gap_position_start[1]
        )

        idx_N <- pos - start_N + 1
        idx_C <- pos - start_C + 1

        nt_start_N <- (idx_N - 1) * 3 + 1
        nt_start_C <- (idx_C - 1) * 3 + 1

        nt_N <- as.character(
          tail_N[[nt_column]][1]
        )

        nt_C <- as.character(
          tail_C[[nt_column]][1]
        )

        ctx_N <- substr(
          nt_N,
          max(1, nt_start_N - 6),
          min(
            nchar(nt_N),
            nt_start_N + 8
          )
        )

        ctx_C <- substr(
          nt_C,
          max(1, nt_start_C - 6),
          min(
            nchar(nt_C),
            nt_start_C + 8
          )
        )

        valid_snp <-
          !is.na(row$codon_N) &&
          !is.na(row$codon_C) &&
          nchar(row$codon_N) == 3 &&
          nchar(row$codon_C) == 3 &&
          row$codon_N != row$codon_C &&
          row$aa_N != row$aa_C

        tibble(
          pair_id = pair_id_value,
          tac102_pos = pos,
          aa_N = row$aa_N,
          aa_C = row$aa_C,
          codon_N = row$codon_N,
          codon_C = row$codon_C,
          context_NT_N = ctx_N,
          context_NT_C = ctx_C,
          is_valid_snp = valid_snp
        )
      }
    )
  )
}


# ==============================================================
# PRINT SNP RESULTS SAFELY
# ==============================================================

cat("Candidate SNP table:\n\n")

if (nrow(snp_results) == 0) {

  cat(
    "No candidate nonsynonymous nucleotide differences remain.\n"
  )

} else {

  print(
    snp_results %>%
      select(
        pair_id,
        tac102_pos,
        aa_N,
        aa_C,
        codon_N,
        codon_C,
        is_valid_snp
      )
  )
}

cat("\n")


# ==============================================================
# SAVE OUTPUTS
# ==============================================================

write_csv(
  frame_audit,
  "raw_read_analysis/reports/tac102_4B3h_direct_frame_audit.csv"
)

write_csv(
  audit_55_results,
  out_55_audit
)

write_csv(
  snp_results,
  out_4_snps
)


# ==============================================================
# SUMMARY TEXT
# ==============================================================

summary_text <- c(

  "Step 4B.3h: Direct Nucleotide-to-Protein Frame & Consistency Audit",

  "==============================================================",

  "",

  "INPUTS:",

  paste(
    "Internal N/C conflicts:",
    nrow(conflicts_df)
  ),

  paste(
    "Tail observations:",
    nrow(tails_df)
  ),

  paste(
    "04B.3g audit rows:",
    nrow(audit_04g_df)
  ),

  paste(
    "Usable tail observations:",
    nrow(tails_usable)
  ),

  "",

  "DIRECT FRAME AUDIT:",

  paste(
    "Tail sequences matching reported translation:",
    frame_1_matches
  ),

  paste(
    "Tail sequences requiring non-frame-1 forward translation:",
    non_frame_1
  ),

  paste(
    "Tail sequences requiring reverse complement:",
    reverse_required
  ),

  paste(
    "Tail sequences where reported translation cannot be reproduced:",
    not_reproduced
  ),

  "",

  "55-ANOMALY AUDIT:",

  paste(
    "Total identical-nucleotide / different-AA anomalies:",
    nrow(anomalies_55)
  ),

  "",

  "Diagnostic classification:",

  capture.output(
    print(anomaly_summary)
  ),

  "",

  "NONSYNONYMOUS CANDIDATE AUDIT:",

  paste(
    "Candidate nonsynonymous nucleotide differences:",
    nrow(snps_4)
  ),

  paste(
    "Validated nonsynonymous SNP candidates:",
    sum(
      snp_results$is_valid_snp,
      na.rm = TRUE
    )
  ),

  "",

  "INTERPRETATION:",

  "The 55 apparent identical-nucleotide / different-amino-acid conflicts do not constitute independent evidence of nucleotide variation.",

  "The direct nucleotide audit tests whether the reported translations can be reproduced directly from the underlying nucleotide sequences.",

  "Candidate nonsynonymous nucleotide differences are retained only when the underlying codons differ and the amino-acid calls differ.",

  "",

  "OUTPUTS:",

  "raw_read_analysis/reports/tac102_4B3h_direct_frame_audit.csv",

  out_55_audit,

  out_4_snps,

  out_summary
)


writeLines(
  summary_text,
  out_summary
)


# ==============================================================
# FINAL CONSOLE REPORT
# ==============================================================

cat("Saved reports:\n")

cat(
  " raw_read_analysis/reports/tac102_4B3h_direct_frame_audit.csv\n"
)

cat(
  " ",
  out_55_audit,
  "\n"
)

cat(
  " ",
  out_4_snps,
  "\n"
)

cat(
  " ",
  out_summary,
  "\n"
)

cat(
  "\n[COMPLETE] Step 4B.3h direct nucleotide/frame audit finished successfully.\n"
)