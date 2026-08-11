# ==============================================================================
# Step 4B.3m: TAC102 Dual-Position Haplotype & Linkage Analysis
#
# PURPOSE
# -------
# Analyze read pairs that contribute evidence at BOTH candidate positions
# (653 and 698), and phase the amino-acid/nucleotide calls at the individual
# paired-template level.
#
# IMPORTANT INTERPRETATION
# ------------------------
# This analysis establishes read-level haplotype structure among the candidate
# reads. It does NOT by itself establish population fixation.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. LOAD PROJECT ENVIRONMENT
# ------------------------------------------------------------------------------

if (!exists("PATHS") || !exists("LOG_FILE")) {
  source(
    here::here(
      "raw_read_analysis",
      "scripts",
      "00_setup.R"
    )
  )
}


cat("==============================================================\n")
cat(" Step 4B.3m: TAC102 Haplotype Linkage\n")
cat(" Positions 653 + 698\n")
cat("==============================================================\n\n")


# ------------------------------------------------------------------------------
# 1. FILE PATHS
# ------------------------------------------------------------------------------

DETAIL_FILE <- file.path(
  PATHS$reports,
  "tac102_4B3l_position_validation_detail.csv"
)

HAPLO_OUT_CSV <- file.path(
  PATHS$reports,
  "tac102_4B3m_haplotype_pairs.csv"
)

HAPLO_SUM_CSV <- file.path(
  PATHS$reports,
  "tac102_4B3m_haplotype_summary.csv"
)

HAPLO_REPORT <- file.path(
  PATHS$reports,
  "tac102_4B3m_haplotype_report.txt"
)


if (!file.exists(DETAIL_FILE)) {
  stop(
    "Input detail file not found: ",
    DETAIL_FILE
  )
}


# ------------------------------------------------------------------------------
# 2. READ 4B.3l DETAIL TABLE
# ------------------------------------------------------------------------------

detail <- read.csv(
  DETAIL_FILE,
  stringsAsFactors = FALSE
)


cat(
  "Input records: ",
  nrow(detail),
  "\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 3. VALIDATE REQUIRED COLUMNS
# ------------------------------------------------------------------------------

required_columns <- c(
  "pair_id",
  "mate",
  "orientation_role",
  "target_position",
  "hsp_strand",
  "codon",
  "codon_min_qual",
  "translated_aa",
  "dup_family"
)

missing_columns <- setdiff(
  required_columns,
  names(detail)
)

if (length(missing_columns) > 0L) {

  stop(
    "Required columns missing from 4B.3l detail file: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 4. ISOLATE VALID CALLS AT POSITIONS 653 AND 698
# ------------------------------------------------------------------------------

valid_detail <- detail[
  !is.na(detail$translated_aa) &
    nzchar(detail$translated_aa) &
    detail$translated_aa != "X" &
    detail$target_position %in% c(653L, 698L),
]


cat(
  "Valid 653/698 records: ",
  nrow(valid_detail),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 5. IDENTIFY SHARED READ PAIRS
#
# A shared pair_id must have a valid record at BOTH positions.
# ------------------------------------------------------------------------------

pairs_653 <- unique(
  valid_detail$pair_id[
    valid_detail$target_position == 653L
  ]
)

pairs_698 <- unique(
  valid_detail$pair_id[
    valid_detail$target_position == 698L
  ]
)

shared_pairs <- intersect(
  pairs_653,
  pairs_698
)


cat(
  "Position 653 read pairs: ",
  length(pairs_653),
  "\n",
  sep = ""
)

cat(
  "Position 698 read pairs: ",
  length(pairs_698),
  "\n",
  sep = ""
)

cat(
  "Shared read pairs covering BOTH positions: ",
  length(shared_pairs),
  "\n\n",
  sep = ""
)


if (length(shared_pairs) == 0L) {

  stop(
    "No shared read pairs were found between positions 653 and 698."
  )
}


# ------------------------------------------------------------------------------
# 6. FUNCTION TO SELECT THE BEST SITE RECORD
#
# If both mates provide a valid call for a position, select the call with the
# highest minimum Phred quality.
#
# The mate and orientation are retained in the output so provenance remains
# visible.
# ------------------------------------------------------------------------------

get_best_site_record <- function(
    pair_id,
    position
) {

  sub <- valid_detail[
    valid_detail$pair_id == pair_id &
      valid_detail$target_position == position,
  ]

  if (nrow(sub) == 0L) {
    return(NULL)
  }

  sub <- sub[
    order(
      -sub$codon_min_qual
    ),
    ,
    drop = FALSE
  ]

  sub[1, , drop = FALSE]
}


# ------------------------------------------------------------------------------
# 7. BUILD PHASED HAPLOTYPE RECORDS
# ------------------------------------------------------------------------------

haplo_list <- lapply(
  shared_pairs,
  function(pid) {

    rec_653 <- get_best_site_record(
      pid,
      653L
    )

    rec_698 <- get_best_site_record(
      pid,
      698L
    )


    if (
      is.null(rec_653) ||
      is.null(rec_698)
    ) {

      return(NULL)
    }


    data.frame(

      pair_id = pid,

      # ------------------------------------------------------------
      # Position 653
      # ------------------------------------------------------------

      codon_653 =
        rec_653$codon,

      aa_653 =
        rec_653$translated_aa,

      qmin_653 =
        rec_653$codon_min_qual,

      qmean_653 =
        rec_653$codon_mean_qual,

      mate_653 =
        rec_653$mate,

      orientation_653 =
        rec_653$orientation_role,

      strand_653 =
        rec_653$hsp_strand,

      family_653 =
        rec_653$dup_family,


      # ------------------------------------------------------------
      # Position 698
      # ------------------------------------------------------------

      codon_698 =
        rec_698$codon,

      aa_698 =
        rec_698$translated_aa,

      qmin_698 =
        rec_698$codon_min_qual,

      qmean_698 =
        rec_698$codon_mean_qual,

      mate_698 =
        rec_698$mate,

      orientation_698 =
        rec_698$orientation_role,

      strand_698 =
        rec_698$hsp_strand,

      family_698 =
        rec_698$dup_family,


      # ------------------------------------------------------------
      # Phased haplotype
      # ------------------------------------------------------------

      haplotype_aa =
        paste0(
          rec_653$translated_aa,
          "653-",
          rec_698$translated_aa,
          "698"
        ),

      haplotype_codon =
        paste0(
          rec_653$codon,
          "-",
          rec_698$codon
        ),


      # ------------------------------------------------------------
      # Combined quality
      # ------------------------------------------------------------

      min_phred_both =
        min(
          rec_653$codon_min_qual,
          rec_698$codon_min_qual
        ),

      mean_phred_both =
        mean(
          c(
            rec_653$codon_mean_qual,
            rec_698$codon_mean_qual
          )
        ),

      stringsAsFactors = FALSE

    )
  }
)


# Remove NULL records if any occurred
haplo_list <- haplo_list[
  !vapply(
    haplo_list,
    is.null,
    logical(1)
  )
]


if (length(haplo_list) == 0L) {

  stop(
    "No complete phased haplotype records could be constructed."
  )
}


haplo_df <- do.call(
  rbind,
  haplo_list
)


# ------------------------------------------------------------------------------
# 8. VERIFY EXPECTED NUMBER OF SHARED PAIRS
# ------------------------------------------------------------------------------

cat(
  "Complete phased records: ",
  nrow(haplo_df),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 9. HAPLOTYPE SUMMARY
# ------------------------------------------------------------------------------

haplo_summary <- aggregate(
  pair_id ~ haplotype_aa + haplotype_codon,
  data = haplo_df,
  FUN = length
)

names(haplo_summary)[
  names(haplo_summary) == "pair_id"
] <- "pair_count"


# ------------------------------------------------------------------------------
# 10. CALCULATE QUALITY SUMMARY SAFELY
# ------------------------------------------------------------------------------

quality_summary <- aggregate(
  min_phred_both ~ haplotype_aa + haplotype_codon,
  data = haplo_df,
  FUN = function(x) {
    round(
      mean(x),
      2
    )
  }
)

names(quality_summary)[
  names(quality_summary) == "min_phred_both"
] <- "mean_qmin_both"


quality_summary$min_qmin_both <- NA_real_
quality_summary$max_qmin_both <- NA_real_


for (i in seq_len(nrow(quality_summary))) {

  aa <- quality_summary$haplotype_aa[i]
  codon <- quality_summary$haplotype_codon[i]

  x <- haplo_df$min_phred_both[
    haplo_df$haplotype_aa == aa &
      haplo_df$haplotype_codon == codon
  ]

  quality_summary$min_qmin_both[i] <-
    min(x)

  quality_summary$max_qmin_both[i] <-
    max(x)
}


# ------------------------------------------------------------------------------
# 11. MERGE SUMMARY TABLES
# ------------------------------------------------------------------------------

haplo_summary <- merge(
  haplo_summary,
  quality_summary,
  by = c(
    "haplotype_aa",
    "haplotype_codon"
  ),
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------------------------
# 12. CALCULATE PERCENTAGE
# ------------------------------------------------------------------------------

haplo_summary$percent <- round(
  100 *
    haplo_summary$pair_count /
    sum(haplo_summary$pair_count),
  1
)


# ------------------------------------------------------------------------------
# 13. SORT BY ABUNDANCE
# ------------------------------------------------------------------------------

haplo_summary <- haplo_summary[
  order(
    -haplo_summary$pair_count,
    haplo_summary$haplotype_aa
  ),
]


# ------------------------------------------------------------------------------
# 14. SAVE OUTPUT TABLES
# ------------------------------------------------------------------------------

write.csv(
  haplo_df,
  HAPLO_OUT_CSV,
  row.names = FALSE
)

write.csv(
  haplo_summary,
  HAPLO_SUM_CSV,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 15. WRITE TEXT REPORT
# ------------------------------------------------------------------------------

report_lines <- c(

  "==============================================================",
  "TAC102 DUAL-POSITION HAPLOTYPE & LINKAGE ANALYSIS",
  "Step 4B.3m",
  "==============================================================",
  "",
  paste0(
    "Date: ",
    format(
      Sys.Date(),
      "%Y-%m-%d"
    )
  ),
  "",
  "INPUT",
  paste0(
    "4B.3l detail file: ",
    DETAIL_FILE
  ),
  "",
  "SHARED READ-PAIR ANALYSIS",
  paste0(
    "Position 653 read pairs: ",
    length(pairs_653)
  ),
  paste0(
    "Position 698 read pairs: ",
    length(pairs_698)
  ),
  paste0(
    "Shared read pairs: ",
    length(shared_pairs)
  ),
  paste0(
    "Complete phased records: ",
    nrow(haplo_df)
  ),
  "",
  "HAPLOTYPE SUMMARY",
  ""
)


for (i in seq_len(nrow(haplo_summary))) {

  report_lines <- c(
    report_lines,
    paste0(
      haplo_summary$haplotype_aa[i],
      " | ",
      haplo_summary$haplotype_codon[i],
      " | n=",
      haplo_summary$pair_count[i],
      " | ",
      haplo_summary$percent[i],
      "% | mean minimum Phred=",
      haplo_summary$mean_qmin_both[i]
    )
  )
}


report_lines <- c(
  report_lines,
  "",
  "INTERPRETATION",
  "",
  "This analysis phases the two candidate positions among shared",
  "read pairs. It evaluates whether alternative amino-acid calls",
  "at positions 653 and 698 occur on the same paired sequencing",
  "templates.",
  "",
  "The analysis does NOT by itself establish population fixation",
  "of any allele. Candidate-read enrichment, sequencing depth,",
  "duplicate structure, and independent genomic confirmation must",
  "be considered before biological fixation is inferred.",
  "",
  "OUTPUTS",
  paste0(
    "Haplotype records: ",
    HAPLO_OUT_CSV
  ),
  paste0(
    "Haplotype summary: ",
    HAPLO_SUM_CSV
  ),
  paste0(
    "Text report: ",
    HAPLO_REPORT
  ),
  "",
  "=============================================================="
)


writeLines(
  report_lines,
  HAPLO_REPORT
)


# ------------------------------------------------------------------------------
# 16. PRINT RESULTS
# ------------------------------------------------------------------------------

cat(
  "==============================================================\n"
)

cat(
  " PHASED HAPLOTYPE SUMMARY\n"
)

cat(
  "==============================================================\n"
)

print(
  haplo_summary,
  row.names = FALSE
)


cat(
  "\n==============================================================\n"
)

cat(
  " DETAILED PHASED RECORDS\n"
)

cat(
  "==============================================================\n"
)

print(
  haplo_df[
    ,
    c(
      "pair_id",
      "haplotype_aa",
      "haplotype_codon",
      "mate_653",
      "mate_698",
      "qmin_653",
      "qmin_698",
      "min_phred_both"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 17. FINAL LOGGING
# ------------------------------------------------------------------------------

log_info(
  paste0(
    "Step 4B.3m complete. ",
    "Phased ",
    nrow(haplo_df),
    " shared read-pair records. ",
    "Haplotype summaries written to reports."
  ),
  LOG_FILE
)


cat(
  "\nStep 4B.3m complete.\n"
)

cat(
  "Haplotype records: ",
  HAPLO_OUT_CSV,
  "\n",
  sep = ""
)

cat(
  "Haplotype summary: ",
  HAPLO_SUM_CSV,
  "\n",
  sep = ""
)

cat(
  "Report: ",
  HAPLO_REPORT,
  "\n",
  sep = ""
)