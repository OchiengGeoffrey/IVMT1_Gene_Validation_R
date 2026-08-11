# ==============================================================================
# Step 4B.3n: TAC102 Quality-Stratified Allele Persistence
#
# PURPOSE
# -------
# Evaluate persistence of the 653/698 alleles and phased haplotypes across
# progressively stricter joint Phred-quality thresholds.
#
# INPUT
# -----
# tac102_4B3m_haplotype_pairs.csv
#
# OUTPUT
# ------
# tac102_4B3n_quality_stratification.csv
# tac102_4B3n_quality_stratification_report.txt
#
# IMPORTANT
# ---------
# The quality filter is applied to min_phred_both, i.e. the minimum Phred
# quality observed across the two codons for each phased read pair.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------------------------

if (!exists("PATHS")) {
  source(
    here::here(
      "raw_read_analysis",
      "scripts",
      "00_setup.R"
    )
  )
}


cat(
  "==============================================================\n",
  " Step 4B.3n: TAC102 Quality-Stratified Allele Persistence\n",
  "==============================================================\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 1. FILE PATHS
# ------------------------------------------------------------------------------

INPUT_FILE <- file.path(
  PATHS$reports,
  "tac102_4B3m_haplotype_pairs.csv"
)

OUTPUT_CSV <- file.path(
  PATHS$reports,
  "tac102_4B3n_quality_stratification.csv"
)

OUTPUT_REPORT <- file.path(
  PATHS$reports,
  "tac102_4B3n_quality_stratification_report.txt"
)


if (!file.exists(INPUT_FILE)) {
  stop(
    "Required input file not found:\n",
    INPUT_FILE
  )
}


# ------------------------------------------------------------------------------
# 2. READ HAPLOTYPE DATA
# ------------------------------------------------------------------------------

haplo_df <- read.csv(
  INPUT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "pair_id",
  "aa_653",
  "aa_698",
  "haplotype_aa",
  "min_phred_both"
)


missing_columns <- setdiff(
  required_columns,
  names(haplo_df)
)


if (length(missing_columns) > 0L) {

  stop(
    "The following required columns are missing from the 4B.3m file:\n",
    paste(missing_columns, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 3. BASIC INPUT VALIDATION
# ------------------------------------------------------------------------------

haplo_df$min_phred_both <- as.numeric(
  haplo_df$min_phred_both
)

haplo_df$aa_653 <- as.character(
  haplo_df$aa_653
)

haplo_df$aa_698 <- as.character(
  haplo_df$aa_698
)

haplo_df$haplotype_aa <- as.character(
  haplo_df$haplotype_aa
)


if (any(is.na(haplo_df$min_phred_both))) {

  stop(
    "NA values detected in min_phred_both. ",
    "Quality values must be present for every phased record."
  )
}


cat(
  "Input file: ",
  INPUT_FILE,
  "\n",
  sep = ""
)

cat(
  "Phased records: ",
  nrow(haplo_df),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------------------------
# 4. QUALITY THRESHOLDS
# ------------------------------------------------------------------------------

thresholds <- c(
  0,
  10,
  15,
  20,
  25,
  30
)


# ------------------------------------------------------------------------------
# 5. QUALITY-STRATIFIED ANALYSIS
# ------------------------------------------------------------------------------

strat_list <- lapply(
  thresholds,
  function(q) {

    sub_df <- haplo_df[
      haplo_df$min_phred_both >= q,
      ,
      drop = FALSE
    ]


    n_pairs <- nrow(
      sub_df
    )


    # --------------------------------------------------------------------------
    # No records surviving threshold
    # --------------------------------------------------------------------------

    if (n_pairs == 0L) {

      return(
        data.frame(
          min_phred_threshold = q,
          shared_pairs_remaining = 0L,

          site_653_M = 0L,
          site_653_non_M = 0L,

          site_698_L = 0L,
          site_698_P = 0L,
          site_698_other = 0L,

          joint_M653_L698 = 0L,
          joint_alternative_haplotypes = 0L,

          pct_M653 = 0,
          pct_non_M653 = 0,

          pct_L698 = 0,
          pct_P698 = 0,

          pct_joint_M653_L698 = 0,
          pct_joint_alternative = 0,

          stringsAsFactors = FALSE
        )
      )
    }


    # --------------------------------------------------------------------------
    # Site 653
    # --------------------------------------------------------------------------

    m653 <- sum(
      sub_df$aa_653 == "M"
    )

    non_m653 <- sum(
      sub_df$aa_653 != "M"
    )


    # --------------------------------------------------------------------------
    # Site 698
    # --------------------------------------------------------------------------

    l698 <- sum(
      sub_df$aa_698 == "L"
    )

    p698 <- sum(
      sub_df$aa_698 == "P"
    )

    other_698 <- sum(
      !sub_df$aa_698 %in% c(
        "L",
        "P"
      )
    )


    # --------------------------------------------------------------------------
    # Joint haplotypes
    # --------------------------------------------------------------------------

    m653_l698 <- sum(
      sub_df$haplotype_aa == "M653-L698"
    )

    joint_alternative <- sum(
      sub_df$haplotype_aa != "M653-L698"
    )


    # --------------------------------------------------------------------------
    # Return one threshold row
    # --------------------------------------------------------------------------

    data.frame(

      min_phred_threshold = q,

      shared_pairs_remaining = n_pairs,

      site_653_M = m653,

      site_653_non_M = non_m653,

      site_698_L = l698,

      site_698_P = p698,

      site_698_other = other_698,

      joint_M653_L698 = m653_l698,

      joint_alternative_haplotypes = joint_alternative,

      pct_M653 = round(
        100 * m653 / n_pairs,
        1
      ),

      pct_non_M653 = round(
        100 * non_m653 / n_pairs,
        1
      ),

      pct_L698 = round(
        100 * l698 / n_pairs,
        1
      ),

      pct_P698 = round(
        100 * p698 / n_pairs,
        1
      ),

      pct_joint_M653_L698 = round(
        100 * m653_l698 / n_pairs,
        1
      ),

      pct_joint_alternative = round(
        100 * joint_alternative / n_pairs,
        1
      ),

      stringsAsFactors = FALSE
    )
  }
)


strat_matrix <- do.call(
  rbind,
  strat_list
)


# ------------------------------------------------------------------------------
# 6. PRINT RESULTS
# ------------------------------------------------------------------------------

cat(
  "\n==============================================================\n",
  " STEP 4B.3n: QUALITY STRATIFICATION MATRIX\n",
  "==============================================================\n\n",
  sep = ""
)


print(
  strat_matrix,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. SAVE CSV
# ------------------------------------------------------------------------------

write.csv(
  strat_matrix,
  OUTPUT_CSV,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. WRITE HUMAN-READABLE REPORT
# ------------------------------------------------------------------------------

report_lines <- c(

  "==============================================================",
  "TAC102 QUALITY-STRATIFIED ALLELE PERSISTENCE",
  "Step 4B.3n",
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
    "4B.3m haplotype file: ",
    INPUT_FILE
  ),
  "",
  paste0(
    "Input phased records: ",
    nrow(haplo_df)
  ),
  "",
  "METHOD",
  "Phased read pairs were progressively filtered using",
  "the minimum Phred quality observed across the two",
  "candidate codons (min_phred_both).",
  "",
  "The thresholds evaluated were:",
  paste(
    thresholds,
    collapse = ", "
  ),
  "",
  "RESULTS"
)


for (i in seq_len(nrow(strat_matrix))) {

  x <- strat_matrix[i, ]

  report_lines <- c(
    report_lines,
    "",
    paste0(
      "Q >= ",
      x$min_phred_threshold
    ),
    paste0(
      "Shared pairs remaining: ",
      x$shared_pairs_remaining
    ),
    paste0(
      "653 M: ",
      x$site_653_M,
      " (",
      x$pct_M653,
      "%)"
    ),
    paste0(
      "653 non-M: ",
      x$site_653_non_M,
      " (",
      x$pct_non_M653,
      "%)"
    ),
    paste0(
      "698 L: ",
      x$site_698_L,
      " (",
      x$pct_L698,
      "%)"
    ),
    paste0(
      "698 P: ",
      x$site_698_P,
      " (",
      x$pct_P698,
      "%)"
    ),
    paste0(
      "Joint M653-L698: ",
      x$joint_M653_L698,
      " (",
      x$pct_joint_M653_L698,
      "%)"
    ),
    paste0(
      "Joint alternative haplotypes: ",
      x$joint_alternative_haplotypes,
      " (",
      x$pct_joint_alternative,
      "%)"
    )
  )
}


report_lines <- c(
  report_lines,
  "",
  "INTERPRETATION",
  "",
  "This analysis evaluates persistence of candidate allele",
  "and haplotype calls as increasingly stringent quality",
  "thresholds are applied.",
  "",
  "The analysis does not by itself establish population",
  "fixation. Quality filtering is interpreted together with",
  "independent read-start families, strand/orientation",
  "consistency, sequencing depth, and independent genomic",
  "evidence.",
  "",
  "OUTPUT",
  paste0(
    "Quality-stratification CSV: ",
    OUTPUT_CSV
  ),
  paste0(
    "Report: ",
    OUTPUT_REPORT
  ),
  "",
  "=============================================================="
)


writeLines(
  report_lines,
  OUTPUT_REPORT
)


# ------------------------------------------------------------------------------
# 9. LOG
# ------------------------------------------------------------------------------

if (exists("LOG_FILE")) {

  log_info(
    paste0(
      "Step 4B.3n complete. ",
      nrow(haplo_df),
      " phased records evaluated across ",
      length(thresholds),
      " Phred thresholds."
    ),
    LOG_FILE
  )
}


cat(
  "\nStep 4B.3n complete.\n",
  "Quality-stratification CSV: ",
  OUTPUT_CSV,
  "\n",
  "Report: ",
  OUTPUT_REPORT,
  "\n",
  sep = ""
)