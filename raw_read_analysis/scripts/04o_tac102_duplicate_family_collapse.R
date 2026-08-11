# ==============================================================================
# Step 4B.3o: TAC102 Duplicate Read-Start Family Collapse
#
# PURPOSE
# -------
# Evaluate whether the observed 653/698 haplotypes are supported by independent
# read-start families rather than being inflated by repeated sequencing reads
# originating from the same apparent template.
#
# INPUT
# -----
# tac102_4B3m_haplotype_pairs.csv
#
# OUTPUT
# ------
# tac102_4B3o_family_collapse.csv
# tac102_4B3o_independent_family_summary.csv
# tac102_4B3o_family_collapse_report.txt
#
# IMPORTANT
# ---------
# The primary analysis uses Q >= 10, matching the proposed Step 4B.3o design.
#
# A "family combination" is defined by:
#
#     family_653 + family_698
#
# This represents the read-start family observed at each of the two candidate
# positions on the phased read pair.
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
  " Step 4B.3o: TAC102 Duplicate Family Collapse\n",
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

OUTPUT_COLLAPSE_CSV <- file.path(
  PATHS$reports,
  "tac102_4B3o_family_collapse.csv"
)

OUTPUT_FAMILY_CSV <- file.path(
  PATHS$reports,
  "tac102_4B3o_independent_family_summary.csv"
)

OUTPUT_REPORT <- file.path(
  PATHS$reports,
  "tac102_4B3o_family_collapse_report.txt"
)


if (!file.exists(INPUT_FILE)) {

  stop(
    "Required input file not found:\n",
    INPUT_FILE
  )
}


# ------------------------------------------------------------------------------
# 2. READ INPUT
# ------------------------------------------------------------------------------

haplo_df <- read.csv(
  INPUT_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "pair_id",
  "family_653",
  "family_698",
  "haplotype_aa",
  "aa_653",
  "aa_698",
  "min_phred_both"
)


missing_columns <- setdiff(
  required_columns,
  names(haplo_df)
)


if (length(missing_columns) > 0L) {

  stop(
    "The following required columns are missing:\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 3. STANDARDIZE COLUMNS
# ------------------------------------------------------------------------------

haplo_df$min_phred_both <- as.numeric(
  haplo_df$min_phred_both
)

haplo_df$family_653 <- as.character(
  haplo_df$family_653
)

haplo_df$family_698 <- as.character(
  haplo_df$family_698
)

haplo_df$haplotype_aa <- as.character(
  haplo_df$haplotype_aa
)

haplo_df$aa_653 <- as.character(
  haplo_df$aa_653
)

haplo_df$aa_698 <- as.character(
  haplo_df$aa_698
)


# ------------------------------------------------------------------------------
# 4. QUALITY THRESHOLD
# ------------------------------------------------------------------------------

QUALITY_THRESHOLD <- 10


q10_df <- haplo_df[
  haplo_df$min_phred_both >= QUALITY_THRESHOLD,
  ,
  drop = FALSE
]


cat(
  "Input phased records: ",
  nrow(haplo_df),
  "\n",
  sep = ""
)

cat(
  "Records at Q >= ",
  QUALITY_THRESHOLD,
  ": ",
  nrow(q10_df),
  "\n\n",
  sep = ""
)


if (nrow(q10_df) == 0L) {

  stop(
    "No phased records survive Q >= ",
    QUALITY_THRESHOLD,
    "."
  )
}


# ------------------------------------------------------------------------------
# 5. COLLAPSE BY FAMILY COMBINATION
#
# Each unique combination of family_653 + family_698 represents one observed
# paired read-start-family combination.
# ------------------------------------------------------------------------------

family_collapse <- aggregate(
  pair_id ~
    family_653 +
    family_698 +
    haplotype_aa,
  data = q10_df,
  FUN = function(x) length(unique(x))
)


names(family_collapse)[
  names(family_collapse) == "pair_id"
] <- "pair_count"


# ------------------------------------------------------------------------------
# 6. SORT FAMILY COLLAPSE TABLE
# ------------------------------------------------------------------------------

family_collapse <- family_collapse[
  order(
    -family_collapse$pair_count,
    family_collapse$haplotype_aa,
    family_collapse$family_653,
    family_collapse$family_698
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# 7. ADD PERCENTAGE
# ------------------------------------------------------------------------------

total_family_observations <- sum(
  family_collapse$pair_count
)


family_collapse$percent <- round(
  100 *
    family_collapse$pair_count /
    total_family_observations,
  1
)


# ------------------------------------------------------------------------------
# 8. IDENTIFY UNIQUE FAMILIES SUPPORTING EACH SITE
# ------------------------------------------------------------------------------

families_653 <- unique(
  q10_df$family_653
)

families_698 <- unique(
  q10_df$family_698
)


# ------------------------------------------------------------------------------
# 9. IDENTIFY UNIQUE FAMILY COMBINATIONS SUPPORTING M653-L698
# ------------------------------------------------------------------------------

m653_l698_df <- q10_df[
  q10_df$haplotype_aa == "M653-L698",
  ,
  drop = FALSE
]


m653_l698_family_combinations <- unique(
  m653_l698_df[
    ,
    c(
      "family_653",
      "family_698"
    ),
    drop = FALSE
  ]
)


# ------------------------------------------------------------------------------
# 10. SITE-LEVEL FAMILY SUPPORT
# ------------------------------------------------------------------------------

m653_families <- unique(
  q10_df$family_653[
    q10_df$aa_653 == "M"
  ]
)


non_m653_families <- unique(
  q10_df$family_653[
    q10_df$aa_653 != "M"
  ]
)


l698_families <- unique(
  q10_df$family_698[
    q10_df$aa_698 == "L"
  ]
)


p698_families <- unique(
  q10_df$family_698[
    q10_df$aa_698 == "P"
  ]
)


# ------------------------------------------------------------------------------
# 11. INDEPENDENT FAMILY SUMMARY TABLE
# ------------------------------------------------------------------------------

family_summary <- data.frame(

  metric = c(

    "Q threshold",

    "Q-qualified phased pairs",

    "Unique 653 read-start families",

    "Unique 698 read-start families",

    "653 M read-start families",

    "653 non-M read-start families",

    "698 L read-start families",

    "698 P read-start families",

    "M653-L698 unique family combinations",

    "M653-L698 phased pairs",

    "Alternative haplotype phased pairs"

  ),

  value = c(

    as.character(
      QUALITY_THRESHOLD
    ),

    as.character(
      nrow(q10_df)
    ),

    as.character(
      length(families_653)
    ),

    as.character(
      length(families_698)
    ),

    as.character(
      length(m653_families)
    ),

    as.character(
      length(non_m653_families)
    ),

    as.character(
      length(l698_families)
    ),

    as.character(
      length(p698_families)
    ),

    as.character(
      nrow(m653_l698_family_combinations)
    ),

    as.character(
      nrow(m653_l698_df)
    ),

    as.character(
      sum(
        q10_df$haplotype_aa != "M653-L698"
      )
    )

  ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 12. PRINT FAMILY COLLAPSE
# ------------------------------------------------------------------------------

cat(
  "\n==============================================================\n",
  " STEP 4B.3o: DUPLICATE FAMILY COLLAPSE\n",
  " Q >= ",
  QUALITY_THRESHOLD,
  "\n",
  "==============================================================\n\n",
  sep = ""
)


print(
  family_collapse,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. PRINT INDEPENDENT FAMILY SUMMARY
# ------------------------------------------------------------------------------

cat(
  "\n==============================================================\n",
  " INDEPENDENT READ-START FAMILY SUMMARY\n",
  "==============================================================\n\n",
  sep = ""
)


print(
  family_summary,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 14. PRINT M653-L698 FAMILY COMBINATIONS
# ------------------------------------------------------------------------------

cat(
  "\n==============================================================\n",
  " M653-L698 INDEPENDENT FAMILY COMBINATIONS\n",
  "==============================================================\n\n",
  sep = ""
)


if (
  nrow(
    m653_l698_family_combinations
  ) > 0L
) {

  print(
    m653_l698_family_combinations,
    row.names = FALSE
  )

} else {

  cat(
    "No M653-L698 family combinations detected at Q >= ",
    QUALITY_THRESHOLD,
    ".\n",
    sep = ""
  )
}


# ------------------------------------------------------------------------------
# 15. SAVE OUTPUT TABLES
# ------------------------------------------------------------------------------

write.csv(
  family_collapse,
  OUTPUT_COLLAPSE_CSV,
  row.names = FALSE
)


write.csv(
  family_summary,
  OUTPUT_FAMILY_CSV,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 16. HUMAN-READABLE REPORT
# ------------------------------------------------------------------------------

report_lines <- c(

  "==============================================================",
  "TAC102 DUPLICATE READ-START FAMILY COLLAPSE",
  "Step 4B.3o",
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
  "METHOD",
  paste0(
    "Joint Phred threshold: Q >= ",
    QUALITY_THRESHOLD
  ),
  "",
  "Read pairs were grouped by the combination of the",
  "independent read-start family at position 653 and",
  "the independent read-start family at position 698.",
  "",
  "RESULTS",
  paste0(
    "Input phased records: ",
    nrow(haplo_df)
  ),
  paste0(
    "Q-qualified phased records: ",
    nrow(q10_df)
  ),
  paste0(
    "Unique 653 read-start families: ",
    length(families_653)
  ),
  paste0(
    "Unique 698 read-start families: ",
    length(families_698)
  ),
  paste0(
    "653 M-supporting families: ",
    length(m653_families)
  ),
  paste0(
    "653 non-M families: ",
    length(non_m653_families)
  ),
  paste0(
    "698 L-supporting families: ",
    length(l698_families)
  ),
  paste0(
    "698 P-supporting families: ",
    length(p698_families)
  ),
  paste0(
    "M653-L698 phased pairs: ",
    nrow(m653_l698_df)
  ),
  paste0(
    "M653-L698 unique family combinations: ",
    nrow(m653_l698_family_combinations)
  ),
  "",
  "INTERPRETATION",
  "",
  "Family collapse is used to assess whether observations",
  "are concentrated within repeated read-start families.",
  "",
  "A haplotype supported by multiple distinct family",
  "combinations provides stronger evidence of independent",
  "molecular representation than repeated observations",
  "within one family.",
  "",
  "This analysis does not independently establish allele",
  "fixation or population frequency.",
  "",
  "OUTPUTS",
  paste0(
    "Family-collapse table: ",
    OUTPUT_COLLAPSE_CSV
  ),
  paste0(
    "Independent-family summary: ",
    OUTPUT_FAMILY_CSV
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
# 17. LOG
# ------------------------------------------------------------------------------

if (exists("LOG_FILE")) {

  log_info(
    paste0(
      "Step 4B.3o complete. ",
      nrow(q10_df),
      " Q-qualified phased records collapsed by read-start family."
    ),
    LOG_FILE
  )
}


cat(
  "\nStep 4B.3o complete.\n",
  "Family-collapse CSV: ",
  OUTPUT_COLLAPSE_CSV,
  "\n",
  "Independent-family summary: ",
  OUTPUT_FAMILY_CSV,
  "\n",
  "Report: ",
  OUTPUT_REPORT,
  "\n",
  sep = ""
)