# ==============================================================================
# Step 4B.3l: TAC102 Direct FASTQ Base-Level Validation
#
# PROJECT
# -------
# Genomic & Transcriptomic Limits on kDNA in T. equiperdum
#
# PURPOSE
# -------
# Independently validate the two recurrent TAC102 candidate substitutions
# identified in Step 4B.3k:
#
#     TAC102 position 653
#     TAC102 position 698
#
# The nucleotide sequence and Phred quality are re-extracted directly from
# the ORIGINAL R1/R2 FASTQ files.
#
# IMPORTANT
# ---------
# This script does NOT use the following values to make the nucleotide call:
#
#     gap_facing_sequence
#     reconstructed_codon
#     reconstructed_aa
#     translated_candidate
#     4B.3j reconstructed nucleotide values
#     4B.3k provenance nucleotide values
#
# Previous analysis is used ONLY to identify candidate reads that should be
# examined.
#
# CORRECTION
# ----------
# TAC102 positions 653 and 698 are INTERNAL protein positions.
#
# Biostrings::translate() defaults to:
#
#     no.init.codon = FALSE
#
# Under that setting, CTG and TTG can be treated as alternative initiation
# codons and returned as methionine.
#
# That behaviour is inappropriate for internal TAC102 codons.
#
# Therefore all codons in this script are translated with:
#
#     no.init.codon = TRUE
#
# Consequently:
#
#     ATG -> M
#     CTG -> L
#     TTG -> L
#     GTG -> V
#
# under the standard genetic code.
#
# ADDITIONAL VALIDATION
# ---------------------
# 1. Original FASTQ sequence extraction
# 2. Original FASTQ Phred extraction
# 3. Independent coordinate reconstruction
# 4. Strand/orientation handling
# 5. Reverse-complement quality handling
# 6. Codon-level quality
# 7. Independent read-start families
# 8. Cross-position shared-read-pair assessment
# 9. Explicit translation sanity check
# 10. Complete output-column consistency
#
# ==============================================================================


# ==============================================================================
# 0. INITIALIZATION
# ==============================================================================

cat(
  "==============================================================\n",
  " Step 4B.3l: TAC102 Direct FASTQ Base-Level Validation\n",
  " Positions: 653, 698\n",
  "==============================================================\n\n",
  sep = ""
)


if (
  !exists("PATHS") ||
  !exists("LOG_FILE") ||
  !exists("CONFIG")
) {

  source(
    here::here(
      "raw_read_analysis",
      "scripts",
      "00_setup.R"
    )
  )
}


log_info(
  "Starting Step 4B.3l: direct FASTQ base-level TAC102 validation.",
  LOG_FILE
)


# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

TARGET_POSITIONS <- c(
  653L,
  698L
)

DUP_WINDOW <- 3L

R1_FASTQ <- PATHS$fastq_r1
R2_FASTQ <- PATHS$fastq_r2

TAILS_FILE <- file.path(
  PATHS$reports,
  "tac102_41pairs_unaligned_tails.csv"
)

DETAIL_OUT <- file.path(
  PATHS$reports,
  "tac102_4B3l_position_validation_detail.csv"
)

SUMMARY_OUT <- file.path(
  PATHS$reports,
  "tac102_4B3l_position_validation_summary.csv"
)

REPORT_OUT <- file.path(
  PATHS$reports,
  "tac102_4B3l_position_validation_report.txt"
)


# ==============================================================================
# 2. VERIFY INPUT FILES
# ==============================================================================

required_files <- c(
  R1_FASTQ,
  R2_FASTQ,
  TAILS_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {

  stop(
    paste0(
      "Required input file(s) not found:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


cat(
  "Input files verified:\n",
  "R1 FASTQ: ",
  R1_FASTQ,
  "\n",
  "R2 FASTQ: ",
  R2_FASTQ,
  "\n",
  "Candidate table: ",
  TAILS_FILE,
  "\n\n",
  sep = ""
)


log_info(
  paste0(
    "Input files verified. ",
    "R1=", R1_FASTQ,
    "; R2=", R2_FASTQ,
    "; candidates=", TAILS_FILE
  ),
  LOG_FILE
)


# ==============================================================================
# 3. LOAD CANDIDATE READ TABLE
#
# The candidate table supplies:
#
#   pair identity
#   mate
#   orientation
#   original BLAST HSP coordinates
#   gap-position coverage
#
# It does NOT supply the nucleotide call used for final validation.
# ==============================================================================

tails <- read.csv(
  TAILS_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "pair_id",
  "mate",
  "orientation_role",
  "hsp_qstart",
  "hsp_sstart",
  "hsp_send",
  "gap_position_start",
  "gap_position_end"
)


missing_columns <- setdiff(
  required_columns,
  names(tails)
)


if (length(missing_columns) > 0L) {

  stop(
    paste0(
      "Missing required columns in candidate table:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}


# Convert required coordinate columns explicitly to numeric/integer.
numeric_columns <- c(
  "hsp_qstart",
  "hsp_sstart",
  "hsp_send",
  "gap_position_start",
  "gap_position_end"
)


for (nm in numeric_columns) {

  tails[[nm]] <- suppressWarnings(
    as.numeric(
      tails[[nm]]
    )
  )
}


cat(
  "Candidate HSP/read entries loaded: ",
  nrow(tails),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 4. IDENTIFY READ/POSITION COMBINATIONS COVERING 653 OR 698
#
# A read can cover both positions. Therefore a separate candidate record is
# retained for each target position.
# ==============================================================================

covers_position <- function(
    row,
    position) {

  if (
    is.na(row$gap_position_start) ||
    is.na(row$gap_position_end)
  ) {

    return(FALSE)
  }


  position >= row$gap_position_start &&
    position <= row$gap_position_end
}


candidate_rows <- list()


for (i in seq_len(nrow(tails))) {

  row <- tails[i, ]


  for (position in TARGET_POSITIONS) {

    if (
      covers_position(
        row,
        position
      )
    ) {

      candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(

        pair_id =
          as.character(
            row$pair_id
          ),

        mate =
          as.character(
            row$mate
          ),

        orientation_role =
          as.character(
            row$orientation_role
          ),

        target_position =
          as.integer(
            position
          ),

        hsp_qstart_raw =
          as.integer(
            row$hsp_qstart
          ),

        hsp_sstart_raw =
          as.integer(
            row$hsp_sstart
          ),

        hsp_send_raw =
          as.integer(
            row$hsp_send
          ),

        stringsAsFactors = FALSE
      )
    }
  }
}


if (
  length(candidate_rows) == 0L
) {

  stop(
    "No candidate reads cover TAC102 positions 653 or 698."
  )
}


candidates <- do.call(
  rbind,
  candidate_rows
)


row.names(candidates) <- NULL


n_653 <- sum(
  candidates$target_position == 653L
)

n_698 <- sum(
  candidates$target_position == 698L
)


cat(
  "Candidate read/position combinations:\n",
  "  Position 653: ",
  n_653,
  "\n",
  "  Position 698: ",
  n_698,
  "\n\n",
  sep = ""
)


log_info(
  paste0(
    "Identified ",
    nrow(candidates),
    " candidate-read/position combinations. ",
    "653=",
    n_653,
    "; 698=",
    n_698
  ),
  LOG_FILE
)


target_read_ids <- unique(
  as.character(
    candidates$pair_id
  )
)


cat(
  "Unique candidate read-pair IDs: ",
  length(target_read_ids),
  "\n\n",
  sep = ""
)


# ==============================================================================
# 5. FASTQ READ IDENTIFIER NORMALIZATION
# ==============================================================================

clean_read_id <- function(x) {

  x <- as.character(x)


  # Remove everything after first whitespace.
  x <- sub(
    "\\s.*$",
    "",
    x
  )


  # Remove /1 or /2 suffix.
  x <- sub(
    "/[12]$",
    "",
    x
  )


  x
}


# ==============================================================================
# 6. STANDARD OUTPUT ROW CONSTRUCTOR
#
# This function is critical.
#
# Every extraction path, including error paths, returns exactly the same
# columns and compatible data types. This prevents the rbind/data.frame error
# encountered in the previous run.
# ==============================================================================

empty_result_row <- function(
    pair_id,
    mate,
    orientation_role,
    target_position,
    hsp_strand = NA_character_,
    read_start_coord = NA_integer_,
    read_length = NA_integer_,
    codon_nt_start = NA_integer_,
    codon_nt_end = NA_integer_,
    codon = NA_character_,
    codon_base1 = NA_character_,
    codon_base2 = NA_character_,
    codon_base3 = NA_character_,
    codon_qual_pos1 = NA_integer_,
    codon_qual_pos2 = NA_integer_,
    codon_qual_pos3 = NA_integer_,
    codon_min_qual = NA_integer_,
    codon_mean_qual = NA_real_,
    translated_aa = NA_character_,
    note = "",
    dup_family = NA_character_) {

  data.frame(

    pair_id =
      as.character(pair_id),

    mate =
      as.character(mate),

    orientation_role =
      as.character(orientation_role),

    target_position =
      as.integer(target_position),

    hsp_strand =
      as.character(hsp_strand),

    read_start_coord =
      as.integer(read_start_coord),

    read_length =
      as.integer(read_length),

    codon_nt_start =
      as.integer(codon_nt_start),

    codon_nt_end =
      as.integer(codon_nt_end),

    codon =
      as.character(codon),

    codon_base1 =
      as.character(codon_base1),

    codon_base2 =
      as.character(codon_base2),

    codon_base3 =
      as.character(codon_base3),

    codon_qual_pos1 =
      as.integer(codon_qual_pos1),

    codon_qual_pos2 =
      as.integer(codon_qual_pos2),

    codon_qual_pos3 =
      as.integer(codon_qual_pos3),

    codon_min_qual =
      as.integer(codon_min_qual),

    codon_mean_qual =
      as.numeric(codon_mean_qual),

    translated_aa =
      as.character(translated_aa),

    note =
      as.character(note),

    dup_family =
      as.character(dup_family),

    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 7. STREAM ORIGINAL FASTQ FILE
#
# IMPORTANT
# ---------
# Do NOT use:
#
#     ShortRead::quality()
#
# because quality() is not exported from the ShortRead namespace in the
# installed version.
#
# The correct accessor is:
#
#     quality(reads)
# ==============================================================================

stream_candidate_fastq <- function(
    fastq_file,
    target_ids,
    mate_label) {

  target_ids <- unique(
    as.character(target_ids)
  )


  sequences <- vector(
    "list",
    length(target_ids)
  )

  qualities <- vector("list", length(target_ids))


  names(sequences) <- target_ids
  names(qualities) <- target_ids


  found <- character(0)


  cat(
    "Streaming ",
    mate_label,
    " FASTQ: ",
    fastq_file,
    "\n",
    sep = ""
  )


  fq <- ShortRead::FastqStreamer(
    fastq_file,
    n = 100000L
  )


  on.exit(
    try(
      close(fq),
      silent = TRUE
    ),
    add = TRUE
  )


  chunk_number <- 0L


  repeat {

    chunk_number <- chunk_number + 1L


    reads <- ShortRead::yield(
      fq
    )


    if (
      length(reads) == 0L
    ) {

      break
    }


    read_names <- clean_read_id(
      ShortRead::id(
        reads
      )
    )


    keep <- read_names %in% target_ids


    if (
      any(keep)
    ) {

      selected <- reads[keep]

      selected_names <- read_names[keep]


      # ------------------------------------------------------------------------
      # Sequence
      # ------------------------------------------------------------------------

      selected_sequences <- ShortRead::sread(selected)


      # ------------------------------------------------------------
      # Extract Phred quality scores directly from the ShortReadQ
      # object.
      #
      # Do NOT use:
      #   ShortRead::quality()
      #   quality()
      #
      # ShortRead::quality is not an exported function in the
      # installed ShortRead version. Direct slot access is explicit
      # and avoids namespace/method-dispatch ambiguity.
      # ------------------------------------------------------------

      selected_quality <- slot(
        selected,
        "quality"
      )

      quality_matrix <- as(
        selected_quality,
        "matrix"
      )

      selected_quality_int <- lapply(
        seq_len(
          nrow(
            quality_matrix
          )
        ),
        function(i) {

          as.integer(
            quality_matrix[i, ]
          )
        }
      )

      names(selected_quality_int) <- selected_names


      # ------------------------------------------------------------------------
      # Store each recovered read
      # ------------------------------------------------------------------------

      for (
        k in seq_along(selected_names)
      ) {

        read_id <- selected_names[k]


        sequences[[read_id]] <-
          selected_sequences[k]


        qualities[[read_id]] <-
          selected_quality_int[[k]]


        found <- c(
          found,
          read_id
        )
      }
    }


    if (
      chunk_number %% 10L == 0L
    ) {

      cat(
        "  Processed FASTQ chunks: ",
        chunk_number,
        "\n",
        sep = ""
      )
    }
  }


  found <- unique(
    found
  )


  missing <- setdiff(
    target_ids,
    found
  )


  if (
    length(missing) > 0L
  ) {

    log_info(
      paste0(
        mate_label,
        ": ",
        length(missing),
        " target reads not found in FASTQ."
      ),
      LOG_FILE
    )

    cat(
      "  WARNING: ",
      length(missing),
      " target reads not found.\n",
      sep = ""
    )
  }


  list(
    sequences = sequences,
    qualities = qualities,
    found = found
  )
}


# ==============================================================================
# 8. RE-STREAM ORIGINAL FASTQ FILES
# ==============================================================================

cat(
  "Re-streaming ORIGINAL FASTQ files...\n\n"
)


r1_data <- stream_candidate_fastq(
  R1_FASTQ,
  target_read_ids,
  "R1"
)


r2_data <- stream_candidate_fastq(
  R2_FASTQ,
  target_read_ids,
  "R2"
)


cat(
  "\nFASTQ recovery summary:\n",
  "  R1: ",
  length(r1_data$found),
  "/",
  length(target_read_ids),
  "\n",
  "  R2: ",
  length(r2_data$found),
  "/",
  length(target_read_ids),
  "\n\n",
  sep = ""
)


log_info(
  paste0(
    "Original FASTQ re-extraction complete. ",
    "R1=",
    length(r1_data$found),
    "/",
    length(target_read_ids),
    "; R2=",
    length(r2_data$found),
    "/",
    length(target_read_ids)
  ),
  LOG_FILE
)


# ==============================================================================
# 9. CODON WINDOW CALCULATION
#
# Derived ONLY from:
#
#     hsp_qstart_raw
#     hsp_sstart_raw
#     hsp_send_raw
#     TAC102 protein position
#
# Coordinates are 1-based FASTQ coordinates.
# ==============================================================================

get_codon_window <- function(
    qstart_raw,
    sstart_raw,
    send_raw,
    protein_position) {

  qstart_raw <- as.integer(
    qstart_raw
  )

  sstart_raw <- as.integer(
    sstart_raw
  )

  send_raw <- as.integer(
    send_raw
  )

  protein_position <- as.integer(
    protein_position
  )


  if (
    any(
      is.na(
        c(
          qstart_raw,
          sstart_raw,
          send_raw,
          protein_position
        )
      )
    )
  ) {

    return(
      list(
        nt_start = NA_integer_,
        nt_end = NA_integer_,
        forward = NA
      )
    )
  }


  forward <- (
    sstart_raw < send_raw
  )


  if (
    forward
  ) {

    nt_start <- (
      sstart_raw +
        (protein_position - qstart_raw) * 3L
    )

    nt_end <- (
      nt_start + 2L
    )

  } else {

    nt_high <- (
      sstart_raw -
        (protein_position - qstart_raw) * 3L
    )

    nt_start <- (
      nt_high - 2L
    )

    nt_end <- nt_high
  }


  list(
    nt_start =
      as.integer(
        nt_start
      ),

    nt_end =
      as.integer(
        nt_end
      ),

    forward =
      forward
  )
}


# ==============================================================================
# 10. INTERNAL-CDS CODON TRANSLATION
#
# CRITICAL:
#
# TAC102 positions 653 and 698 are internal protein positions.
#
# Therefore no.init.codon = TRUE is mandatory.
# ==============================================================================

translate_codon <- function(
    codon) {

  codon <- toupper(
    as.character(codon)
  )


  if (
    length(codon) != 1L ||
    nchar(codon) != 3L ||
    grepl(
      "[^ACGT]",
      codon
    )
  ) {

    return("X")
  }


  tryCatch(

    as.character(
      Biostrings::translate(
        Biostrings::DNAString(codon),
        no.init.codon = TRUE,
        if.fuzzy.codon = "X"
      )
    ),

    error = function(e) {

      "X"
    }
  )
}


# ============================================================
# TRANSLATION SANITY CHECK
# Internal CDS translation: alternative start-codon behavior
# must NOT be applied.
# ============================================================

translation_test_codons <- c(
  "ATG",
  "CTG",
  "TTG",
  "GTG",
  "AAA",
  "TTT",
  "GGG"
)

translation_test_expected <- c(
  "M",
  "L",
  "L",
  "V",
  "K",
  "F",
  "G"
)

translation_test_observed <- vapply(
  translation_test_codons,
  function(codon) {

    as.character(
      Biostrings::translate(
        Biostrings::DNAString(codon),
        no.init.codon = TRUE,
        if.fuzzy.codon = "X"
      )
    )

  },
  character(1)
)

names(translation_test_observed) <- translation_test_codons

if (
  !identical(
    unname(translation_test_observed),
    translation_test_expected
  )
) {

  stop(
    paste0(
      "Translation sanity check FAILED.\n",
      "Expected internal-CDS standard-code translation.\n",
      "Observed: ",
      paste(
        paste0(
          names(translation_test_observed),
          "=",
          translation_test_observed
        ),
        collapse = ", "
      )
    )
  )
}

cat(
  "Translation sanity check PASSED.\n",
  "Internal CDS standard-code translation confirmed:\n",
  paste(
    paste0(
      names(translation_test_observed),
      "=",
      translation_test_observed
    ),
    collapse = ", "
  ),
  "\n"
)


# ==============================================================================
# 12. PER-READ, PER-POSITION BASE-LEVEL EXTRACTION
# ==============================================================================

extract_base_evidence <- function(
    row) {

  pair_id <- as.character(
    row$pair_id
  )

  mate <- as.character(
    row$mate
  )

  orientation_role <- as.character(
    row$orientation_role
  )

  position <- as.integer(
    row$target_position
  )


  # --------------------------------------------------------------------------
  # Retrieve sequence + quality from original FASTQ stream
  # --------------------------------------------------------------------------

  if (
    mate == "R1"
  ) {

    seq_obj <- r1_data$sequences[[pair_id]]

    qual_obj <- r1_data$qualities[[pair_id]]

  } else if (
    mate == "R2"
  ) {

    seq_obj <- r2_data$sequences[[pair_id]]

    qual_obj <- r2_data$qualities[[pair_id]]

  } else {

    return(
      empty_result_row(
        pair_id =
          pair_id,

        mate =
          mate,

        orientation_role =
          orientation_role,

        target_position =
          position,

        note =
          "unknown mate"
      )
    )
  }


  # --------------------------------------------------------------------------
  # Check whether read was recovered
  # --------------------------------------------------------------------------

  if (
    is.null(seq_obj) ||
    is.null(qual_obj)
  ) {

    return(
      empty_result_row(
        pair_id =
          pair_id,

        mate =
          mate,

        orientation_role =
          orientation_role,

        target_position =
          position,

        note =
          "read or quality not found in original FASTQ"
      )
    )
  }


  # --------------------------------------------------------------------------
  # Read length
  # --------------------------------------------------------------------------

  read_length <- length(
    seq_obj
  )


  # --------------------------------------------------------------------------
  # Re-derive codon coordinates
  # --------------------------------------------------------------------------

  window <- get_codon_window(

    qstart_raw =
      row$hsp_qstart_raw,

    sstart_raw =
      row$hsp_sstart_raw,

    send_raw =
      row$hsp_send_raw,

    protein_position =
      position
  )


  if (
    is.na(window$nt_start) ||
    is.na(window$nt_end)
  ) {

    return(
      empty_result_row(
        pair_id =
          pair_id,

        mate =
          mate,

        orientation_role =
          orientation_role,

        target_position =
          position,

        hsp_strand =
          NA_character_,

        read_length =
          read_length,

        note =
          "invalid HSP coordinate(s)"
      )
    )
  }


  # --------------------------------------------------------------------------
  # Verify codon window is inside read
  # --------------------------------------------------------------------------

  if (
    window$nt_start < 1L ||
    window$nt_end > read_length ||
    window$nt_start > window$nt_end
  ) {

    return(
      empty_result_row(

        pair_id =
          pair_id,

        mate =
          mate,

        orientation_role =
          orientation_role,

        target_position =
          position,

        hsp_strand =
          ifelse(
            isTRUE(window$forward),
            "forward",
            "reverse"
          ),

        read_start_coord =
          window$nt_start,

        read_length =
          read_length,

        codon_nt_start =
          window$nt_start,

        codon_nt_end =
          window$nt_end,

        note =
          "codon window outside read bounds"
      )
    )
  }


  # --------------------------------------------------------------------------
  # Extract nucleotide codon directly from FASTQ
  # --------------------------------------------------------------------------

  raw_codon_nt <- Biostrings::subseq(
    seq_obj,
    start =
      window$nt_start,
    end =
      window$nt_end
  )


  # --------------------------------------------------------------------------
  # Extract Phred scores directly from FASTQ
  # --------------------------------------------------------------------------

  raw_codon_qual <- as.integer(
    qual_obj[
      window$nt_start:
        window$nt_end
    ]
  )


  if (
    length(raw_codon_qual) != 3L
  ) {

    return(
      empty_result_row(

        pair_id =
          pair_id,

        mate =
          mate,

        orientation_role =
          orientation_role,

        target_position =
          position,

        hsp_strand =
          ifelse(
            isTRUE(window$forward),
            "forward",
            "reverse"
          ),

        read_start_coord =
          window$nt_start,

        read_length =
          read_length,

        codon_nt_start =
          window$nt_start,

        codon_nt_end =
          window$nt_end,

        note =
          "failed to obtain three Phred scores"
      )
    )
  }


  # --------------------------------------------------------------------------
  # Convert reverse-strand reads into biological orientation.
  #
  # IMPORTANT:
  # Quality scores must be reversed together with the nucleotide positions.
  # --------------------------------------------------------------------------

  if (
    isTRUE(window$forward)
  ) {

    final_codon <- raw_codon_nt

    final_quality <- raw_codon_qual

  } else {

    final_codon <-
      Biostrings::reverseComplement(
        raw_codon_nt
      )

    final_quality <-
      rev(
        raw_codon_qual
      )
  }


  # --------------------------------------------------------------------------
  # Convert codon to character string
  # --------------------------------------------------------------------------

  codon_string <- toupper(
    as.character(
      final_codon
    )
  )


  if (
    nchar(codon_string) != 3L
  ) {

    return(
      empty_result_row(

        pair_id =
          pair_id,

        mate =
          mate,

        orientation_role =
          orientation_role,

        target_position =
          position,

        hsp_strand =
          ifelse(
            isTRUE(window$forward),
            "forward",
            "reverse"
          ),

        read_start_coord =
          window$nt_start,

        read_length =
          read_length,

        codon_nt_start =
          window$nt_start,

        codon_nt_end =
          window$nt_end,

        codon =
          codon_string,

        note =
          "extracted codon is not three nucleotides"
      )
    )
  }


  # --------------------------------------------------------------------------
  # Translate internal TAC102 codon
  # --------------------------------------------------------------------------

  translated_aa <- translate_codon(
    codon_string
  )


  # --------------------------------------------------------------------------
  # Individual codon bases
  # --------------------------------------------------------------------------

  codon_base1 <- substr(
    codon_string,
    1L,
    1L
  )

  codon_base2 <- substr(
    codon_string,
    2L,
    2L
  )

  codon_base3 <- substr(
    codon_string,
    3L,
    3L
  )


  # --------------------------------------------------------------------------
  # Final result row
  # --------------------------------------------------------------------------

  empty_result_row(

    pair_id =
      pair_id,

    mate =
      mate,

    orientation_role =
      orientation_role,

    target_position =
      position,

    hsp_strand =
      ifelse(
        isTRUE(window$forward),
        "forward",
        "reverse"
      ),

    read_start_coord =
      window$nt_start,

    read_length =
      read_length,

    codon_nt_start =
      window$nt_start,

    codon_nt_end =
      window$nt_end,

    codon =
      codon_string,

    codon_base1 =
      codon_base1,

    codon_base2 =
      codon_base2,

    codon_base3 =
      codon_base3,

    codon_qual_pos1 =
      final_quality[1L],

    codon_qual_pos2 =
      final_quality[2L],

    codon_qual_pos3 =
      final_quality[3L],

    codon_min_qual =
      min(
        final_quality
      ),

    codon_mean_qual =
      round(
        mean(
          final_quality
        ),
        2
      ),

    translated_aa =
      translated_aa,

    note =
      ""
  )
}


# ==============================================================================
# 13. RUN BASE-LEVEL EXTRACTION
# ==============================================================================

cat(
  "Extracting base-level evidence from recovered FASTQ reads...\n\n"
)


result_list <- vector(
  "list",
  nrow(candidates)
)


for (
  i in seq_len(
    nrow(candidates)
  )
) {

  result_list[[i]] <- tryCatch(

    extract_base_evidence(
      candidates[
        i,
        ,
        drop = FALSE
      ]
    ),

    error = function(e) {

      empty_result_row(

        pair_id =
          as.character(
            candidates$pair_id[i]
          ),

        mate =
          as.character(
            candidates$mate[i]
          ),

        orientation_role =
          as.character(
            candidates$orientation_role[i]
          ),

        target_position =
          as.integer(
            candidates$target_position[i]
          ),

        note =
          paste0(
            "error during extraction: ",
            conditionMessage(e)
          )
      )
    }
  )


  if (
    i %% 10L == 0L ||
    i == nrow(candidates)
  ) {

    cat(
      "  Processed ",
      i,
      "/",
      nrow(candidates),
      " candidate records.\n",
      sep = ""
    )
  }
}


# ==============================================================================
# 14. COMBINE RESULTS
# ==============================================================================

results <- do.call(
  rbind,
  result_list
)


row.names(results) <- NULL


# ==============================================================================
# 15. ENFORCE OUTPUT COLUMN TYPES
# ==============================================================================

integer_columns <- c(
  "target_position",
  "read_start_coord",
  "read_length",
  "codon_nt_start",
  "codon_nt_end",
  "codon_qual_pos1",
  "codon_qual_pos2",
  "codon_qual_pos3",
  "codon_min_qual"
)


for (
  column in integer_columns
) {

  if (
    column %in% names(results)
  ) {

    results[[column]] <-
      suppressWarnings(
        as.integer(
          results[[column]]
        )
      )
  }
}


results$codon_mean_qual <-
  suppressWarnings(
    as.numeric(
      results$codon_mean_qual
    )
  )


# ==============================================================================
# 16. FINAL TRANSLATION CONSISTENCY CHECK
#
# Every extracted codon is translated again independently and compared with
# the stored translated_aa value.
# ==============================================================================

valid_translation_rows <- which(

  !is.na(
    results$codon
  ) &

    results$codon != "" &

    !is.na(
      results$translated_aa
    )
)


if (
  length(valid_translation_rows) > 0L
) {

  independently_translated <- vapply(

    results$codon[
      valid_translation_rows
    ],

    translate_codon,

    character(1)
  )


  stored_translation <- as.character(

    results$translated_aa[
      valid_translation_rows
    ]
  )


  mismatch_indices <- which(

    independently_translated !=
      stored_translation
  )


  if (
    length(mismatch_indices) > 0L
  ) {

    bad_rows <- valid_translation_rows[
      mismatch_indices
    ]


    stop(
      paste0(
        "Final translation consistency check FAILED.\n",
        "Rows with translation mismatch: ",
        length(bad_rows),
        "\n",
        "First affected rows: ",
        paste(
          head(
            bad_rows,
            10L
          ),
          collapse = ", "
        )
      )
    )
  }
}


cat(
  "Final translation consistency check: PASSED\n\n"
)


# ==============================================================================
# 17. SAVE DETAILED OUTPUT
# ==============================================================================

write.csv(
  results,
  DETAIL_OUT,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


cat(
  "Detailed validation output saved:\n  ",
  DETAIL_OUT,
  "\n\n",
  sep = ""
)


# ==============================================================================
# 18. POSITION / ALLELE SUMMARY
#
# Summary is deliberately constructed row-by-row rather than through a chain
# of aggregate()/merge() calls. This guarantees stable column structure.
# ==============================================================================

usable_results <- results[
  results$note == "" &
    !is.na(results$translated_aa) &
    results$translated_aa != "X",
  ,
  drop = FALSE
]


summary_rows <- list()


for (
  position in TARGET_POSITIONS
) {

  df <- usable_results[
    usable_results$target_position == position,
    ,
    drop = FALSE
  ]


  if (
    nrow(df) == 0L
  ) {

    next
  }


  alleles <- sort(
    unique(
      as.character(
        df$translated_aa
      )
    )
  )


  for (
    allele in alleles
  ) {

    allele_df <- df[
      df$translated_aa == allele,
      ,
      drop = FALSE
    ]


    family_values <- allele_df$dup_family

    family_values <- family_values[
      !is.na(family_values) &
        family_values != ""
    ]


    independent_family_count <-
      length(
        unique(
          family_values
        )
      )


    min_quality <- allele_df$codon_min_qual

    min_quality <- min_quality[
      !is.na(min_quality)
    ]


    mean_quality <- allele_df$codon_min_qual

    mean_quality <- mean_quality[
      !is.na(mean_quality)
    ]


    if (
      length(min_quality) == 0L
    ) {

      minimum_min_phred <- NA_integer_

    } else {

      minimum_min_phred <-
        as.integer(
          min(
            min_quality
          )
        )
    }


    if (
      length(mean_quality) == 0L
    ) {

      mean_min_phred <- NA_real_

    } else {

      mean_min_phred <-
        round(
          mean(
            mean_quality
          ),
          2
        )
    }


    summary_rows[[length(summary_rows) + 1L]] <- data.frame(

      target_position =
        as.integer(
          position
        ),

      translated_aa =
        as.character(
          allele
        ),

      read_count =
        as.integer(
          nrow(allele_df)
        ),

      mean_min_phred =
        as.numeric(
          mean_min_phred
        ),

      minimum_min_phred =
        as.integer(
          minimum_min_phred
        ),

      independent_families =
        as.integer(
          independent_family_count
        ),

      stringsAsFactors = FALSE
    )
  }
}


if (
  length(summary_rows) > 0L
) {

  summary_df <- do.call(
    rbind,
    summary_rows
  )

} else {

  summary_df <- data.frame(

    target_position =
      integer(0),

    translated_aa =
      character(0),

    read_count =
      integer(0),

    mean_min_phred =
      numeric(0),

    minimum_min_phred =
      integer(0),

    independent_families =
      integer(0),

    stringsAsFactors = FALSE
  )
}


summary_df <- summary_df[
  order(
    summary_df$target_position,
    summary_df$translated_aa
  ),
  ,
  drop = FALSE
]


write.csv(
  summary_df,
  SUMMARY_OUT,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


cat(
  "Position/allele summary saved:\n  ",
  SUMMARY_OUT,
  "\n\n",
  sep = ""
)


# ==============================================================================
# 19. DUPLICATE / READ-START FAMILY CHECK
#
# This is a HEURISTIC.
#
# Sharing a read-start coordinate does NOT prove PCR duplication.
#
# Reads with the same:
#
#     mate
#     strand
#     approximate start coordinate
#
# are grouped into the same potential read-start family.
# ==============================================================================

flag_duplicate_families <- function(
    df) {

  df$dup_family <- NA_character_


  if (
    nrow(df) == 0L ||
    !"read_start_coord" %in% names(df)
  ) {

    return(df)
  }


  valid <- (

    !is.na(
      df$read_start_coord
    ) &

      !is.na(
        df$mate
      ) &

      !is.na(
        df$hsp_strand
      )
  )


  if (
    !any(valid)
  ) {

    return(df)
  }


  sub <- df[
    valid,
    ,
    drop = FALSE
  ]


  sub <- sub[
    order(
      sub$mate,
      sub$hsp_strand,
      sub$read_start_coord
    ),
    ,
    drop = FALSE
  ]


  family_number <- 0L


  family_ids <- character(
    nrow(sub)
  )


  for (
    i in seq_len(
      nrow(sub)
    )
  ) {

    new_family <- (

      i == 1L ||

        sub$mate[i] !=
          sub$mate[i - 1L] ||

        sub$hsp_strand[i] !=
          sub$hsp_strand[i - 1L] ||

        abs(
          sub$read_start_coord[i] -
            sub$read_start_coord[i - 1L]
        ) > DUP_WINDOW
    )


    if (
      new_family
    ) {

      family_number <-
        family_number + 1L
    }


    family_ids[i] <-
      paste0(
        sub$mate[i],
        "_",
        sub$hsp_strand[i],
        "_F",
        family_number
      )
  }


  sub$dup_family <-
    family_ids


  # --------------------------------------------------------------------------
  # Restore family IDs to original row order.
  #
  # pair_id + mate + target_position uniquely identify a candidate record.
  # --------------------------------------------------------------------------

  key_original <- paste(

    df$pair_id,
    df$mate,
    df$target_position,

    sep = "|"
  )


  key_sub <- paste(

    sub$pair_id,
    sub$mate,
    sub$target_position,

    sep = "|"
  )


  match_index <- match(
    key_original[valid],
    key_sub
  )


  df$dup_family[valid] <-
    sub$dup_family[
      match_index
    ]


  df
}


# Apply family assignment separately within each target position.
results_split <- split(
  results,
  results$target_position
)


results_split <- lapply(
  results_split,
  flag_duplicate_families
)


results <- do.call(
  rbind,
  results_split
)


row.names(results) <- NULL


# Re-save detailed output now that duplicate families have been assigned.
write.csv(
  results,
  DETAIL_OUT,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


# ==============================================================================
# 20. REBUILD SUMMARY AFTER DUPLICATE FAMILY ASSIGNMENT
#
# This ensures independent_families in the final summary reflects the actual
# family assignments.
# ==============================================================================

usable_results <- results[
  results$note == "" &
    !is.na(results$translated_aa) &
    results$translated_aa != "X",
  ,
  drop = FALSE
]


summary_rows <- list()


for (
  position in TARGET_POSITIONS
) {

  df <- usable_results[
    usable_results$target_position == position,
    ,
    drop = FALSE
  ]


  if (
    nrow(df) == 0L
  ) {

    next
  }


  alleles <- sort(
    unique(
      as.character(
        df$translated_aa
      )
    )
  )


  for (
    allele in alleles
  ) {

    allele_df <- df[
      df$translated_aa == allele,
      ,
      drop = FALSE
    ]


    family_values <- allele_df$dup_family

    family_values <- family_values[
      !is.na(family_values) &
        family_values != ""
    ]


    quality_values <-
      allele_df$codon_min_qual

    quality_values <-
      quality_values[
        !is.na(
          quality_values
        )
      ]


    if (
      length(quality_values) == 0L
    ) {

      mean_min_phred <- NA_real_

      minimum_min_phred <- NA_integer_

    } else {

      mean_min_phred <-
        round(
          mean(
            quality_values
          ),
          2
        )

      minimum_min_phred <-
        as.integer(
          min(
            quality_values
          )
        )
    }


    summary_rows[[length(summary_rows) + 1L]] <- data.frame(

      target_position =
        as.integer(
          position
        ),

      translated_aa =
        as.character(
          allele
        ),

      read_count =
        as.integer(
          nrow(allele_df)
        ),

      mean_min_phred =
        as.numeric(
          mean_min_phred
        ),

      minimum_min_phred =
        as.integer(
          minimum_min_phred
        ),

      independent_families =
        as.integer(
          length(
            unique(
              family_values
            )
          )
        ),

      stringsAsFactors = FALSE
    )
  }
}


if (
  length(summary_rows) > 0L
) {

  summary_df <- do.call(
    rbind,
    summary_rows

  )

} else {

  summary_df <- data.frame(

    target_position =
      integer(0),

    translated_aa =
      character(0),

    read_count =
      integer(0),

    mean_min_phred =
      numeric(0),

    minimum_min_phred =
      integer(0),

    independent_families =
      integer(0),

    stringsAsFactors = FALSE
  )
}


summary_df <- summary_df[
  order(
    summary_df$target_position,
    summary_df$translated_aa
  ),
  ,
  drop = FALSE
]


write.csv(
  summary_df,
  SUMMARY_OUT,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)


# ==============================================================================
# 21. POSITION SUMMARY STATISTICS
# ==============================================================================

make_position_summary <- function(
    position) {

  df <- results[
    results$target_position == position &
      results$note == "",
    ,
    drop = FALSE
  ]


  if (
    nrow(df) == 0L
  ) {

    return(
      list(

        candidate_records =
          sum(
            candidates$target_position ==
              position
          ),

        usable_records =
          0L,

        independent_read_start_families =
          0L,

        mean_codon_min_phred =
          NA_real_,

        minimum_codon_min_phred =
          NA_integer_
      )
    )
  }


  usable <- df[
    !is.na(df$translated_aa) &
      df$translated_aa != "X",
    ,
    drop = FALSE
  ]


  quality_values <-
    usable$codon_min_qual

  quality_values <-
    quality_values[
      !is.na(
        quality_values
      )
    ]


  family_values <-
    usable$dup_family

  family_values <-
    family_values[
      !is.na(
        family_values
      ) &
        family_values != ""
    ]


  if (
    length(quality_values) == 0L
  ) {

    mean_quality <- NA_real_

    minimum_quality <- NA_integer_

  } else {

    mean_quality <-
      round(
        mean(
          quality_values
        ),
        2
      )

    minimum_quality <-
      as.integer(
        min(
          quality_values
        )
      )
  }


  list(

    candidate_records =
      sum(
        candidates$target_position ==
          position
      ),

    usable_records =
      nrow(usable),

    independent_read_start_families =
      length(
        unique(
          family_values
        )
      ),

    mean_codon_min_phred =
      mean_quality,

    minimum_codon_min_phred =
      minimum_quality
  )
}


summary_653 <- make_position_summary(
  653L
)


summary_698 <- make_position_summary(
  698L
)


# ==============================================================================
# 22. CROSS-POSITION SHARED READ-PAIR CHECK
# ==============================================================================

ids_653 <- unique(
  as.character(
    candidates$pair_id[
      candidates$target_position ==
        653L
    ]
  )
)


ids_698 <- unique(
  as.character(
    candidates$pair_id[
      candidates$target_position ==
        698L
    ]
  )
)


shared_pair_ids <- intersect(
  ids_653,
  ids_698
)


cat(
  "\n==============================================================\n",
  " CROSS-POSITION READ-PAIR INDEPENDENCE\n",
  "==============================================================\n\n",
  sep = ""
)


cat(
  "Read pairs contributing to both positions: ",
  length(shared_pair_ids),
  "\n",
  sep = ""
)


if (
  length(shared_pair_ids) > 0L
) {

  cat(
    "\nShared read-pair IDs:\n"
  )

  print(
    shared_pair_ids
  )

} else {

  cat(
    "\nNo read pairs contribute to both positions.\n",
    "The two candidate sites are supported by distinct\n",
    "candidate physical fragments.\n",
    sep = ""
  )
}


# ==============================================================================
# 23. ALLELE / FAMILY STRUCTURE
# ==============================================================================

cat(
  "\n==============================================================\n",
  " ALLELE / FAMILY STRUCTURE\n",
  "==============================================================\n",
  sep = ""
)


for (
  position in TARGET_POSITIONS
) {

  cat(
    "\nPosition ",
    position,
    ":\n",
    sep = ""
  )


  df <- results[
    results$target_position == position &
      results$note == "" &
      !is.na(results$translated_aa) &
      results$translated_aa != "X",
    ,
    drop = FALSE
  ]


  if (
    nrow(df) == 0L
  ) {

    cat(
      "  No usable nucleotide evidence.\n"
    )

    next
  }


  family_table <- aggregate(
    dup_family ~ translated_aa,
    data = df,
    FUN = function(x) {

      length(
        unique(
          x[
            !is.na(x) &
              x != ""
          ]
        )
      )
    }
  )


  names(family_table)[2] <-
    "independent_families"


  print(
    family_table,
    row.names = FALSE
  )
}


# ==============================================================================
# 24. TEXT REPORT
# ==============================================================================

report_lines <- c(

  "==============================================================",

  "TAC102 DIRECT FASTQ BASE-LEVEL VALIDATION",

  "Step 4B.3l",

  "==============================================================",

  "",

  paste0(
    "Date: ",
    Sys.Date()
  ),

  paste0(
    "Time: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  "",

  "SOFTWARE / TRANSLATION",

  paste0(
    "  R version: ",
    R.version.string
  ),

  paste0(
    "  Biostrings version: ",
    as.character(
      packageVersion(
        "Biostrings"
      )
    )
  ),

  "",

  "  Internal TAC102 codons were translated using:",

  "    Biostrings::translate(..., no.init.codon=TRUE)",

  "",

  "  Translation sanity check:",

  "    ATG -> M",

  "    CTG -> L",

  "    TTG -> L",

  "    GTG -> V",

  "    AAA -> K",

  "    TTT -> F",

  "    GGG -> G",

  "",

  "METHOD",

  "  Candidate reads were identified using the TAC102 unaligned-tail",

  "  candidate table and original BLAST HSP coordinates.",

  "",

  "  Nucleotide coordinates were independently re-derived from:",

  "    hsp_qstart_raw",

  "    hsp_sstart_raw",

  "    hsp_send_raw",

  "",

  "  Sequence and Phred quality were re-extracted directly from the",

  "  ORIGINAL R1/R2 FASTQ files.",

  "",

  "  The following reconstructed values were NOT used to make the",

  "  nucleotide call:",

  "    gap_facing_sequence",

  "    reconstructed_codon",

  "    reconstructed_aa",

  "    translated_candidate",

  "    4B.3j reconstructed values",

  "    4B.3k provenance nucleotide values",

  "",

  "POSITION 653",

  paste0(
    "  Candidate records: ",
    summary_653$candidate_records
  ),

  paste0(
    "  Usable records: ",
    summary_653$usable_records
  ),

  paste0(
    "  Independent read-start families: ",
    summary_653$independent_read_start_families
  ),

  paste0(
    "  Mean codon minimum Phred: ",
    summary_653$mean_codon_min_phred
  ),

  paste0(
    "  Minimum codon minimum Phred: ",
    summary_653$minimum_codon_min_phred
  ),

  "",

  "POSITION 698",

  paste0(
    "  Candidate records: ",
    summary_698$candidate_records
  ),

  paste0(
    "  Usable records: ",
    summary_698$usable_records
  ),

  paste0(
    "  Independent read-start families: ",
    summary_698$independent_read_start_families
  ),

  paste0(
    "  Mean codon minimum Phred: ",
    summary_698$mean_codon_min_phred
  ),

  paste0(
    "  Minimum codon minimum Phred: ",
    summary_698$minimum_codon_min_phred
  ),

  "",

  "CROSS-POSITION INDEPENDENCE",

  paste0(
    "  Read pairs contributing to BOTH positions: ",
    length(shared_pair_ids)
  ),

  "",

  "ALLELE SUMMARY",

  "  See tac102_4B3l_position_validation_summary.csv",

  "",

  "INTERPRETATION GUIDE",

  "  This script does not automatically declare a biological SNP.",

  "",

  "  Strong candidate evidence would consist of:",

  "    - reproducible alternative base/codon calls",

  "    - high Phred quality",

  "    - multiple independent read-start families",

  "    - consistent orientation handling",

  "    - no concentration in a single duplicate family",

  "    - no evidence that the signal is restricted to one problematic read",

  "",

  "  Low-quality or isolated alternative calls should remain suspect.",

  "",

  "OUTPUTS",

  paste0(
    "  Detailed per-read validation: ",
    DETAIL_OUT
  ),

  paste0(
    "  Position/allele summary: ",
    SUMMARY_OUT
  ),

  paste0(
    "  Text report: ",
    REPORT_OUT
  ),

  "",

  "=============================================================="
)


writeLines(
  report_lines,
  REPORT_OUT
)


# ==============================================================================
# 25. FINAL LOGGING
# ==============================================================================

log_info(
  paste0(
    "Step 4B.3l direct FASTQ base-level validation ",
    "completed successfully."
  ),
  LOG_FILE
)


file.copy(
  LOG_FILE,
  file.path(
    PATHS$logs,
    "latest_execution.log"
  ),
  overwrite = TRUE
)


# ==============================================================================
# 26. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n==============================================================\n",
  " Step 4B.3l COMPLETED SUCCESSFULLY\n",
  "==============================================================\n\n",
  sep = ""
)


cat(
  "Translation validation: PASSED\n",
  "  ATG -> M\n",
  "  CTG -> L\n",
  "  TTG -> L\n",
  "  GTG -> V\n",
  "  AAA -> K\n",
  "  TTT -> F\n",
  "  GGG -> G\n\n",
  sep = ""
)


cat(
  "Candidate records:\n",
  "  Position 653: ",
  n_653,
  "\n",
  "  Position 698: ",
  n_698,
  "\n\n",
  sep = ""
)


cat(
  "Usable records:\n",
  "  Position 653: ",
  summary_653$usable_records,
  "\n",
  "  Position 698: ",
  summary_698$usable_records,
  "\n\n",
  sep = ""
)


cat(
  "Independent read-start families:\n",
  "  Position 653: ",
  summary_653$independent_read_start_families,
  "\n",
  "  Position 698: ",
  summary_698$independent_read_start_families,
  "\n\n",
  sep = ""
)


cat(
  "Shared physical read pairs across positions: ",
  length(shared_pair_ids),
  "\n\n",
  sep = ""
)


cat(
  "Output files:\n",
  "  Detailed: ",
  DETAIL_OUT,
  "\n",
  "  Summary:  ",
  SUMMARY_OUT,
  "\n",
  "  Report:   ",
  REPORT_OUT,
  "\n\n",
  sep = ""
)


cat(
  "[SUCCESS] Corrected TAC102 Step 4B.3l completed.\n"
)