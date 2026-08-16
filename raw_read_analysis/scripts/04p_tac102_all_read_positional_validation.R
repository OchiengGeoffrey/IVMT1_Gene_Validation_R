# ==============================================================================
# Step 4B.3p: TAC102 All-Read Positional Validation
# Positions 653 and 698
#
# PURPOSE
# -------
# Independently determine which TAC102 alleles are present among ALL sequencing
# reads capable of covering protein positions 653 and/or 698.
#
# IMPORTANT DESIGN PRINCIPLE
# --------------------------
# This analysis is intentionally independent of the previous candidate-read
# recovery pipeline.
#
# It does NOT use:
#   - BLAST HSP tables
#   - raw_read_recovery_* files
#   - tac102_41pairs_* files
#   - gap_facing_sequence
#   - translated_candidate
#   - 4B.3j outputs
#   - 4B.3k outputs
#   - 4B.3l outputs
#   - 4B.3m outputs
#   - 4B.3n outputs
#   - 4B.3o outputs
#
# INPUTS
# -------
#   1. TAC102 reference sequence from reference_sequences.rds
#   2. Original R1 FASTQ
#   3. Original R2 FASTQ
#
# METHOD
# ------
# Every sequencing read is independently examined in all six reading frames:
#
#   forward frames:      +1, +2, +3
#   reverse-complement:  -1, -2, -3
#
# A short TAC102 reference peptide window surrounding the target position is
# searched directly within each translated read.
#
# The target amino acid itself is EXCLUDED from the anchor identity calculation
# so that a genuine variant at position 653 or 698 does not prevent its own
# detection.
#
# BEFORE ANY REAL FASTQ PROCESSING
# --------------------------------
# A synthetic known-answer test suite must pass.
#
# The script will STOP if any pre-flight test fails.
#
# PILOT MODE
# ----------
# MAX_READS controls the maximum number of reads processed per FASTQ.
#
# Examples:
#
#   MAX_READS <- 100000L
#       Process first 100,000 R1 reads and first 100,000 R2 reads.
#
#   MAX_READS <- NA_integer_
#       Process the complete FASTQ files.
#
# IMPORTANT:
# MAX_READS is applied independently to R1 and R2.
#
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

if (!exists("PATHS") ||
    !exists("LOG_FILE") ||
    !exists("CONFIG")) {

  source(
    here::here(
      "raw_read_analysis",
      "scripts",
      "00_setup.R"
    )
  )
}


# ==============================================================================
# 1. USER-CONTROLLED PARAMETERS
# ==============================================================================

TARGET_POSITIONS <- c(653L, 698L)

# Number of amino acids on either side of target position used for anchoring.
#
# 8 aa on either side means:
#
#   653 -> reference positions 645:661
#   698 -> reference positions 690:706
#
# The target amino acid itself is excluded from identity scoring.
FLANK_WINDOW <- 8L

# Minimum anchor identity.
#
# Example:
# 17-aa anchor with target excluded = 16 informative residues.
#
# 0.75 therefore requires at least 12/16 matching residues.
MIN_ANCHOR_IDENTITY <- 0.75

# ----------------------------------------------------------------------
# PILOT CONTROL
# ----------------------------------------------------------------------
#
# Set to 100000L for the recommended first pilot.
#
# Set to NA_integer_ to process the complete FASTQ.
#
# MAX_READS applies independently to R1 and R2.
#
MAX_READS <- 100000L

# FASTQ streaming chunk size.
#
# 500,000 reads was used in the earlier design, but a smaller chunk is
# intentionally used here to keep memory behaviour conservative during
# the pilot.
STREAM_CHUNK_SIZE <- 100000L

# Progress reporting interval in number of reads.
PROGRESS_INTERVAL <- 1000000L

# Whether to retain all position-hit records.
#
# TRUE is recommended for the pilot and final analysis because the complete
# evidence table is scientifically useful.
SAVE_ALL_HITS <- TRUE


# ==============================================================================
# 2. REFERENCE
# ==============================================================================

reference_sequences <- readRDS(
  file.path(
    PATHS$intermediate,
    "reference_sequences.rds"
  )
)

if (!"TAC102" %in% names(reference_sequences)) {
  stop(
    "TAC102 was not found in reference_sequences.rds.\n",
    "Available names: ",
    paste(names(reference_sequences), collapse = ", ")
  )
}

tac_ref_chr <- as.character(
  reference_sequences[["TAC102"]]
)

if (length(tac_ref_chr) != 1L || is.na(tac_ref_chr)) {
  stop("TAC102 reference sequence is invalid.")
}

tac_ref_chr <- toupper(tac_ref_chr)

cat("\n")
cat("==============================================================\n")
cat("Step 4B.3p: TAC102 ALL-READ POSITIONAL VALIDATION\n")
cat("==============================================================\n\n")

cat(
  "TAC102 reference length:",
  nchar(tac_ref_chr),
  "aa\n"
)

cat(
  "Target positions:",
  paste(TARGET_POSITIONS, collapse = ", "),
  "\n"
)

cat(
  "Flank window:",
  FLANK_WINDOW,
  "aa\n"
)

cat(
  "Minimum anchor identity:",
  MIN_ANCHOR_IDENTITY,
  "\n"
)

if (is.na(MAX_READS)) {

  cat(
    "MAX_READS: FULL DATASET\n"
  )

} else {

  cat(
    "MAX_READS per FASTQ:",
    format(MAX_READS, big.mark = ","),
    "\n"
  )
}

cat(
  "Streaming chunk size:",
  format(STREAM_CHUNK_SIZE, big.mark = ","),
  "\n\n"
)


# ==============================================================================
# 3. DEPENDENCIES
# ==============================================================================

required_packages <- c(
  "Biostrings",
  "ShortRead"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {

  stop(
    "Required packages are missing: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running 4B.3p."
  )
}


# ==============================================================================
# 4. GENETIC CODE
# ==============================================================================

GENCODE <- Biostrings::getGeneticCode("1")


# ==============================================================================
# 5. OUTPUT FILES
# ==============================================================================

output_hits <- file.path(
  PATHS$reports,
  "tac102_4B3p_all_read_positional_evidence.csv"
)

output_summary <- file.path(
  PATHS$reports,
  "tac102_4B3p_all_read_positional_summary.csv"
)

output_haplotype <- file.path(
  PATHS$reports,
  "tac102_4B3p_all_read_haplotype_summary.csv"
)

output_report <- file.path(
  PATHS$reports,
  "tac102_4B3p_all_read_positional_report.txt"
)

output_preflight <- file.path(
  PATHS$reports,
  "tac102_4B3p_preflight_results.csv"
)


# ==============================================================================
# 6. BASIC VALIDATION OF TARGET POSITIONS
# ==============================================================================

if (any(TARGET_POSITIONS < 1L) ||
    any(TARGET_POSITIONS > nchar(tac_ref_chr))) {

  stop(
    "One or more TARGET_POSITIONS falls outside the TAC102 reference."
  )
}


# ==============================================================================
# 7. REFERENCE ANCHOR GENERATION
# ==============================================================================

make_anchor <- function(
    ref_chr,
    target_pos,
    flank = FLANK_WINDOW) {

  ref_len <- nchar(ref_chr)

  win_start <- max(
    1L,
    target_pos - flank
  )

  win_end <- min(
    ref_len,
    target_pos + flank
  )

  anchor <- substr(
    ref_chr,
    win_start,
    win_end
  )

  target_offset <- target_pos - win_start + 1L

  list(
    anchor = anchor,
    win_start = win_start,
    win_end = win_end,
    target_offset = target_offset
  )
}


anchor_653 <- make_anchor(
  tac_ref_chr,
  653L,
  FLANK_WINDOW
)

anchor_698 <- make_anchor(
  tac_ref_chr,
  698L,
  FLANK_WINDOW
)


cat(
  "653 anchor:",
  anchor_653$win_start,
  "-",
  anchor_653$win_end,
  "|",
  anchor_653$anchor,
  "\n"
)

cat(
  "698 anchor:",
  anchor_698$win_start,
  "-",
  anchor_698$win_end,
  "|",
  anchor_698$anchor,
  "\n\n"
)


# ==============================================================================
# 8. CODON TABLE FOR SYNTHETIC TESTS
# ==============================================================================

TEST_CODON_TABLE <- c(
  A = "GCT",
  R = "CGT",
  N = "AAT",
  D = "GAT",
  C = "TGT",
  Q = "CAA",
  E = "GAA",
  G = "GGT",
  H = "CAT",
  I = "ATT",
  L = "CTT",
  K = "AAA",
  M = "ATG",
  F = "TTT",
  P = "CCT",
  S = "TCT",
  T = "ACT",
  W = "TGG",
  Y = "TAT",
  V = "GTT",
  X = "NNN"
)


# ==============================================================================
# 9. TRANSLATION HELPER
# ==============================================================================

translate_read_frame <- function(
    dna,
    frame) {

    if (!inherits(dna, "DNAString")) {
    stop(
      "translate_read_frame() requires a DNAString."
    )
  }

  if (frame < 0L || frame > 2L) {
    stop(
      "frame must be 0, 1, or 2."
    )
  }

  dna_len <- length(dna)

  usable_len <- dna_len - frame

  if (usable_len < 3L) {
    return(NULL)
  }

  usable_len <- usable_len -
    (usable_len %% 3L)

  if (usable_len < 3L) {
    return(NULL)
  }

  frame_seq <- Biostrings::subseq(
    dna,
    start = frame + 1L,
    width = usable_len
  )

  pep <- tryCatch(

    Biostrings::translate(
      frame_seq,
      genetic.code = GENCODE,
      no.init.codon = TRUE,
      if.fuzzy.codon = "X"
    ),

    error = function(e) {
      NULL
    }
  )

  if (is.null(pep)) {
    return(NULL)
  }

  list(
    peptide = as.character(pep),
    frame_seq = frame_seq,
    frame = frame
  )
}


# ==============================================================================
# 10. POSITION LOCATOR
# ==============================================================================

locate_position_in_read <- function(
    read_seq,
    target_pos,
    ref_chr,
    flank = FLANK_WINDOW,
    min_identity = MIN_ANCHOR_IDENTITY) {

  # --------------------------------------------------------------------
  # HARD TYPE GUARD
  # --------------------------------------------------------------------

    if (!inherits(read_seq, "DNAString")) {

    stop(
      "locate_position_in_read(): read_seq must be a DNAString.\n",
      "This guard intentionally prevents DNAStringSet/DNAString indexing bugs."
    )
  }

  # --------------------------------------------------------------------
  # REFERENCE WINDOW
  # --------------------------------------------------------------------

  ref_len <- nchar(ref_chr)

  win_start <- max(
    1L,
    target_pos - flank
  )

  win_end <- min(
    ref_len,
    target_pos + flank
  )

  anchor_peptide <- substr(
    ref_chr,
    win_start,
    win_end
  )

  target_offset_in_anchor <-
    target_pos - win_start + 1L

  anchor_chars <- strsplit(
    anchor_peptide,
    "",
    fixed = TRUE
  )[[1L]]

  # Target excluded from anchor scoring.
  compare_idx <- setdiff(
    seq_along(anchor_chars),
    target_offset_in_anchor
  )

  read_len <- length(read_seq)

  if (read_len < 3L) {

    return(
      list(
        found = FALSE,
        reason = "read too short"
      )
    )
  }

  # --------------------------------------------------------------------
  # SIX READING FRAMES
  # --------------------------------------------------------------------

  orientations <- list(

    fwd = read_seq,

    rev = Biostrings::reverseComplement(
      read_seq
    )

  )

  best <- list(
    found = FALSE,
    score = -Inf
  )

  for (orientation_name in names(orientations)) {

    oriented_seq <- orientations[[orientation_name]]

    for (frame in 0:2) {

      translated <- translate_read_frame(
        oriented_seq,
        frame
      )

      if (is.null(translated)) {
        next
      }

      pep <- translated$peptide

      if (nchar(pep) < nchar(anchor_peptide)) {
        next
      }

      n_positions <-
        nchar(pep) -
        nchar(anchor_peptide) +
        1L

      if (n_positions < 1L) {
        next
      }

      # --------------------------------------------------------------
      # SLIDING ANCHOR SEARCH
      # --------------------------------------------------------------

      for (start_pos in seq_len(n_positions)) {

        candidate_window <- substr(
          pep,
          start_pos,
          start_pos + nchar(anchor_peptide) - 1L
        )

        candidate_chars <- strsplit(
          candidate_window,
          "",
          fixed = TRUE
        )[[1L]]

        if (length(candidate_chars) != length(anchor_chars)) {
          next
        }

        matches <- candidate_chars[compare_idx] ==
          anchor_chars[compare_idx]

        identity <- mean(
          matches
        )

        if (is.na(identity)) {
          next
        }

        # ------------------------------------------------------------
        # ACCEPT ONLY IF ABOVE THRESHOLD
        # ------------------------------------------------------------

        if (identity >= min_identity &&
            identity > best$score) {

          target_aa <- candidate_chars[
            target_offset_in_anchor
          ]

          # Position of target codon within the translated frame.
          #
          # start_pos is 1-based amino-acid position.
          #
          # target_offset_in_anchor is 1-based amino-acid offset.
          #
          # Therefore:
          #
          # target_aa_position =
          #   start_pos + target_offset - 1
          #
          # Codon start within frame_seq:
          #
          # ((aa_position - 1) * 3) + 1
          #
          target_aa_position <-
            start_pos +
            target_offset_in_anchor -
            1L

          codon_start_in_frame <-
            ((target_aa_position - 1L) * 3L) + 1L

          codon_end_in_frame <-
            codon_start_in_frame + 2L

          if (codon_end_in_frame >
              length(translated$frame_seq)) {

            next
          }

          target_codon <- as.character(
            Biostrings::subseq(
              translated$frame_seq,
              start = codon_start_in_frame,
              width = 3L
            )
          )

          best <- list(

            found = TRUE,

            score = identity,

            orientation = orientation_name,

            frame = frame,

            target_aa = target_aa,

            target_aa_position = target_aa_position,

            codon_start_in_frame =
              codon_start_in_frame,

            codon_end_in_frame =
              codon_end_in_frame,

            target_codon =
              target_codon,

            anchor_start_reference =
              win_start,

            anchor_end_reference =
              win_end,

            anchor_peptide =
              anchor_peptide,

            anchor_candidate =
              candidate_window

          )
        }
      }
    }
  }

  if (!best$found) {

    return(
      list(
        found = FALSE,
        reason = "no anchor match above threshold"
      )
    )
  }

  best
}


# ==============================================================================
# 11. CONVERT FRAME-RELATIVE CODON TO ORIGINAL READ COORDINATES
# ==============================================================================

get_original_read_codon_coordinates <- function(
    read_length,
    orientation,
    frame,
    codon_start_in_frame) {

  # --------------------------------------------------------------
  # FORWARD ORIENTATION
  # --------------------------------------------------------------

  if (orientation == "fwd") {

    read_start <-
      frame +
      codon_start_in_frame

    read_end <-
      read_start + 2L

    return(
      c(
        start = read_start,
        end = read_end
      )
    )
  }

  # --------------------------------------------------------------
  # REVERSE-COMPLEMENT ORIENTATION
  # --------------------------------------------------------------
  #
  # The translated sequence is based on reverseComplement(read).
  #
  # A codon occupying positions:
  #
  #   codon_start_in_frame : codon_start_in_frame + 2
  #
  # in the reverse-complement corresponds to:
  #
  #   read_length - codon_end + 1
  #   through
  #   read_length - codon_start + 1
  #
  # in the original read.
  # --------------------------------------------------------------

  if (orientation == "rev") {

    rc_start <-
      frame +
      codon_start_in_frame

    rc_end <-
      rc_start + 2L

    read_start <-
      read_length -
      rc_end +
      1L

    read_end <-
      read_length -
      rc_start +
      1L

    return(
      c(
        start = read_start,
        end = read_end
      )
    )
  }

  stop(
    "Unknown orientation: ",
    orientation
  )
}


# ==============================================================================
# 12. PHRED EXTRACTION
# ==============================================================================

extract_codon_quality <- function(
    quality_vector,
    start,
    end,
    orientation) {

  if (start < 1L ||
      end > length(quality_vector) ||
      start > end) {

    stop(
      "Invalid quality-coordinate request: ",
      start,
      "-",
      end,
      " for quality vector of length ",
      length(quality_vector)
    )
  }

  q <- as.integer(
    quality_vector[start:end]
  )

  # The nucleotide order is reversed in the reverse-complement representation.
  #
  # For reporting the original-read codon quality, preserve the original
  # read order. The minimum is invariant to reversal, but the vector itself
  # is useful for debugging.
  if (orientation == "rev") {
    q_for_codon <- rev(q)
  } else {
    q_for_codon <- q
  }

  list(
    values = q_for_codon,
    min = min(q_for_codon),
    mean = mean(q_for_codon)
  )
}


# ==============================================================================
# 13. SYNTHETIC READ BUILDER
# ==============================================================================

make_synthetic_read <- function(
    center_pos,
    offset_from_read_start,
    read_nt_length = 90L,
    reverse = FALSE,
    mutate_pos = NULL,
    mutate_to = NULL) {

  # --------------------------------------------------------------
  # Convert requested read length into an amino-acid window.
  # --------------------------------------------------------------

  aa_length <- floor(
    read_nt_length / 3L
  )

  # offset_from_read_start is interpreted as amino-acid offset:
  # read_start = target - offset
  win_start <-
    center_pos -
    offset_from_read_start

  win_start <- max(
    1L,
    win_start
  )

  win_end <-
    win_start +
    aa_length -
    1L

  win_end <- min(
    nchar(tac_ref_chr),
    win_end
  )

  pep <- substr(
    tac_ref_chr,
    win_start,
    win_end
  )

  if (!is.null(mutate_pos)) {

    idx <-
      mutate_pos -
      win_start +
      1L

    if (idx < 1L ||
        idx > nchar(pep)) {

      stop(
        "Requested synthetic mutation is outside synthetic peptide."
      )
    }

    chars <- strsplit(
      pep,
      "",
      fixed = TRUE
    )[[1L]]

    chars[idx] <- mutate_to

    pep <- paste(
      chars,
      collapse = ""
    )
  }

  aa <- strsplit(
    pep,
    "",
    fixed = TRUE
  )[[1L]]

  if (any(is.na(TEST_CODON_TABLE[aa]))) {

    stop(
      "Synthetic peptide contains amino acid without test codon."
    )
  }

  nt <- paste(
    TEST_CODON_TABLE[aa],
    collapse = ""
  )

  seq <- Biostrings::DNAString(
    nt
  )

  if (reverse) {

    seq <- Biostrings::reverseComplement(
      seq
    )
  }

  expected_aa <- substr(
    pep,
    center_pos - win_start + 1L,
    center_pos - win_start + 1L
  )

  expected_codon <- TEST_CODON_TABLE[
    expected_aa
  ]

  list(

    seq = seq,

    expected_start = win_start,

    expected_end = win_end,

    expected_aa =
      expected_aa,

    expected_codon =
      unname(expected_codon)

  )
}


# ==============================================================================
# 14. SYNTHETIC QUALITY VECTOR BUILDER
# ==============================================================================

make_test_quality <- function(
    read_length,
    high_quality = 35L) {

  rep(
    high_quality,
    read_length
  )
}


# ==============================================================================
# 15. PRE-FLIGHT TEST SUITE
# ==============================================================================

run_preflight_tests <- function() {

  results <- list()

  # --------------------------------------------------------------------
  # TEST 1
  # DNAString type guard (Fixed: replaced Biostrings::is.* with inherits)
  # --------------------------------------------------------------------

  t1 <- make_synthetic_read(
    center_pos = 653L,
    offset_from_read_start = 10L,
    read_nt_length = 90L
  )

  results[["01_DNAString_type"]] <-
    inherits(t1$seq, "DNAString") &&
    !inherits(t1$seq, "DNAStringSet")


  # --------------------------------------------------------------------
  # TEST 2
  # Forward orientation
  # --------------------------------------------------------------------

  r2 <- locate_position_in_read(
    t1$seq,
    653L,
    tac_ref_chr
  )

  results[["02_forward_orientation"]] <-
    isTRUE(r2$found) &&
    r2$orientation == "fwd" &&
    r2$target_aa ==
      t1$expected_aa


  # --------------------------------------------------------------------
  # TEST 3
  # Reverse-complement orientation
  # --------------------------------------------------------------------

  t3 <- make_synthetic_read(
    center_pos = 698L,
    offset_from_read_start = 10L,
    read_nt_length = 90L,
    reverse = TRUE
  )

  r3 <- locate_position_in_read(
    t3$seq,
    698L,
    tac_ref_chr
  )

  results[["03_reverse_orientation"]] <-
    isTRUE(r3$found) &&
    r3$orientation == "rev" &&
    r3$target_aa ==
      t3$expected_aa


  # --------------------------------------------------------------------
  # TEST 4
  # Frame offset / coordinate recovery
  # --------------------------------------------------------------------

  t4 <- make_synthetic_read(
    center_pos = 653L,
    offset_from_read_start = 10L,
    read_nt_length = 90L
  )

  # Add one nucleotide before the translated region.
  # This forces a non-zero reading frame.
  t4_seq <- Biostrings::DNAString(
    paste0(
      "A",
      as.character(t4$seq)
    )
  )

  r4 <- locate_position_in_read(
    t4_seq,
    653L,
    tac_ref_chr
  )

  results[["04_nonzero_frame"]] <-
    isTRUE(r4$found) &&
    r4$target_aa ==
      t4$expected_aa


  # --------------------------------------------------------------------
  # TEST 5
  # True variant at target position must be detected
  # --------------------------------------------------------------------

  t5 <- make_synthetic_read(
    center_pos = 653L,
    offset_from_read_start = 10L,
    read_nt_length = 90L,
    mutate_pos = 653L,
    mutate_to = "L"
  )

  r5 <- locate_position_in_read(
    t5$seq,
    653L,
    tac_ref_chr
  )

  results[["05_detects_653_variant"]] <-
    isTRUE(r5$found) &&
    r5$target_aa == "L"


  # --------------------------------------------------------------------
  # TEST 6
  # True variant at position 698 must be detected
  # --------------------------------------------------------------------

  t6 <- make_synthetic_read(
    center_pos = 698L,
    offset_from_read_start = 10L,
    read_nt_length = 90L,
    mutate_pos = 698L,
    mutate_to = "P"
  )

  r6 <- locate_position_in_read(
    t6$seq,
    698L,
    tac_ref_chr
  )

  results[["06_detects_698_variant"]] <-
    isTRUE(r6$found) &&
    r6$target_aa == "P"


  # --------------------------------------------------------------------
  # TEST 7
  # Non-covering unrelated read must be rejected
  # --------------------------------------------------------------------

  unrelated_read <- Biostrings::DNAString(
    paste(
      rep(
        "ATGGGTAAACCTTTAGCG",
        5L
      ),
      collapse = ""
    )
  )

  r7 <- locate_position_in_read(
    unrelated_read,
    653L,
    tac_ref_chr
  )

  results[["07_rejects_unrelated_read"]] <-
    !isTRUE(r7$found)


  # --------------------------------------------------------------------
  # TEST 8
  # Forward codon-coordinate recovery
  # --------------------------------------------------------------------

  coord8 <-
    get_original_read_codon_coordinates(
      read_length = length(t1$seq),
      orientation = r2$orientation,
      frame = r2$frame,
      codon_start_in_frame =
        r2$codon_start_in_frame
    )

  results[["08_forward_coordinate_bounds"]] <-
    coord8["start"] >= 1L &&
    coord8["end"] <= length(t1$seq) &&
    coord8["end"] -
      coord8["start"] + 1L == 3L


  # --------------------------------------------------------------------
  # TEST 9
  # Reverse codon-coordinate recovery
  # --------------------------------------------------------------------

  coord9 <-
    get_original_read_codon_coordinates(
      read_length = length(t3$seq),
      orientation = r3$orientation,
      frame = r3$frame,
      codon_start_in_frame =
        r3$codon_start_in_frame
    )

  results[["09_reverse_coordinate_bounds"]] <-
    coord9["start"] >= 1L &&
    coord9["end"] <= length(t3$seq) &&
    coord9["end"] -
      coord9["start"] + 1L == 3L


  # --------------------------------------------------------------------
  # TEST 10
  # Forward recovered nucleotide codon must translate to recovered AA
  # --------------------------------------------------------------------

  codon10 <- as.character(
    Biostrings::subseq(
      t1$seq,
      start = coord8["start"],
      width = 3L
    )
  )

  aa10 <- as.character(
    Biostrings::translate(
      Biostrings::DNAString(codon10),
      genetic.code = GENCODE,
      no.init.codon = TRUE,
      if.fuzzy.codon = "X"
    )
  )

  results[["10_forward_codon_translation"]] <-
    aa10 == r2$target_aa


  # --------------------------------------------------------------------
  # TEST 11
  # Reverse recovered nucleotide codon must translate to recovered AA
  # --------------------------------------------------------------------

  original_rc_codon <- Biostrings::reverseComplement(
    Biostrings::subseq(
      t3$seq,
      start = coord9["start"],
      width = 3L
    )
  )

  aa11 <- as.character(
    Biostrings::translate(
      original_rc_codon,
      genetic.code = GENCODE,
      no.init.codon = TRUE,
      if.fuzzy.codon = "X"
    )
  )

  results[["11_reverse_codon_translation"]] <-
    aa11 == r3$target_aa


  # --------------------------------------------------------------------
  # TEST 12
  # Phred coordinate extraction
  # --------------------------------------------------------------------

  q12 <- make_test_quality(
    length(t1$seq),
    high_quality = 35L
  )

  q_result12 <- extract_codon_quality(
    q12,
    start = coord8["start"],
    end = coord8["end"],
    orientation = r2$orientation
  )

  results[["12_forward_phred_coordinates"]] <-
    length(q_result12$values) == 3L &&
    all(q_result12$values == 35L) &&
    q_result12$min == 35L


  # --------------------------------------------------------------------
  # TEST 13
  # Reverse Phred coordinate extraction
  # --------------------------------------------------------------------

  q13 <- make_test_quality(
    length(t3$seq),
    high_quality = 37L
  )

  q_result13 <- extract_codon_quality(
    q13,
    start = coord9["start"],
    end = coord9["end"],
    orientation = r3$orientation
  )

  results[["13_reverse_phred_coordinates"]] <-
    length(q_result13$values) == 3L &&
    all(q_result13$values == 37L) &&
    q_result13$min == 37L


  # --------------------------------------------------------------------
  # TEST 14
  # Deliberately heterogeneous Phred values
  # --------------------------------------------------------------------

  q14 <- rep(
    40L,
    length(t1$seq)
  )

  q14[coord8["start"]] <- 11L
  q14[coord8["start"] + 1L] <- 23L
  q14[coord8["start"] + 2L] <- 37L

  q_result14 <- extract_codon_quality(
    q14,
    start = coord8["start"],
    end = coord8["end"],
    orientation = "fwd"
  )

  results[["14_exact_phred_window"]] <-
    identical(
      as.integer(q_result14$values),
      c(11L, 23L, 37L)
    ) &&
    q_result14$min == 11L &&
    abs(q_result14$mean - 23.6666667) < 1e-6


  # --------------------------------------------------------------------
  # TEST 15
  # Reverse-complement codon coordinates must correspond to the same
  # biological codon.
  # --------------------------------------------------------------------

  forward_equivalent <- Biostrings::reverseComplement(
    Biostrings::subseq(
      t3$seq,
      start = coord9["start"],
      width = 3L
    )
  )

  reverse_expected <- TEST_CODON_TABLE[
    t3$expected_aa
  ]

  results[["15_reverse_biological_codon"]] <-
    as.character(forward_equivalent) ==
    reverse_expected


  # --------------------------------------------------------------------
  # TEST 16
  # Initiation codon must NOT alter translation.
  # --------------------------------------------------------------------

  initiation_test <- Biostrings::DNAString(
    "ATGCTG"
  )

  initiation_translation <- as.character(
    Biostrings::translate(
      initiation_test,
      genetic.code = GENCODE,
      no.init.codon = TRUE,
      if.fuzzy.codon = "X"
    )
  )

  results[["16_no_init_translation"]] <-
    initiation_translation == "ML"


  # --------------------------------------------------------------------
  # TEST 17
  # DNAString extraction using [[ must produce a single DNAString.
  # (Fixed: replaced Biostrings::is.* with inherits)
  # --------------------------------------------------------------------

  synthetic_set <- Biostrings::DNAStringSet(
    c(
      as.character(t1$seq),
      as.character(t3$seq)
    )
  )

  extracted_single <- synthetic_set[[1L]]

  results[["17_double_bracket_extraction"]] <-
    inherits(extracted_single, "DNAString") &&
    !inherits(extracted_single, "DNAStringSet")


  # --------------------------------------------------------------------
  # TEST 18
  # Reference AA at 653 must be recovered exactly via the standard
  # synthetic-read pathway (mirrors Tests 2/3 but asserts against the
  # TAC102 reference directly rather than the synthetic read's own
  # expected_aa, as an independent cross-check).
  # --------------------------------------------------------------------

  expected_aa_18 <- substr(tac_ref_chr, 653L, 653L)

  t18 <- make_synthetic_read(
    center_pos = 653L,
    offset_from_read_start = 10L,
    read_nt_length = 90L,
    reverse = FALSE
  )

  r18 <- locate_position_in_read(
    t18$seq,
    653L,
    tac_ref_chr
  )

  results[["18_reference_653"]] <-
    isTRUE(r18$found) &&
    identical(
      as.character(r18$target_aa),
      as.character(expected_aa_18)
    )

  # --------------------------------------------------------------------
  # TEST 19
  # target position 698 reference amino acid consistency
  # --------------------------------------------------------------------

  results[["19_reference_698"]] <-
    substr(
      tac_ref_chr,
      698L,
      698L
    ) == "L"


  # --------------------------------------------------------------------
  # TEST 20
  # Synthetic 653 mutation must have expected codon
  # --------------------------------------------------------------------

  results[["20_653_variant_codon"]] <-
    r5$target_codon ==
    unname(TEST_CODON_TABLE["L"])


  # --------------------------------------------------------------------
  # TEST 21
  # Synthetic 698 mutation must have expected codon
  # --------------------------------------------------------------------

  results[["21_698_variant_codon"]] <-
    r6$target_codon ==
    unname(TEST_CODON_TABLE["P"])


  # --------------------------------------------------------------------
  # TEST 22
  # Read length guard
  # --------------------------------------------------------------------

  results[["22_length_guard"]] <-
    length(t1$seq) ==
    nchar(as.character(t1$seq))


  # --------------------------------------------------------------------
  # TEST 23
  # Anchor identity threshold
  # --------------------------------------------------------------------

  t23 <- make_synthetic_read(
    center_pos = 653L,
    offset_from_read_start = 10L,
    read_nt_length = 90L
  )

  t23_pep <- as.character(
    Biostrings::translate(
      t23$seq,
      genetic.code = GENCODE,
      no.init.codon = TRUE,
      if.fuzzy.codon = "X"
    )
  )

  pep_chars23 <- strsplit(
    t23_pep,
    "",
    fixed = TRUE
  )[[1L]]

  target_idx23 <- 10L + 1L

  mutate_idx23 <- target_idx23 - 2L

  if (mutate_idx23 >= 1L &&
      mutate_idx23 <= length(pep_chars23)) {

    pep_chars23[mutate_idx23] <-
      ifelse(
        pep_chars23[mutate_idx23] == "A",
        "G",
        "A"
      )

    mutated_pep23 <- paste(
      pep_chars23,
      collapse = ""
    )

    nt23 <- paste(
      TEST_CODON_TABLE[
        strsplit(
          mutated_pep23,
          "",
          fixed = TRUE
        )[[1L]]
      ],
      collapse = ""
    )

    mutated_read23 <- Biostrings::DNAString(
      nt23
    )

    r23 <- locate_position_in_read(
      mutated_read23,
      653L,
      tac_ref_chr
    )

    results[["23_anchor_tolerates_one_mismatch"]] <-
      isTRUE(r23$found)

  } else {

    results[["23_anchor_tolerates_one_mismatch"]] <-
      FALSE
  }


  # --------------------------------------------------------------------
  # FINAL
  # --------------------------------------------------------------------

  result_df <- data.frame(
    test = names(results),
    passed = unlist(results),
    stringsAsFactors = FALSE
  )

  list(
    all_pass = all(result_df$passed),
    table = result_df
  )
}
# ==============================================================================
# 16. RUN PRE-FLIGHT
# ==============================================================================

cat("\n")
cat("==============================================================\n")
cat("RUNNING 4B.3p PRE-FLIGHT TEST SUITE\n")
cat("==============================================================\n\n")

preflight <- run_preflight_tests()

print(
  preflight$table,
  row.names = FALSE
)

write.csv(
  preflight$table,
  output_preflight,
  row.names = FALSE
)

failed_tests <- preflight$table$test[
  !preflight$table$passed
]

cat("\n")

if (!preflight$all_pass) {

  cat(
    "==============================================================\n"
  )

  cat(
    "PRE-FLIGHT FAILED\n"
  )

  cat(
    "==============================================================\n"
  )

  cat(
    "Failed tests:\n"
  )

  cat(
    paste(
      failed_tests,
      collapse = "\n"
    ),
    "\n\n"
  )

  stop(
    "4B.3p REFUSED TO START FASTQ STREAMING.\n",
    "Fix the failed pre-flight tests before processing the dataset."
  )
}


cat(
  "==============================================================\n"
)

cat(
  "PRE-FLIGHT PASSED: ALL TESTS\n"
)

cat(
  "FASTQ streaming is now permitted.\n"
)

cat(
  "==============================================================\n\n"
)


if (exists("log_info")) {

  log_info(
    "4B.3p pre-flight synthetic tests: ALL PASSED.",
    LOG_FILE
  )
}


# ==============================================================================
# 17. FASTQ INPUT VALIDATION
# ==============================================================================

if (!file.exists(PATHS$fastq_r1)) {

  stop(
    "R1 FASTQ does not exist:\n",
    PATHS$fastq_r1
  )
}

if (!file.exists(PATHS$fastq_r2)) {

  stop(
    "R2 FASTQ does not exist:\n",
    PATHS$fastq_r2
  )
}


# ==============================================================================
# 18. PROCESS ONE FASTQ
# ==============================================================================

process_fastq_for_positions <- function(
    fastq_file,
    mate_label,
    target_positions,
    max_reads = MAX_READS) {

  cat("\n")
  cat(
    "==============================================================\n"
  )

  cat(
    "PROCESSING:",
    mate_label,
    "\n"
  )

  cat(
    "FASTQ:",
    fastq_file,
    "\n"
  )

  if (is.na(max_reads)) {

    cat(
      "Mode: FULL DATASET\n"
    )

  } else {

    cat(
      "Mode: PILOT\n"
    )

    cat(
      "Maximum reads:",
      format(max_reads, big.mark = ","),
      "\n"
    )
  }

  cat(
    "==============================================================\n"
  )


  start_time <- Sys.time()

  fq <- ShortRead::FastqStreamer(
    fastq_file,
    n = STREAM_CHUNK_SIZE
  )

  on.exit(
    close(fq),
    add = TRUE
  )


  result_list <- list()

  chunk_num <- 0L
  total_reads <- 0L
  total_hits <- 0L

  last_report_reads <- 0L


  repeat {

    # --------------------------------------------------------------
    # STOP IF MAX_READS REACHED
    # --------------------------------------------------------------

    if (!is.na(max_reads) &&
        total_reads >= max_reads) {

      break
    }


    chunk_num <- chunk_num + 1L

    reads <- ShortRead::yield(
      fq
    )

    if (length(reads) == 0L) {
      break
    }


    # --------------------------------------------------------------
    # LIMIT LAST CHUNK IN PILOT MODE
    # --------------------------------------------------------------

    if (!is.na(max_reads)) {

      remaining <-
        max_reads -
        total_reads

      if (length(reads) > remaining) {

        reads <- reads[seq_len(remaining)]
      }
    }


    n_chunk <- length(reads)

    if (n_chunk == 0L) {
      break
    }

    total_reads <-
      total_reads +
      n_chunk


    # --------------------------------------------------------------
    # READ SEQUENCES
    # --------------------------------------------------------------

    seqs <- ShortRead::sread(
      reads
    )

    ids <- sub(
      "\\s.*$",
      "",
      as.character(
        ShortRead::id(reads)
      )
    )


    # --------------------------------------------------------------
    # QUALITY MATRIX
    # --------------------------------------------------------------

    quals <- slot(
      reads,
      "quality"
    )

    qual_matrix <- as(quals, "matrix")


    # --------------------------------------------------------------
    # PROCESS EACH READ
    # --------------------------------------------------------------

    for (i in seq_len(n_chunk)) {

      # ------------------------------------------------------------
      # IMPORTANT:
      #
      # [[ extracts ONE DNAString.
      #
      # This is intentional and protects against the earlier
      # DNAString/DNAStringSet length bug.
      # ------------------------------------------------------------

      read_seq <- seqs[[i]]

            if (!inherits(read_seq, "DNAString")) {

        stop(
          "Internal type failure: seqs[[i]] is not DNAString."
        )
      }


      read_length <- length(
        read_seq
      )


      read_position_hits <- list()


      for (pos in target_positions) {

        hit <- locate_position_in_read(
          read_seq = read_seq,
          target_pos = pos,
          ref_chr = tac_ref_chr,
          flank = FLANK_WINDOW,
          min_identity = MIN_ANCHOR_IDENTITY
        )


        if (!isTRUE(hit$found)) {
          next
        }


        # ----------------------------------------------------------
        # ORIGINAL READ CODON COORDINATES
        # ----------------------------------------------------------

        coords <-
          get_original_read_codon_coordinates(
            read_length = read_length,
            orientation = hit$orientation,
            frame = hit$frame,
            codon_start_in_frame =
              hit$codon_start_in_frame
          )


        read_pos_start <-
          as.integer(coords["start"])

        read_pos_end <-
          as.integer(coords["end"])


        # ----------------------------------------------------------
        # HARD COORDINATE GUARD
        # ----------------------------------------------------------

        if (read_pos_start < 1L ||
            read_pos_end > ncol(qual_matrix) ||
            read_pos_end -
            read_pos_start +
            1L != 3L) {

          stop(
            "Invalid recovered codon coordinates.\n",
            "Read: ",
            ids[i],
            "\n",
            "Mate: ",
            mate_label,
            "\n",
            "Position: ",
            pos,
            "\n",
            "Coordinates: ",
            read_pos_start,
            "-",
            read_pos_end,
            "\n",
            "Read length: ",
            read_length
          )
        }


        # ----------------------------------------------------------
        # QUALITY
        # ----------------------------------------------------------

        q_result <-
          extract_codon_quality(
            quality_vector =
              qual_matrix[i, ],
            start =
              read_pos_start,
            end =
              read_pos_end,
            orientation =
              hit$orientation
          )


        # ----------------------------------------------------------
        # INDEPENDENT CODON RECONSTRUCTION
        #
        # This is an additional runtime guard:
        #
        # recovered nucleotide codon
        #         ↓
        # reverse-complement if required
        #         ↓
        # translate
        #         ↓
        # recovered amino acid
        #
        # must equal the amino acid reported by the locator.
        # ----------------------------------------------------------

        original_codon <- as.character(
          Biostrings::subseq(
            read_seq,
            start = read_pos_start,
            width = 3L
          )
        )


        biological_codon <- original_codon

        if (hit$orientation == "rev") {

          biological_codon <-
            as.character(
              Biostrings::reverseComplement(
                Biostrings::DNAString(
                  original_codon
                )
              )
            )
        }


        reconstructed_aa <- as.character(
          Biostrings::translate(
            Biostrings::DNAString(
              biological_codon
            ),
            genetic.code = GENCODE,
            no.init.codon = TRUE,
            if.fuzzy.codon = "X"
          )
        )


        if (reconstructed_aa !=
            hit$target_aa) {

          stop(
            "RUNTIME CODON/AA CONSISTENCY FAILURE.\n",
            "Read: ",
            ids[i],
            "\n",
            "Mate: ",
            mate_label,
            "\n",
            "Position: ",
            pos,
            "\n",
            "Reported AA: ",
            hit$target_aa,
            "\n",
            "Reconstructed AA: ",
            reconstructed_aa,
            "\n",
            "Codon: ",
            biological_codon
          )
        }


        # ----------------------------------------------------------
        # STORE HIT
        # ----------------------------------------------------------

        hit_record <- data.frame(

          read_id = ids[i],

          mate = mate_label,

          target_position = pos,

          orientation = hit$orientation,

          frame = hit$frame,

          anchor_identity =
            round(hit$score, 4),

          anchor_reference_start =
            hit$anchor_start_reference,

          anchor_reference_end =
            hit$anchor_end_reference,

          anchor_peptide =
            hit$anchor_peptide,

          anchor_candidate =
            hit$anchor_candidate,

          target_aa =
            hit$target_aa,

          target_codon =
            biological_codon,

          read_codon_start =
            read_pos_start,

          read_codon_end =
            read_pos_end,

          codon_min_qual =
            q_result$min,

          codon_mean_qual =
            round(q_result$mean, 2),

          read_length =
            read_length,

          stringsAsFactors = FALSE
        )


        read_position_hits[[length(read_position_hits) + 1L]] <- hit_record
      }


      if (length(read_position_hits) > 0L) {

        result_list[[length(result_list) + 1L]] <- do.call(
          rbind,
          read_position_hits
        )

        total_hits <-
          total_hits +
          length(read_position_hits)
      }
    }


    # --------------------------------------------------------------
    # PROGRESS
    # --------------------------------------------------------------

    if (
      total_reads - last_report_reads >=
      PROGRESS_INTERVAL
    ) {

      elapsed <-
        as.numeric(
          difftime(
            Sys.time(),
            start_time,
            units = "secs"
          )
        )

      rate <-
        ifelse(
          elapsed > 0,
          total_reads / elapsed,
          NA_real_
        )


      cat(
        sprintf(
          "  %s | reads=%s | hits=%s | rate=%.1f reads/sec | elapsed=%.1f min\n",
          mate_label,
          format(total_reads, big.mark = ","),
          format(total_hits, big.mark = ","),
          rate,
          elapsed / 60
        )
      )


      if (exists("log_info")) {

        log_info(
          sprintf(
            "4B.3p %s: %s reads processed; %s hits.",
            mate_label,
            format(total_reads, big.mark = ","),
            format(total_hits, big.mark = ",")
          ),
          LOG_FILE
        )
      }


      last_report_reads <-
        total_reads
    }
  }


  # ------------------------------------------------------------------
  # FINAL TIMING
  # ------------------------------------------------------------------

  elapsed <-
    as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "secs"
      )
    )


  rate <-
    ifelse(
      elapsed > 0,
      total_reads / elapsed,
      NA_real_
    )


  cat(
    "\n",
    mate_label,
    "COMPLETE\n",
    sep = ""
  )

  cat(
    "Reads processed:",
    format(total_reads, big.mark = ","),
    "\n"
  )

  cat(
    "Position hits:",
    format(total_hits, big.mark = ","),
    "\n"
  )

  cat(
    "Elapsed:",
    round(elapsed / 60, 2),
    "minutes\n"
  )

  cat(
    "Read rate:",
    round(rate, 2),
    "reads/sec\n"
  )


  if (length(result_list) == 0L) {

    hits_df <- data.frame(

      read_id = character(),

      mate = character(),

      target_position = integer(),

      orientation = character(),

      frame = integer(),

      anchor_identity = numeric(),

      anchor_reference_start = integer(),

      anchor_reference_end = integer(),

      anchor_peptide = character(),

      anchor_candidate = character(),

      target_aa = character(),

      target_codon = character(),

      read_codon_start = integer(),

      read_codon_end = integer(),

      codon_min_qual = numeric(),

      codon_mean_qual = numeric(),

      read_length = integer(),

      stringsAsFactors = FALSE
    )

  } else {

    hits_df <- do.call(
      rbind,
      result_list
    )
  }


  list(

    hits = hits_df,

    total_reads = total_reads,

    total_hits = total_hits,

    elapsed_seconds = elapsed,

    reads_per_second = rate
  )
}


# ==============================================================================
# 19. RUN R1
# ==============================================================================

cat(
  "\nStarting R1 positional scan...\n"
)

r1_all <- process_fastq_for_positions(
  fastq_file = PATHS$fastq_r1,
  mate_label = "R1",
  target_positions = TARGET_POSITIONS,
  max_reads = MAX_READS
)


# ==============================================================================
# 20. RUN R2
# ==============================================================================

cat(
  "\nStarting R2 positional scan...\n"
)

r2_all <- process_fastq_for_positions(
  fastq_file = PATHS$fastq_r2,
  mate_label = "R2",
  target_positions = TARGET_POSITIONS,
  max_reads = MAX_READS
)


# ==============================================================================
# 21. COMBINE RESULTS
# ==============================================================================

all_hits <- rbind(
  r1_all$hits,
  r2_all$hits
)


# ==============================================================================
# 22. REMOVE IMPOSSIBLE DUPLICATE POSITION RECORDS
# ==============================================================================
#
# A read should contribute at most one confidently located instance of each
# target position.
#
# If the locator somehow returns duplicated records for the same read/mate/
# position, this is treated as an internal methodological warning rather than
# silently accepted.
# ==============================================================================

duplicate_key <- paste(
  all_hits$read_id,
  all_hits$mate,
  all_hits$target_position,
  sep = "|"
)

duplicate_counts <- table(
  duplicate_key
)

duplicate_keys <- names(
  duplicate_counts[
    duplicate_counts > 1L
  ]
)

if (length(duplicate_keys) > 0L) {

  warning(
    "Multiple hits were generated for the same read/mate/position.\n",
    "Inspect duplicate_key records before interpreting results."
  )
}


# ==============================================================================
# 23. SAVE RAW POSITIONAL EVIDENCE
# ==============================================================================

if (SAVE_ALL_HITS) {

  write.csv(
    all_hits,
    output_hits,
    row.names = FALSE
  )
}


# ==============================================================================
# 24. POSITION SUMMARY
# ==============================================================================

make_position_summary <- function(
    hits,
    position) {

  sub <- hits[
    hits$target_position == position,
    ,
    drop = FALSE
  ]


  if (nrow(sub) == 0L) {

    return(
      data.frame(

        target_position = position,

        total_hits = 0L,

        unique_reads = 0L,

        M = 0L,

        L = 0L,

        P = 0L,

        other = 0L,

        pct_reference_aa = NA_real_,

        mean_anchor_identity = NA_real_,

        median_anchor_identity = NA_real_,

        mean_codon_min_qual = NA_real_,

        median_codon_min_qual = NA_real_,

        stringsAsFactors = FALSE
      )
    )
  }


  ref_aa <- substr(
    tac_ref_chr,
    position,
    position
  )


  total <- nrow(sub)

  ref_count <-
    sum(sub$target_aa == ref_aa)

  M_count <-
    sum(sub$target_aa == "M")

  L_count <-
    sum(sub$target_aa == "L")

  P_count <-
    sum(sub$target_aa == "P")

  other_count <-
    sum(
      !(sub$target_aa %in% c(
        "M",
        "L",
        "P"
      ))
    )


  data.frame(

    target_position = position,

    total_hits = total,

    unique_reads =
      length(unique(sub$read_id)),

    M = M_count,

    L = L_count,

    P = P_count,

    other = other_count,

    pct_reference_aa =
      round(
        100 *
          ref_count /
          total,
        2
      ),

    mean_anchor_identity =
      round(
        mean(sub$anchor_identity),
        4
      ),

    median_anchor_identity =
      round(
        median(sub$anchor_identity),
        4
      ),

    mean_codon_min_qual =
      round(
        mean(sub$codon_min_qual),
        2
      ),

    median_codon_min_qual =
      round(
        median(sub$codon_min_qual),
        2
      ),

    stringsAsFactors = FALSE
  )
}


position_summary <- do.call(
  rbind,
  lapply(
    TARGET_POSITIONS,
    function(x) {
      make_position_summary(
        all_hits,
        x
      )
    }
  )
)


# ==============================================================================
# 25. DETAILED ALLELE COUNTS
# ==============================================================================

allele_summary <- aggregate(
  read_id ~ mate + target_position + target_aa + target_codon,
  data = all_hits,
  FUN = length
)

names(allele_summary)[
  names(allele_summary) == "read_id"
] <- "hit_count"


if (nrow(allele_summary) > 0L) {

  allele_summary$percent <- NA_real_

  for (i in seq_len(nrow(allele_summary))) {

    denom <- sum(
      allele_summary$hit_count[
        allele_summary$mate ==
          allele_summary$mate[i] &
        allele_summary$target_position ==
          allele_summary$target_position[i]
      ]
    )

    allele_summary$percent[i] <-
      round(
        100 *
          allele_summary$hit_count[i] /
          denom,
        2
      )
  }
}


# ==============================================================================
# 26. QUALITY-STRATIFIED ALL-READ SUMMARY
# ==============================================================================

quality_thresholds <- c(
  0L,
  10L,
  15L,
  20L,
  25L,
  30L
)

quality_summary <- list()

for (q in quality_thresholds) {

  sub <- all_hits[
    all_hits$codon_min_qual >= q,
    ,
    drop = FALSE
  ]

  if (nrow(sub) == 0L) {

    quality_summary[[length(quality_summary) + 1L]] <- data.frame(

      min_phred = q,

      total_hits = 0L,

      position_653_hits = 0L,

      position_653_M = 0L,

      position_653_non_M = 0L,

      position_698_hits = 0L,

      position_698_L = 0L,

      position_698_P = 0L,

      position_698_other = 0L,

      stringsAsFactors = FALSE
    )

    next
  }


  sub653 <- sub[
    sub$target_position == 653L,
    ,
    drop = FALSE
  ]

  sub698 <- sub[
    sub$target_position == 698L,
    ,
    drop = FALSE
  ]


  quality_summary[[length(quality_summary) + 1L]] <- data.frame(

    min_phred = q,

    total_hits = nrow(sub),

    position_653_hits =
      nrow(sub653),

    position_653_M =
      sum(sub653$target_aa == "M"),

    position_653_non_M =
      sum(sub653$target_aa != "M"),

    position_698_hits =
      nrow(sub698),

    position_698_L =
      sum(sub698$target_aa == "L"),

    position_698_P =
      sum(sub698$target_aa == "P"),

    position_698_other =
      sum(
        !(sub698$target_aa %in% c(
          "L",
          "P"
        ))
      ),

    stringsAsFactors = FALSE
  )
}

quality_summary <- do.call(
  rbind,
  quality_summary
)


# ==============================================================================
# 27. HAPLOTYPE ANALYSIS AMONG ALL-READ DATA
# ==============================================================================
#
# Important:
#
# We can only phase 653 and 698 when the SAME read covers both positions.
#
# This is not paired-read phasing yet.
#
# This is single-read physical linkage.
#
# A future extension could identify R1/R2 pairs, but this section deliberately
# stays within the independent positional validation mechanism.
# ==============================================================================

hits_by_read <- split(
  all_hits,
  paste(
    all_hits$mate,
    all_hits$read_id,
    sep = "|"
  )
)

haplotype_records <- list()

for (key in names(hits_by_read)) {

  sub <- hits_by_read[[key]]

  has653 <-
    any(
      sub$target_position == 653L
    )

  has698 <-
    any(
      sub$target_position == 698L
    )

  if (!has653 || !has698) {
    next
  }


  row653 <- sub[
    sub$target_position == 653L,
    ,
    drop = FALSE
  ]

  row698 <- sub[
    sub$target_position == 698L,
    ,
    drop = FALSE
  ]


  # If multiple records somehow exist, do not silently choose one.
  if (nrow(row653) != 1L ||
      nrow(row698) != 1L) {

    next
  }


  haplotype_records[[length(haplotype_records) + 1L]] <- data.frame(

    read_id = row653$read_id,

    mate = row653$mate,

    aa_653 = row653$target_aa,

    codon_653 = row653$target_codon,

    qmin_653 = row653$codon_min_qual,

    aa_698 = row698$target_aa,

    codon_698 = row698$target_codon,

    qmin_698 = row698$codon_min_qual,

    min_phred_both =
      min(
        row653$codon_min_qual,
        row698$codon_min_qual
      ),

    haplotype_aa =
      paste(
        row653$target_aa,
        row698$target_aa,
        sep = "653-"
      ),

    haplotype_codon =
      paste(
        row653$target_codon,
        row698$target_codon,
        sep = "-"
      ),

    stringsAsFactors = FALSE
  )
}


if (length(haplotype_records) > 0L) {

  haplotype_df <- do.call(
    rbind,
    haplotype_records
  )

} else {

  haplotype_df <- data.frame(

    read_id = character(),

    mate = character(),

    aa_653 = character(),

    codon_653 = character(),

    qmin_653 = numeric(),

    aa_698 = character(),

    codon_698 = character(),

    qmin_698 = numeric(),

    min_phred_both = numeric(),

    haplotype_aa = character(),

    haplotype_codon = character(),

    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 28. CORRECT HAPLOTYPE LABEL
# ==============================================================================

if (nrow(haplotype_df) > 0L) {

  haplotype_df$haplotype_aa <-
    paste0(
      haplotype_df$aa_653,
      "653-",
      haplotype_df$aa_698,
      "698"
    )
}


# ==============================================================================
# 29. HAPLOTYPE SUMMARY
# ==============================================================================

if (nrow(haplotype_df) > 0L) {

  haplotype_summary <- aggregate(
    read_id ~ haplotype_aa + haplotype_codon,
    data = haplotype_df,
    FUN = length
  )

  names(haplotype_summary)[
    names(haplotype_summary) == "read_id"
  ] <- "read_count"

  haplotype_summary$percent <-
    round(
      100 *
        haplotype_summary$read_count /
        sum(haplotype_summary$read_count),
      2
    )

} else {

  haplotype_summary <- data.frame(

    haplotype_aa = character(),

    haplotype_codon = character(),

    read_count = integer(),

    percent = numeric(),

    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 30. SAVE SUMMARIES
# ==============================================================================

write.csv(
  position_summary,
  output_summary,
  row.names = FALSE
)

write.csv(
  haplotype_summary,
  output_haplotype,
  row.names = FALSE
)


# ==============================================================================
# 31. WRITE TEXT REPORT
# ==============================================================================

run_mode <- ifelse(
  is.na(MAX_READS),
  "FULL DATASET",
  paste0(
    "PILOT (MAX_READS = ",
    format(MAX_READS, big.mark = ","),
    " per FASTQ)"
  )
)


report_lines <- c(

  "==============================================================",

  "TAC102 ALL-READ POSITIONAL VALIDATION",

  "Step 4B.3p",

  "==============================================================",

  "",

  paste(
    "Date:",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  paste(
    "Run mode:",
    run_mode
  ),

  "",

  "DESIGN",

  "This analysis independently examines the original FASTQ files",

  "for reads capable of covering TAC102 protein positions 653 and 698.",

  "",

  "It does not use the previous BLAST-HSP candidate-read selection",

  "pipeline or any 4B.3j-4B.3o output.",

  "",

  "Each read is examined in six possible translation frames",

  "using direct peptide-anchor matching to the TAC102 reference.",

  "",

  paste(
    "Flank window:",
    FLANK_WINDOW,
    "aa"
  ),

  paste(
    "Minimum anchor identity:",
    MIN_ANCHOR_IDENTITY
  ),

  "",

  "PRE-FLIGHT",

  paste(
    "All pre-flight tests passed:",
    preflight$all_pass
  ),

  "",

  "READ PROCESSING",

  paste(
    "R1 reads processed:",
    format(
      r1_all$total_reads,
      big.mark = ","
    )
  ),

  paste(
    "R2 reads processed:",
    format(
      r2_all$total_reads,
      big.mark = ","
    )
  ),

  paste(
    "R1 positional hits:",
    format(
      r1_all$total_hits,
      big.mark = ","
    )
  ),

  paste(
    "R2 positional hits:",
    format(
      r2_all$total_hits,
      big.mark = ","
    )
  ),

  "",

  "POSITION SUMMARY",

  capture.output(
    print(
      position_summary,
      row.names = FALSE
    )
  ),

  "",

  "ALLELE SUMMARY",

  if (nrow(allele_summary) > 0L)
    capture.output(
      print(
        allele_summary,
        row.names = FALSE
      )
    )
  else
    "No positional hits detected.",

  "",

  "QUALITY-STRATIFIED SUMMARY",

  capture.output(
    print(
      quality_summary,
      row.names = FALSE
    )
  ),

  "",

  "SINGLE-READ HAPLOTYPE SUMMARY",

  if (nrow(haplotype_summary) > 0L)
    capture.output(
      print(
        haplotype_summary,
        row.names = FALSE
      )
    )
  else
    "No reads independently covering both positions were detected.",

  "",

  "INTERPRETATION CAUTION",

  "This analysis identifies sequence calls among reads that can be",

  "independently anchored to the TAC102 reference around the target",

  "positions. It does not by itself establish genomic fixation,",

  "population allele frequency, or biological phenotype.",

  "",

  "High-quality support should be interpreted together with",

  "coverage depth, strand/orientation consistency, independent",

  "read-start families, paired-read evidence, and independent",

  "genomic confirmation.",

  "",

  "OUTPUTS",

  paste(
    "Positional evidence:",
    output_hits
  ),

  paste(
    "Position summary:",
    output_summary
  ),

  paste(
    "All-read haplotype summary:",
    output_haplotype
  ),

  paste(
    "Pre-flight results:",
    output_preflight
  ),

  paste(
    "Report:",
    output_report
  ),

  "",

  "=============================================================="
)


writeLines(
  report_lines,
  output_report
)


# ==============================================================================
# 32. FINAL CONSOLE SUMMARY
# ==============================================================================

cat("\n")
cat("==============================================================\n")
cat("4B.3p COMPLETE\n")
cat("==============================================================\n\n")

cat(
  "Run mode:",
  run_mode,
  "\n"
)

cat(
  "R1 reads:",
  format(r1_all$total_reads, big.mark = ","),
  "\n"
)

cat(
  "R2 reads:",
  format(r2_all$total_reads, big.mark = ","),
  "\n"
)

cat(
  "Total positional hits:",
  format(nrow(all_hits), big.mark = ","),
  "\n\n"
)

cat(
  "POSITION SUMMARY\n"
)

print(
  position_summary,
  row.names = FALSE
)

cat(
  "\nALL-READ HAPLOTYPES\n"
)

if (nrow(haplotype_summary) > 0L) {

  print(
    haplotype_summary,
    row.names = FALSE
  )

} else {

  cat(
    "No single reads covering both positions detected.\n"
  )
}


cat("\n")
cat(
  "Outputs written to:\n"
)

cat(
  output_hits,
  "\n"
)

cat(
  output_summary,
  "\n"
)

cat(
  output_haplotype,
  "\n"
)

cat(
  output_report,
  "\n"
)

cat("\n")
cat(
  "==============================================================\n"
)


if (exists("log_info")) {

  log_info(
    sprintf(
      paste0(
        "4B.3p complete. R1=%s reads; R2=%s reads; ",
        "total positional hits=%s."
      ),
      format(r1_all$total_reads, big.mark = ","),
      format(r2_all$total_reads, big.mark = ","),
      format(nrow(all_hits), big.mark = ",")
    ),
    LOG_FILE
  )
}