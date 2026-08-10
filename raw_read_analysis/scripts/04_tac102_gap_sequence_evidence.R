# ==============================================================================
# Step 4B.3b: TAC102 Gap-Facing Sequence Evidence
#
# Project:
#   Genomic & Transcriptomic Limits on kDNA in T. equiperdum
#
# Script:
#   raw_read_analysis/scripts/04_tac102_gap_sequence_evidence.R
#
# PURPOSE
# -------
# Recover the full nucleotide sequences of the 41 TAC102 both-flank read pairs
# identified in Step 4B.3a and examine the gap-facing portions of those reads.
#
# This is an R/Bioconductor-only analysis.
#
# NO:
#   - Conda
#   - WSL
#   - seqtk
#   - Python
#   - external executables
#
# FASTQ INPUTS
# ------------
#   data/reads/SRR7910035_1.fastq
#   data/reads/SRR7910035_2.fastq
#
# PRE-REGISTERED THRESHOLDS
# -------------------------
# Minimum candidate nucleotide length:
#   30 nt
#
# Class A:
#   >=10 aa matched
#   >=50% amino-acid identity
#   original-HSP-frame, codon-anchored comparison
#
# Class B:
#   >=6 aa matched
#   >=30% amino-acid identity
#   fails Class A
#
# Class C:
#   >=30 nt interpretable candidate
#   does not satisfy A or B
#
# Class D:
#   <30 nt or technically uninterpretable
#
# Six-frame translation and composition similarity are descriptive only.
# They NEVER determine A/B/C/D classification.
#
# IMPORTANT PROTOCOL DEVIATION
# ----------------------------
# The original frozen Class A protocol referred to a nucleotide-level
# comparison against the expected TAC102 gap sequence.
#
# TAC102 in this independent raw-read pipeline has a protein reference only.
# No TAC102 nucleotide reference is introduced here.
#
# Therefore the nucleotide criterion is explicitly substituted by:
#
#   original-HSP-frame, codon-anchored nucleotide extraction followed by
#   direct positional amino-acid comparison against TAC102 reference
#   positions 649-732.
#
# This is flagged in all relevant outputs.
#
# ==============================================================================


cat(
  "==============================================================\n",
  " Step 4B.3b: TAC102 Gap-Facing Sequence Evidence            \n",
  "==============================================================\n\n",
  sep = ""
)


# ==============================================================================
# 0. INITIALISE PROJECT
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


log_info(
  "Starting Step 4B.3b: TAC102 gap-facing sequence evidence.",
  LOG_FILE
)


# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

MINIMUM_CANDIDATE_NT <- 30L

# TAC102 reference gap established in previous analysis.
GAP_START <- 649L
GAP_END   <- 732L

# N-terminal flank immediately preceding the gap.
N_FLANK_MIN <- 550L
N_FLANK_MAX <- 648L

# C-terminal flank immediately following the gap.
C_FLANK_MIN <- 733L
C_FLANK_MAX <- 830L


# Raw FASTQ files.
#
# IMPORTANT:
# These are SOURCE files, not derivative files.
#
R1_FASTQ <- here::here(
  "data",
  "reads",
  "SRR7910035_1.fastq"
)

R2_FASTQ <- here::here(
  "data",
  "reads",
  "SRR7910035_2.fastq"
)


# Output directory.
DERIVATIVE_DIR <- PATHS$raw_recovery_derivative

if (!dir.exists(DERIVATIVE_DIR)) {
  dir.create(
    DERIVATIVE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 2. CHECK REQUIRED INPUTS
# ==============================================================================

cat("Checking input files...\n\n")

if (!file.exists(R1_FASTQ)) {
  stop(
    paste0(
      "ERROR: R1 FASTQ was not found:\n",
      R1_FASTQ,
      "\n\nExpected location:\n",
      here::here("data", "reads", "SRR7910035_1.fastq")
    )
  )
}

if (!file.exists(R2_FASTQ)) {
  stop(
    paste0(
      "ERROR: R2 FASTQ was not found:\n",
      R2_FASTQ,
      "\n\nExpected location:\n",
      here::here("data", "reads", "SRR7910035_2.fastq")
    )
  )
}

cat("R1 FASTQ: FOUND\n")
cat("R2 FASTQ: FOUND\n\n")


# ==============================================================================
# 3. LOAD REFERENCE AND PREVIOUS HSP RESULTS
# ==============================================================================

cat("Loading TAC102 reference and previous HSP results...\n\n")


reference_sequences <- readRDS(
  file.path(
    PATHS$intermediate,
    "reference_sequences.rds"
  )
)

if (!"TAC102" %in% names(reference_sequences)) {
  stop(
    "TAC102 was not found in reference_sequences.rds."
  )
}

tac_ref <- as.character(
  reference_sequences[["TAC102"]]
)

tac_ref_length <- nchar(tac_ref)


per_read_file <- file.path(
  PATHS$reports,
  "raw_read_recovery_per_read.csv"
)

all_hsps_file <- file.path(
  PATHS$reports,
  "raw_read_recovery_all_hsps.csv"
)


if (!file.exists(per_read_file)) {
  stop(
    paste0(
      "Missing required file:\n",
      per_read_file
    )
  )
}

if (!file.exists(all_hsps_file)) {
  stop(
    paste0(
      "Missing required file:\n",
      all_hsps_file
    )
  )
}


per_read_df <- read.csv(
  per_read_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all_hsps <- read.csv(
  all_hsps_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


tac_pr <- per_read_df[
  per_read_df$gene == "TAC102",
  ,
  drop = FALSE
]

tac_hsps <- all_hsps[
  all_hsps$gene == "TAC102",
  ,
  drop = FALSE
]


cat(
  sprintf(
    "TAC102 reference length: %d aa\n",
    tac_ref_length
  )
)

cat(
  sprintf(
    "TAC102 per-read rows: %d\n",
    nrow(tac_pr)
  )
)

cat(
  sprintf(
    "TAC102 raw HSP rows: %d\n\n",
    nrow(tac_hsps)
  )
)


# ==============================================================================
# 4. CHECK REQUIRED HSP COLUMNS
# ==============================================================================

required_hsp_columns <- c(
  "sseqid",
  "mate",
  "qstart",
  "qend",
  "sstart",
  "send",
  "qstart_raw",
  "qend_raw",
  "sstart_raw",
  "send_raw",
  "bitscore"
)


missing_hsp_columns <- setdiff(
  required_hsp_columns,
  names(tac_hsps)
)

if (length(missing_hsp_columns) > 0) {

  stop(
    paste0(
      "Missing required HSP columns:\n",
      paste(
        missing_hsp_columns,
        collapse = ", "
      )
    )
  )
}


required_per_read_columns <- c(
  "gene",
  "mate",
  "read_id",
  "pair_id",
  "qstart_union",
  "qend_union"
)


missing_per_read_columns <- setdiff(
  required_per_read_columns,
  names(tac_pr)
)

if (length(missing_per_read_columns) > 0) {

  stop(
    paste0(
      "Missing required per-read columns:\n",
      paste(
        missing_per_read_columns,
        collapse = ", "
      )
    )
  )
}


# Ensure numeric HSP fields are numeric.
numeric_hsp_columns <- c(
  "qstart",
  "qend",
  "sstart",
  "send",
  "qstart_raw",
  "qend_raw",
  "sstart_raw",
  "send_raw",
  "bitscore"
)


for (nm in numeric_hsp_columns) {
  tac_hsps[[nm]] <- as.numeric(
    tac_hsps[[nm]]
  )
}


# ==============================================================================
# 5. RECONSTRUCT THE 41 BOTH-FLANK PAIRS
#
# IMPORTANT:
# Mate role is determined dynamically.
#
# We do NOT assume:
#
#   R1 = N
#   R2 = C
#
# because Step 4B.3a demonstrated that this assumption fails for at least
# two pairs.
# ==============================================================================


cat(
  "Reconstructing both-flank pairs with dynamic N/C assignment...\n\n"
)


tac_pr$qstart_union <- as.character(
  tac_pr$qstart_union
)

tac_pr$qend_union <- as.character(
  tac_pr$qend_union
)


parse_coordinate_vector <- function(x) {

  if (is.na(x) || !nzchar(x)) {
    return(integer(0))
  }

  pieces <- strsplit(
    x,
    ";",
    fixed = TRUE
  )[[1]]

  pieces <- pieces[
    nzchar(pieces)
  ]

  as.integer(pieces)
}


pair_ids <- unique(
  tac_pr$pair_id
)


pair_role_rows <- vector(
  "list",
  length(pair_ids)
)


for (i in seq_along(pair_ids)) {

  pid <- pair_ids[i]

  sub_df <- tac_pr[
    tac_pr$pair_id == pid,
    ,
    drop = FALSE
  ]

  role_n <- NA_character_
  role_c <- NA_character_


  for (j in seq_len(nrow(sub_df))) {

    starts <- parse_coordinate_vector(
      sub_df$qstart_union[j]
    )

    ends <- parse_coordinate_vector(
      sub_df$qend_union[j]
    )


    if (length(starts) > 0 &&
        length(ends) > 0) {

      n_overlap <- any(
        starts <= N_FLANK_MAX &
        ends >= N_FLANK_MIN
      )

      c_overlap <- any(
        starts <= C_FLANK_MAX &
        ends >= C_FLANK_MIN
      )


      if (n_overlap) {
        role_n <- sub_df$mate[j]
      }

      if (c_overlap) {
        role_c <- sub_df$mate[j]
      }
    }
  }


  pair_role_rows[[i]] <- data.frame(
    pair_id = pid,
    n_mate = role_n,
    c_mate = role_c,
    stringsAsFactors = FALSE
  )
}


pair_roles <- do.call(
  rbind,
  pair_role_rows
)


both_flank <- pair_roles[
  !is.na(pair_roles$n_mate) &
  !is.na(pair_roles$c_mate) &
  pair_roles$n_mate != pair_roles$c_mate,
  ,
  drop = FALSE
]


cat(
  sprintf(
    "Both-flank pairs identified: %d\n\n",
    nrow(both_flank)
  )
)


if (nrow(both_flank) != 41) {

  warning(
    paste0(
      "Expected 41 both-flank pairs based on Step 4B.3a, ",
      "but reconstructed ",
      nrow(both_flank),
      "."
    )
  )
}


# Save dynamic role table.
write.csv(
  both_flank,
  file.path(
    DERIVATIVE_DIR,
    "tac102_41pairs_dynamic_roles.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 6. DETERMINE THE FLANK-RELEVANT RAW HSP
#
# For each mate we choose the highest-bitscore HSP overlapping the appropriate
# flank.
#
# This uses the actual raw HSP rows, not the union coordinates.
# ==============================================================================


get_flank_hsp <- function(
    read_id,
    mate,
    side
) {

  sub <- tac_hsps[
    tac_hsps$sseqid == read_id &
    tac_hsps$mate == mate,
    ,
    drop = FALSE
  ]


  if (nrow(sub) == 0) {
    return(NULL)
  }


  if (side == "N") {

    keep <- (
      sub$qstart <= N_FLANK_MAX &
      sub$qend >= N_FLANK_MIN
    )

  } else {

    keep <- (
      sub$qstart <= C_FLANK_MAX &
      sub$qend >= C_FLANK_MIN
    )
  }


  sub <- sub[
    keep,
    ,
    drop = FALSE
  ]


  if (nrow(sub) == 0) {
    return(NULL)
  }


  sub[
    which.max(sub$bitscore),
    ,
    drop = FALSE
  ]
}


# ==============================================================================
# 7. FASTQ READ EXTRACTION
#
# Entirely R/Bioconductor.
#
# We stream each FASTQ once and retain only the requested read IDs.
#
# The original FASTQ contains ~28.3 million reads per mate, so the entire
# FASTQ is NEVER loaded into memory.
# ==============================================================================


cat(
  "Preparing FASTQ read extraction...\n\n"
)


# We need one nucleotide sequence for each member of each pair.
#
# pair_id is the read identifier in the current dataset.
#
# For every both-flank pair, both R1 and R2 are requested because the N/C role
# may be either R1/R2.


target_ids <- unique(
  both_flank$pair_id
)


cat(
  sprintf(
    "Target pair/read IDs: %d\n",
    length(target_ids)
  )
)


extract_fastq_reads <- function(
    fastq_file,
    target_ids,
    mate_label,
    output_fasta
) {

  cat(
    sprintf(
      "\nStreaming %s FASTQ:\n%s\n",
      mate_label,
      fastq_file
    )
  )


  target_ids <- unique(
    as.character(target_ids)
  )


  found <- character(0)

  extracted_sequences <- DNAStringSet()

  extracted_names <- character(0)

  malformed_chunks <- integer(0)

  chunk_number <- 0L


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


  repeat {

    chunk_number <- chunk_number + 1L

    reads <- tryCatch(
      ShortRead::yield(fq),
      error = function(e) {

        stop(
          paste0(
            "\nFASTQ parsing failed while streaming ",
            mate_label,
            ".\n\n",
            "File: ",
            fastq_file,
            "\n",
            "Chunk: ",
            chunk_number,
            "\n\n",
            "Original error:\n",
            conditionMessage(e),
            "\n\n",
            "No read sequences from this failed chunk were used."
          ),
          call. = FALSE
        )
      }
    )


    if (length(reads) == 0) {
      break
    }


    read_names <- as.character(
      ShortRead::id(reads)
    )


    # Remove common FASTQ read-name suffixes.
    #
    # This preserves the core SRR7910035.<read_id> identifier.
    read_names_clean <- sub(
      "\\s.*$",
      "",
      read_names
    )

    read_names_clean <- sub(
      "/[12]$",
      "",
      read_names_clean
    )


    keep <- read_names_clean %in% target_ids


    if (any(keep)) {

      selected <- reads[keep]

      selected_names <- read_names_clean[keep]

      selected_sequences <- ShortRead::sread(
        selected
      )


      extracted_sequences <- c(
        extracted_sequences,
        DNAStringSet(selected_sequences)
      )

      extracted_names <- c(
        extracted_names,
        selected_names
      )

      found <- unique(
        c(
          found,
          selected_names
        )
      )
    }


    if (
      length(found) == length(target_ids)
    ) {

      cat(
        sprintf(
          "All %d target IDs found in %s.\n",
          length(target_ids),
          mate_label
        )
      )

      break
    }


    if (
      chunk_number %% 25L == 0L
    ) {

      cat(
        sprintf(
          "  %s: processed %d chunks; found %d/%d targets\n",
          mate_label,
          chunk_number,
          length(found),
          length(target_ids)
        )
      )
    }
  }


  # --------------------------------------------------------------------------
  # Remove duplicate occurrences if any.
  #
  # A proper FASTQ should contain each read ID once. If duplicate IDs occur,
  # retain the first sequence and record the duplicate condition.
  # --------------------------------------------------------------------------

  duplicate_ids <- unique(
    extracted_names[
      duplicated(extracted_names)
    ]
  )


  if (length(duplicate_ids) > 0) {

    warning(
      sprintf(
        "%s FASTQ contained duplicate target read IDs: %d",
        mate_label,
        length(duplicate_ids)
      )
    )


    keep_first <- !duplicated(
      extracted_names
    )

    extracted_sequences <- extracted_sequences[
      keep_first
    ]

    extracted_names <- extracted_names[
      keep_first
    ]
  }


  names(extracted_sequences) <- extracted_names


  # --------------------------------------------------------------------------
  # Write extracted sequences as FASTA.
  # --------------------------------------------------------------------------

  Biostrings::writeXStringSet(
    extracted_sequences,
    filepath = output_fasta
  )


  missing <- setdiff(
    target_ids,
    found
  )


  data.frame(
    mate = mate_label,
    target_count = length(target_ids),
    recovered_count = length(found),
    missing_count = length(missing),
    duplicate_target_count = length(duplicate_ids),
    missing_ids = if (
      length(missing) > 0
    ) {
      paste(missing, collapse = ";")
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# Extract R1
# ------------------------------------------------------------------------------

r1_fasta <- file.path(
  DERIVATIVE_DIR,
  "tac102_41pair_R1.fasta"
)


r1_recovery <- extract_fastq_reads(
  fastq_file = R1_FASTQ,
  target_ids = target_ids,
  mate_label = "R1",
  output_fasta = r1_fasta
)


# ------------------------------------------------------------------------------
# Extract R2
# ------------------------------------------------------------------------------

r2_fasta <- file.path(
  DERIVATIVE_DIR,
  "tac102_41pair_R2.fasta"
)


r2_recovery <- extract_fastq_reads(
  fastq_file = R2_FASTQ,
  target_ids = target_ids,
  mate_label = "R2",
  output_fasta = r2_fasta
)


fastq_recovery_summary <- rbind(
  r1_recovery,
  r2_recovery
)


write.csv(
  fastq_recovery_summary,
  file.path(
    DERIVATIVE_DIR,
    "tac102_fastq_recovery_summary.csv"
  ),
  row.names = FALSE
)


cat("\n")
print(fastq_recovery_summary)
cat("\n")


# ==============================================================================
# 8. LOAD EXTRACTED READ SEQUENCES
# ==============================================================================

r1_seqs <- Biostrings::readDNAStringSet(
  r1_fasta
)

r2_seqs <- Biostrings::readDNAStringSet(
  r2_fasta
)


# Normalize names.
names(r1_seqs) <- sub(
  "\\s.*$",
  "",
  names(r1_seqs)
)

names(r2_seqs) <- sub(
  "\\s.*$",
  "",
  names(r2_seqs)
)


# ==============================================================================
# 9. ORIGINAL-HSP-FRAME CODON ANCHORING
#
# The original HSP establishes:
#
#   qstart_raw / qend_raw
#   sstart_raw / send_raw
#
# The candidate sequence must be extended from that HSP without six-frame
# searching.
#
# Forward:
#
#   reference protein position p
#   -> read nucleotide coordinate:
#
#      sstart + (p - qstart) * 3
#
#
# Reverse:
#
# The read sequence is reverse-oriented relative to the reference. The codon
# corresponding to protein position p is recovered from the opposite direction
# and reverse-complemented.
# ==============================================================================


translate_codon <- function(
    codon
) {

  tryCatch(
    as.character(
      Biostrings::translate(
        Biostrings::DNAString(codon),
        if.fuzzy.codon = "X"
      )
    ),
    error = function(e) {
      "X"
    }
  )
}


extract_codon_at_reference_position <- function(
    read_seq,
    qstart_raw,
    sstart_raw,
    send_raw,
    protein_position
) {

  read_len <- length(read_seq)

  forward <- (
    sstart_raw < send_raw
  )


  if (forward) {

    nt_start <- (
      sstart_raw +
      (protein_position - qstart_raw) * 3L
    )

    nt_end <- nt_start + 2L


    if (
      nt_start < 1L ||
      nt_end > read_len
    ) {

      return(
        list(
          success = FALSE,
          codon = NA_character_,
          aa = NA_character_,
          nt_start = NA_integer_,
          nt_end = NA_integer_
        )
      )
    }


    codon <- as.character(
      Biostrings::subseq(
        read_seq,
        start = nt_start,
        end = nt_end
      )
    )


  } else {

    # On the reverse strand, sstart is the highest-coordinate end of the
    # first aligned codon in the read-relative sequence.
    #
    # Move toward lower read coordinates as protein position increases.

    nt_high <- (
      sstart_raw -
      (protein_position - qstart_raw) * 3L
    )

    nt_start <- nt_high - 2L
    nt_end <- nt_high


    if (
      nt_start < 1L ||
      nt_end > read_len
    ) {

      return(
        list(
          success = FALSE,
          codon = NA_character_,
          aa = NA_character_,
          nt_start = NA_integer_,
          nt_end = NA_integer_
        )
      )
    }


    codon_raw <- Biostrings::subseq(
      read_seq,
      start = nt_start,
      end = nt_end
    )


    codon <- as.character(
      Biostrings::reverseComplement(
        codon_raw
      )
    )
  }


  aa <- translate_codon(
    codon
  )


  list(
    success = TRUE,
    codon = codon,
    aa = aa,
    nt_start = nt_start,
    nt_end = nt_end
  )
}


# ==============================================================================
# 10. DETERMINE GAP-FACING READ SEGMENT
#
# This is deliberately done in READ coordinates.
#
# For an N-side read:
#
#   Its HSP terminates at the N-flank boundary.
#   The relevant unaligned segment is the portion of the read immediately
#   beyond that HSP, toward the TAC102 gap.
#
# For a C-side read:
#
#   Its HSP begins at the C-flank boundary.
#   The relevant unaligned segment is the portion immediately before that HSP,
#   toward the TAC102 gap.
#
# Strand orientation determines whether "toward the gap" corresponds to the
# nucleotide coordinates above or below the HSP.
# ==============================================================================


get_gap_facing_nt <- function(
    read_seq,
    hsp,
    side
) {

  read_len <- length(read_seq)

  sstart <- as.integer(
    hsp$sstart_raw
  )

  send <- as.integer(
    hsp$send_raw
  )


  forward <- (
    sstart < send
  )


  if (side == "N") {

    # N-side reference coordinates increase toward the gap.
    #
    # Forward HSP:
    #   gap-facing read segment is after send.
    #
    # Reverse HSP:
    #   gap-facing read segment is before send in read coordinates.

    if (forward) {

      candidate_start <- send + 1L
      candidate_end <- read_len

      if (
        candidate_start > candidate_end
      ) {

        return(
          list(
            sequence = DNAString(""),
            nt_start = NA_integer_,
            nt_end = NA_integer_,
            orientation = "forward"
          )
        )
      }

      seq <- Biostrings::subseq(
        read_seq,
        start = candidate_start,
        end = candidate_end
      )

    } else {

      candidate_start <- 1L
      candidate_end <- send - 1L

      if (
        candidate_end < candidate_start
      ) {

        return(
          list(
            sequence = DNAString(""),
            nt_start = NA_integer_,
            nt_end = NA_integer_,
            orientation = "reverse"
          )
        )
      }

      seq <- Biostrings::reverseComplement(
        Biostrings::subseq(
          read_seq,
          start = candidate_start,
          end = candidate_end
        )
      )

      # After reverse-complementing, sequence orientation is now toward the
      # reference/gap direction.
    }


  } else {

    # C-side reference coordinates increase away from the gap.
    #
    # Forward HSP:
    #   gap-facing sequence is before sstart.
    #
    # Reverse HSP:
    #   gap-facing sequence is after sstart.

    if (forward) {

      candidate_start <- 1L
      candidate_end <- sstart - 1L

      if (
        candidate_end < candidate_start
      ) {

        return(
          list(
            sequence = DNAString(""),
            nt_start = NA_integer_,
            nt_end = NA_integer_,
            orientation = "forward"
          )
        )
      }

      seq <- Biostrings::reverseComplement(
        Biostrings::subseq(
          read_seq,
          start = candidate_start,
          end = candidate_end
        )
      )

    } else {

      candidate_start <- send + 1L
      candidate_end <- read_len

      if (
        candidate_start > candidate_end
      ) {

        return(
          list(
            sequence = DNAString(""),
            nt_start = NA_integer_,
            nt_end = NA_integer_,
            orientation = "reverse"
          )
        )
      }

      seq <- Biostrings::subseq(
        read_seq,
        start = candidate_start,
        end = candidate_end
      )
    }
  }


  list(
    sequence = seq,
    nt_start = candidate_start,
    nt_end = candidate_end,
    orientation = if (forward) {
      "forward"
    } else {
      "reverse"
    }
  )
}


# ==============================================================================
# 11. CODON-ANCHORED GAP EXTENSION
#
# Instead of blindly translating the entire unaligned nucleotide tail, we
# identify which TAC102 protein positions can be reached from the original
# HSP without inventing a new frame.
#
# N-side:
#
#   qend -> GAP_END
#
# C-side:
#
#   qstart -> GAP_START
#
# Only codons physically present in the read are retained.
# ==============================================================================


extend_into_gap <- function(
    read_seq,
    hsp,
    side
) {

  qstart <- as.integer(
    hsp$qstart_raw
  )

  qend <- as.integer(
    hsp$qend_raw
  )

  sstart <- as.integer(
    hsp$sstart_raw
  )

  send <- as.integer(
    hsp$send_raw
  )


  forward <- (
    sstart < send
  )


  if (side == "N") {

    protein_positions <- seq.int(
      from = qend + 1L,
      to = GAP_END,
      by = 1L
    )

  } else {

    protein_positions <- seq.int(
      from = qstart - 1L,
      to = GAP_START,
      by = -1L
    )
  }


  protein_positions <- protein_positions[
    protein_positions >= GAP_START &
    protein_positions <= GAP_END
  ]


  if (length(protein_positions) == 0) {

    return(
      data.frame(
        protein_position = integer(0),
        codon = character(0),
        aa = character(0),
        nt_start = integer(0),
        nt_end = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }


  rows <- vector(
    "list",
    length(protein_positions)
  )


  for (i in seq_along(protein_positions)) {

    p <- protein_positions[i]


    x <- extract_codon_at_reference_position(
      read_seq = read_seq,
      qstart_raw = qstart,
      sstart_raw = sstart,
      send_raw = send,
      protein_position = p
    )


    if (!x$success) {
      break
    }


    rows[[i]] <- data.frame(
      protein_position = p,
      codon = x$codon,
      aa = x$aa,
      nt_start = x$nt_start,
      nt_end = x$nt_end,
      stringsAsFactors = FALSE
    )
  }


  rows <- rows[
    !vapply(
      rows,
      is.null,
      logical(1)
    )
  ]


  if (length(rows) == 0) {

    return(
      data.frame(
        protein_position = integer(0),
        codon = character(0),
        aa = character(0),
        nt_start = integer(0),
        nt_end = integer(0),
        stringsAsFactors = FALSE
      )
    )
  }


  do.call(
    rbind,
    rows
  )
}


# ==============================================================================
# 12. DESCRIPTIVE SIX-FRAME TRANSLATION
#
# This does NOT affect classification.
#
# It is retained solely as descriptive information about the candidate
# nucleotide segment.
# ==============================================================================


six_frame_descriptive <- function(
    seq
) {

  if (
    length(seq) == 0
  ) {
    return(
      list(
        best_frame = NA_character_,
        best_length = 0L
      )
    )
  }


  seq <- Biostrings::DNAString(
    as.character(seq)
  )


  orientations <- list(
    forward = seq,
    reverse = Biostrings::reverseComplement(seq)
  )


  best_length <- 0L
  best_frame <- NA_character_


  for (
    orientation_name in names(orientations)
  ) {

    s <- orientations[[orientation_name]]


    for (
      frame in 0:2
    ) {

      if (
        length(s) <= frame
      ) {
        next
      }


      usable <- length(s) - frame

      usable <- usable -
        (usable %% 3L)


      if (
        usable < 3L
      ) {
        next
      }


      frame_seq <- Biostrings::subseq(
        s,
        start = frame + 1L,
        width = usable
      )


      tr <- tryCatch(
        as.character(
          Biostrings::translate(
            frame_seq,
            if.fuzzy.codon = "X"
          )
        ),
        error = function(e) ""
      )


      aa_non_x <- sum(
        strsplit(
          tr,
          "",
          fixed = TRUE
        )[[1]] != "X"
      )


      if (
        aa_non_x > best_length
      ) {

        best_length <- aa_non_x

        best_frame <- paste0(
          orientation_name,
          ":frame",
          frame
        )
      }
    }
  }


  list(
    best_frame = best_frame,
    best_length = best_length
  )
}


# ==============================================================================
# 13. COMPOSITION DESCRIPTOR
#
# Descriptive only.
# ==============================================================================


composition_similarity <- function(
    aa_seq
) {

  if (
    is.na(aa_seq) ||
    !nzchar(aa_seq)
  ) {
    return(NA_real_)
  }


  chars <- strsplit(
    aa_seq,
    "",
    fixed = TRUE
  )[[1]]


  k_fraction <- mean(
    chars == "K"
  )


  # Gap-region lysine fraction previously established in Step 4B.1.
  GAP_K_FRACTION <- 0.369


  round(
    k_fraction / GAP_K_FRACTION,
    3
  )
}


# ==============================================================================
# 14. BUILD READ SEQUENCE LOOKUP
# ==============================================================================


get_read_sequence <- function(
    mate,
    read_id
) {

  if (mate == "R1") {

    seq <- r1_seqs[[read_id]]

  } else if (mate == "R2") {

    seq <- r2_seqs[[read_id]]

  } else {

    seq <- NULL
  }


  if (
    is.null(seq)
  ) {
    return(NULL)
  }


  seq
}


# ==============================================================================
# 15. ANALYSE ALL BOTH-FLANK PAIRS
# ==============================================================================


cat(
  "\nBeginning TAC102 gap-facing sequence analysis...\n\n"
)


results <- list()

result_counter <- 0L


for (
  i in seq_len(nrow(both_flank))
) {

  pid <- both_flank$pair_id[i]

  n_mate <- both_flank$n_mate[i]

  c_mate <- both_flank$c_mate[i]


  n_hsp <- get_flank_hsp(
    read_id = pid,
    mate = n_mate,
    side = "N"
  )


  c_hsp <- get_flank_hsp(
    read_id = pid,
    mate = c_mate,
    side = "C"
  )


  for (
    side in c("N", "C")
  ) {

    hsp <- if (
      side == "N"
    ) {
      n_hsp
    } else {
      c_hsp
    }


    mate <- if (
      side == "N"
    ) {
      n_mate
    } else {
      c_mate
    }


    read_seq <- get_read_sequence(
      mate = mate,
      read_id = pid
    )


    result_counter <- result_counter + 1L


    row <- data.frame(
      pair_id = pid,
      mate = mate,
      orientation_role = side,

      full_read_length = if (
        !is.null(read_seq)
      ) {
        length(read_seq)
      } else {
        NA_integer_
      },

      hsp_qstart = if (
        !is.null(hsp)
      ) {
        hsp$qstart_raw
      } else {
        NA_integer_
      },

      hsp_qend = if (
        !is.null(hsp)
      ) {
        hsp$qend_raw
      } else {
        NA_integer_
      },

      hsp_sstart = if (
        !is.null(hsp)
      ) {
        hsp$sstart_raw
      } else {
        NA_integer_
      },

      hsp_send = if (
        !is.null(hsp)
      ) {
        hsp$send_raw
      } else {
        NA_integer_
      },

      hsp_bitscore = if (
        !is.null(hsp)
      ) {
        hsp$bitscore
      } else {
        NA_real_
      },

      hsp_strand = if (
        !is.null(hsp)
      ) {
        if (
          hsp$sstart_raw < hsp$send_raw
        ) {
          "forward"
        } else {
          "reverse"
        }
      } else {
        NA_character_
      },

      gap_facing_nt_start = NA_integer_,
      gap_facing_nt_end = NA_integer_,
      gap_facing_nt_length = 0L,
      gap_facing_sequence = NA_character_,

      translated_candidate = NA_character_,

      gap_position_start = NA_integer_,
      gap_position_end = NA_integer_,

      aligned_aa = 0L,
      pident_aa = NA_real_,

      six_frame_best_frame = NA_character_,
      six_frame_descriptive_length = NA_integer_,

      composition_similarity_ratio = NA_real_,

      nucleotide_criterion_note =
        "No TAC102 nucleotide reference exists; codon-anchored original-HSP-frame extension substituted for literal nucleotide-reference alignment.",

      evidence_class = "D",

      class_reason = "",

      stringsAsFactors = FALSE
    )


    # --------------------------------------------------------------------------
    # Missing HSP
    # --------------------------------------------------------------------------

    if (
      is.null(hsp)
    ) {

      row$class_reason <-
        "No flank-relevant raw HSP found."

      results[[result_counter]] <- row

      next
    }


    # --------------------------------------------------------------------------
    # Missing read sequence
    # --------------------------------------------------------------------------

    if (
      is.null(read_seq)
    ) {

      row$class_reason <-
        "Flank-relevant HSP exists but corresponding FASTQ read was not recovered."

      results[[result_counter]] <- row

      next
    }


    # --------------------------------------------------------------------------
    # Determine actual gap-facing nucleotide segment.
    # --------------------------------------------------------------------------

    gap_nt <- get_gap_facing_nt(
      read_seq = read_seq,
      hsp = hsp,
      side = side
    )


    gap_seq <- gap_nt$sequence


    row$gap_facing_nt_start <-
      gap_nt$nt_start

    row$gap_facing_nt_end <-
      gap_nt$nt_end

    row$gap_facing_nt_length <-
      length(gap_seq)


    if (
      length(gap_seq) > 0
    ) {

      row$gap_facing_sequence <-
        as.character(gap_seq)

    } else {

      row$gap_facing_sequence <-
        ""
    }


    # --------------------------------------------------------------------------
    # Descriptive six-frame analysis.
    # --------------------------------------------------------------------------

    sf <- six_frame_descriptive(
      gap_seq
    )


    row$six_frame_best_frame <-
      sf$best_frame

    row$six_frame_descriptive_length <-
      sf$best_length


    # --------------------------------------------------------------------------
    # Minimum candidate length criterion.
    # --------------------------------------------------------------------------

    if (
      length(gap_seq) < MINIMUM_CANDIDATE_NT
    ) {

      row$class_reason <- sprintf(
        "Gap-facing nucleotide segment is < %d nt (%d nt).",
        MINIMUM_CANDIDATE_NT,
        length(gap_seq)
      )

      results[[result_counter]] <- row

      next
    }


    # --------------------------------------------------------------------------
    # Original-HSP-frame, codon-anchored extension.
    # --------------------------------------------------------------------------

    ext <- extend_into_gap(
      read_seq = read_seq,
      hsp = hsp,
      side = side
    )


    if (
      nrow(ext) == 0
    ) {

      row$class_reason <-
        ">=30 nt candidate exists, but no codon-anchored extension into TAC102 gap could be established from the original HSP."

      results[[result_counter]] <- row

      next
    }


    # --------------------------------------------------------------------------
    # Keep only codon-anchored positions corresponding to the actual gap.
    # --------------------------------------------------------------------------

    ext <- ext[
      ext$protein_position >= GAP_START &
      ext$protein_position <= GAP_END,
      ,
      drop = FALSE
    ]


    if (
      nrow(ext) == 0
    ) {

      row$class_reason <-
        "No codon-anchored positions within TAC102 649-732 were recoverable."

      results[[result_counter]] <- row

      next
    }


    row$gap_position_start <-
      min(ext$protein_position)

    row$gap_position_end <-
      max(ext$protein_position)


    row$translated_candidate <-
      paste(
        ext$aa,
        collapse = ""
      )


    # --------------------------------------------------------------------------
    # Direct positional comparison against TAC102 protein reference.
    #
    # No independent alignment is performed.
    # No six-frame search is used.
    # No cherry-picking is performed.
    # --------------------------------------------------------------------------

    ref_chars <- strsplit(
      tac_ref,
      "",
      fixed = TRUE
    )[[1]]


    valid_positions <- ext$protein_position[
      ext$protein_position >= 1 &
      ext$protein_position <= length(ref_chars)
    ]


    if (
      length(valid_positions) == 0
    ) {

      row$class_reason <-
        "Candidate positions fall outside the TAC102 reference."

      results[[result_counter]] <- row

      next
    }


    ext_valid <- ext[
      ext$protein_position %in% valid_positions,
      ,
      drop = FALSE
    ]


    ref_aa <- ref_chars[
      ext_valid$protein_position
    ]


    read_aa <- ext_valid$aa


    interpretable <- (
      read_aa != "X" &
      !is.na(read_aa) &
      !is.na(ref_aa)
    )


    n_interpretable <- sum(
      interpretable
    )


    if (
      n_interpretable > 0
    ) {

      matches <- (
        read_aa[interpretable] ==
        ref_aa[interpretable]
      )

      n_match <- sum(
        matches
      )

      pct_id <- (
        n_match /
        n_interpretable
      ) * 100

    } else {

      n_match <- 0L
      pct_id <- 0
    }


    row$aligned_aa <-
      n_interpretable

    row$pident_aa <-
      round(
        pct_id,
        1
      )


    row$composition_similarity_ratio <-
      composition_similarity(
        row$translated_candidate
      )


    # --------------------------------------------------------------------------
    # PRE-REGISTERED CLASSIFICATION
    # --------------------------------------------------------------------------

    if (
      n_match >= 10L &&
      pct_id >= 50
    ) {

      row$evidence_class <-
        "A"

      row$class_reason <- sprintf(
        "%d/%d interpretable aa match (%.1f%%), original-HSP-frame codon-anchored extension. NOTE: literal nucleotide-reference criterion unavailable because TAC102 CDS reference is not present.",
        n_match,
        n_interpretable,
        pct_id
      )


    } else if (
      n_match >= 6L &&
      pct_id >= 30
    ) {

      row$evidence_class <-
        "B"

      row$class_reason <- sprintf(
        "%d/%d interpretable aa match (%.1f%%), below Class A threshold. NOTE: literal nucleotide-reference criterion unavailable because TAC102 CDS reference is not present.",
        n_match,
        n_interpretable,
        pct_id
      )


    } else {

      row$evidence_class <-
        "C"

      row$class_reason <- sprintf(
        "%d/%d interpretable aa match (%.1f%%), no qualifying A/B correspondence.",
        n_match,
        n_interpretable,
        pct_id
      )
    }


    results[[result_counter]] <- row
  }


  if (
    i %% 5L == 0L ||
    i == nrow(both_flank)
  ) {

    cat(
      sprintf(
        "Analysed %d/%d both-flank pairs.\n",
        i,
        nrow(both_flank)
      )
    )
  }
}


# ==============================================================================
# 16. COMBINE RESULTS
# ==============================================================================

results_df <- do.call(
  rbind,
  results
)


# ==============================================================================
# 17. OUTPUT FULL RESULTS
# ==============================================================================

full_results_file <- file.path(
  PATHS$reports,
  "tac102_41pairs_unaligned_tails.csv"
)


write.csv(
  results_df,
  full_results_file,
  row.names = FALSE
)


# ==============================================================================
# 18. COMPACT GAP-EVIDENCE TABLE
# ==============================================================================

gap_evidence_columns <- c(
  "pair_id",
  "mate",
  "orientation_role",
  "full_read_length",
  "hsp_sstart",
  "hsp_send",
  "hsp_strand",
  "gap_facing_nt_start",
  "gap_facing_nt_end",
  "gap_facing_nt_length",
  "gap_position_start",
  "gap_position_end",
  "aligned_aa",
  "pident_aa",
  "evidence_class",
  "class_reason"
)


gap_evidence_df <- results_df[
  ,
  gap_evidence_columns,
  drop = FALSE
]


gap_evidence_file <- file.path(
  PATHS$reports,
  "tac102_41pairs_gap_evidence.csv"
)


write.csv(
  gap_evidence_df,
  gap_evidence_file,
  row.names = FALSE
)


# ==============================================================================
# 19. CLASS DISTRIBUTION
#
# IMPORTANT:
# Distribution is reported across ALL 82 candidate reads:
#   41 pairs x 2 flank reads
#
# We do NOT report only the strongest read from each pair.
# ==============================================================================


class_counts <- table(
  factor(
    results_df$evidence_class,
    levels = c(
      "A",
      "B",
      "C",
      "D"
    )
  )
)


pair_class_summary <- aggregate(
  evidence_class ~ pair_id,
  data = results_df,
  FUN = function(x) {
    paste(
      unique(x),
      collapse = ";"
    )
  }
)


pair_class_summary_file <- file.path(
  PATHS$reports,
  "tac102_41pairs_class_summary.csv"
)


write.csv(
  pair_class_summary,
  pair_class_summary_file,
  row.names = FALSE
)


# ==============================================================================
# 20. TECHNICAL SUMMARY
# ==============================================================================

technical_counts <- list(

  both_flank_pairs =
    nrow(both_flank),

  candidate_reads =
    nrow(results_df),

  reads_with_sequence =
    sum(
      !is.na(
        results_df$full_read_length
      )
    ),

  reads_missing_sequence =
    sum(
      is.na(
        results_df$full_read_length
      )
    ),

  candidates_ge_30nt =
    sum(
      results_df$gap_facing_nt_length >=
      MINIMUM_CANDIDATE_NT,
      na.rm = TRUE
    ),

  candidates_lt_30nt =
    sum(
      results_df$gap_facing_nt_length <
      MINIMUM_CANDIDATE_NT,
      na.rm = TRUE
    ),

  class_A =
    as.integer(
      class_counts["A"]
    ),

  class_B =
    as.integer(
      class_counts["B"]
    ),

  class_C =
    as.integer(
      class_counts["C"]
    ),

  class_D =
    as.integer(
      class_counts["D"]
    )
)


# ==============================================================================
# 21. REPORT
# ==============================================================================

report <- c(

  "==============================================================",

  "TAC102 GAP-FACING SEQUENCE EVIDENCE",

  "Step 4B.3b",

  "==============================================================",

  "",

  sprintf(
    "Date: %s",
    Sys.Date()
  ),

  "",

  "PROJECT ARCHITECTURE",

  "  Analysis type: independent R/Bioconductor raw-read analysis",

  "  External executables: none",

  "  Conda/WSL: not used",

  "",

  "INPUT FASTQ",

  sprintf(
    "  R1: %s",
    R1_FASTQ
  ),

  sprintf(
    "  R2: %s",
    R2_FASTQ
  ),

  "",

  "REFERENCE",

  sprintf(
    "  TAC102 protein length: %d aa",
    tac_ref_length
  ),

  sprintf(
    "  Gap region: %d-%d",
    GAP_START,
    GAP_END
  ),

  "",

  "PAIR RECONSTRUCTION",

  sprintf(
    "  Both-flank pairs: %d",
    technical_counts$both_flank_pairs
  ),

  sprintf(
    "  Candidate flank reads: %d",
    technical_counts$candidate_reads
  ),

  "",

  "FASTQ RECOVERY",

  sprintf(
    "  Reads with sequence recovered: %d",
    technical_counts$reads_with_sequence
  ),

  sprintf(
    "  Reads missing sequence: %d",
    technical_counts$reads_missing_sequence
  ),

  "",

  "CANDIDATE LENGTH",

  sprintf(
    "  Candidates >= %d nt: %d",
    MINIMUM_CANDIDATE_NT,
    technical_counts$candidates_ge_30nt
  ),

  sprintf(
    "  Candidates < %d nt: %d",
    MINIMUM_CANDIDATE_NT,
    technical_counts$candidates_lt_30nt
  ),

  "",

  "PRE-REGISTERED CLASSIFICATION",

  sprintf(
    "  Class A: %d",
    technical_counts$class_A
  ),

  sprintf(
    "  Class B: %d",
    technical_counts$class_B
  ),

  sprintf(
    "  Class C: %d",
    technical_counts$class_C
  ),

  sprintf(
    "  Class D: %d",
    technical_counts$class_D
  ),

  "",

  "THRESHOLDS",

  sprintf(
    "  Minimum candidate: %d nt",
    MINIMUM_CANDIDATE_NT
  ),

  "  Class A: >=10 aa matched AND >=50% identity",

  "  Class B: >=6 aa matched AND >=30% identity, failing A",

  "  Class C: >=30 nt interpretable, failing A/B",

  "  Class D: <30 nt or technically uninterpretable",

  "",

  "FRAME / ORIENTATION",

  "  Mate role determined dynamically from reference coordinates.",

  "  R1/R2 fixed-role assumption is NOT used.",

  "  Original HSP strand is retained.",

  "  Primary comparison uses original-HSP-frame codon anchoring.",

  "  Six-frame translation is descriptive only.",

  "  Composition similarity is descriptive only.",

  "",

  "KNOWN PROTOCOL DEVIATION",

  "  TAC102 has no nucleotide reference in this pipeline.",

  "  Therefore the literal nucleotide-reference alignment criterion",

  "  cannot be performed.",

  "  It is substituted by original-HSP-frame, codon-anchored nucleotide",

  "  extraction followed by direct positional protein comparison.",

  "  This limitation is explicitly retained in the output.",

  "",

  "IMPORTANT INTERPRETIVE LIMITATION",

  "  query_gap_aa is an HSP-boundary coordinate measurement, not a",

  "  direct physical insert-size measurement.",

  "  Step 4B.3b tests whether the gap-facing raw-read sequence contains",

  "  sequence compatible with TAC102 positions 649-732.",

  "  It does not by itself establish TAC102 identity at the nucleotide",

  "  level in the absence of a TAC102 nucleotide reference.",

  "",

  "OUTPUTS",

  sprintf(
    "  Full results: %s",
    full_results_file
  ),

  sprintf(
    "  Gap evidence: %s",
    gap_evidence_file
  ),

  sprintf(
    "  FASTQ recovery summary: %s",
    file.path(
      PATHS$raw_recovery_derivative,
      "tac102_fastq_recovery_summary.csv"
    )
  ),

  "",

  "=============================================================="
)


report_file <- file.path(
  PATHS$reports,
  "tac102_4B3b_gap_sequence_report.txt"
)


writeLines(
  report,
  report_file
)


# ==============================================================================
# 22. LOGGING
# ==============================================================================

log_info(
  sprintf(
    "Step 4B.3b completed. Both-flank pairs: %d; candidate reads: %d; A:%d B:%d C:%d D:%d.",
    technical_counts$both_flank_pairs,
    technical_counts$candidate_reads,
    technical_counts$class_A,
    technical_counts$class_B,
    technical_counts$class_C,
    technical_counts$class_D
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
# 23. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n==============================================================\n",
  " Step 4B.3b COMPLETED\n",
  "==============================================================\n\n",
  sep = ""
)


cat(
  sprintf(
    "Both-flank pairs: %d\n",
    technical_counts$both_flank_pairs
  )
)

cat(
  sprintf(
    "Candidate reads: %d\n",
    technical_counts$candidate_reads
  )
)

cat(
  sprintf(
    "Sequences recovered: %d\n",
    technical_counts$reads_with_sequence
  )
)

cat(
  sprintf(
    "Sequences missing: %d\n",
    technical_counts$reads_missing_sequence
  )
)

cat(
  sprintf(
    "Candidates >=30 nt: %d\n",
    technical_counts$candidates_ge_30nt
  )
)

cat("\nEvidence classes:\n")

print(
  class_counts
)


cat("\nOutput files:\n")

cat(
  "  ",
  full_results_file,
  "\n",
  sep = ""
)

cat(
  "  ",
  gap_evidence_file,
  "\n",
  sep = ""
)

cat(
  "  ",
  report_file,
  "\n",
  sep = ""
)

cat("\n[SUCCESS] TAC102 Step 4B.3b completed.\n")