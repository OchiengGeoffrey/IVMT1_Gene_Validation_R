# ==============================================================================
# 04p2_tac102_fast_positional_validation.R
#
# PROJECT
# -------
# Genomic & Transcriptomic Limits on kDNA Retention in Trypanosoma equiperdum
#
# PHASE
# -----
# Phase II: Candidate Gene Recovery from T. equiperdum IVM-t1 raw reads
#
# PURPOSE
# -------
# Fast, exhaustive, resumable positional validation of TAC102 positions:
#
#     653
#     698
#
# This script uses:
#
#     1. fast nucleotide prefilter
#     2. rigorous six-frame positional validation
#     3. direct FASTQ nucleotide + Phred evidence
#     4. paired-read validation
#     5. resumable chunked evidence
#
# IMPORTANT:
# The target codon in the prefilter is represented as NNN so that the
# prefilter remains allele-neutral.
#
# The rigorous validator rejects ambiguous target codons and stop codons.
# ==============================================================================


# ==============================================================================
# 0. INITIALIZATION
# ==============================================================================

cat(
  "\n",
  "====================================================================\n",
  " 04p2: TAC102 FAST ALL-READ POSITIONAL VALIDATION\n",
  "====================================================================\n",
  "\n",
  sep = ""
)


if (
  !exists("PATHS") ||
  !exists("CONFIG") ||
  !exists("TARGETS") ||
  !exists("LOG_FILE")
) {

  source(
    here::here(
      "raw_read_analysis",
      "scripts",
      "00_setup.R"
    )
  )
}


if (
  !exists("PATHS") ||
  !exists("TARGETS")
) {

  stop(
    "\n04p2 initialization failed.\n",
    "PATHS and/or TARGETS are unavailable after 00_setup.R.\n"
  )
}


# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

TARGET_POSITIONS <- c(
  653L,
  698L
)


FLANK_WINDOW <- 8L

ANCHOR_WIDTH <- (
  2L * FLANK_WINDOW + 1L
)

TARGET_OFFSET_IN_ANCHOR <- (
  FLANK_WINDOW + 1L
)

MIN_ANCHOR_IDENTITY <- 0.75


RUN_MODE <- toupper(
  Sys.getenv(
    "TAC102_RUN_MODE",
    unset = "PILOT"
  )
)


if (
  !RUN_MODE %in% c(
    "PILOT",
    "FULL"
  )
) {

  stop(
    "TAC102_RUN_MODE must be PILOT or FULL.\n"
  )
}


MAX_READS <- if (RUN_MODE == "PILOT") {

  as.integer(
    Sys.getenv(
      "TAC102_MAX_READS",
      unset = "100000"
    )
  )

} else {

  Inf
}


CHUNK_SIZE <- as.integer(
  Sys.getenv(
    "TAC102_CHUNK_SIZE",
    unset = "100000"
  )
)


PREFILTER_MAX_MISMATCH <- as.integer(
  Sys.getenv(
    "TAC102_PREFILTER_MISMATCH",
    unset = "1"
  )
)


VALIDATOR_MAX_AA_MISMATCH <- as.integer(
  Sys.getenv(
    "TAC102_VALIDATOR_MAX_MISMATCH",
    unset = "5"
  )
)


R1_FASTQ <- PATHS$fastq_r1
R2_FASTQ <- PATHS$fastq_r2


# ==============================================================================
# 2. OUTPUT ARCHITECTURE
# ==============================================================================

RUN_ROOT <- file.path(
  PATHS$reports,
  "tac102_4B3p_fast",
  tolower(RUN_MODE)
)


DIRS <- list(

  root = RUN_ROOT,

  manifest = file.path(
    RUN_ROOT,
    "manifest"
  ),

  patterns = file.path(
    RUN_ROOT,
    "patterns"
  ),

  candidates = file.path(
    RUN_ROOT,
    "candidates"
  ),

  validated = file.path(
    RUN_ROOT,
    "validated"
  ),

  paired = file.path(
    RUN_ROOT,
    "paired"
  ),

  summaries = file.path(
    RUN_ROOT,
    "summaries"
  ),

  logs = file.path(
    RUN_ROOT,
    "logs"
  )
)


invisible(
  lapply(
    DIRS,
    function(x) {

      if (!dir.exists(x)) {

        dir.create(
          x,
          recursive = TRUE,
          showWarnings = FALSE
        )
      }
    }
  )
)


MANIFEST_PATH <- file.path(
  DIRS$manifest,
  "chunk_manifest.csv"
)


RUN_METADATA_PATH <- file.path(
  DIRS$manifest,
  "run_metadata.csv"
)


PREFLIGHT_PATH <- file.path(
  DIRS$summaries,
  "preflight_results.csv"
)


PATTERN_CATALOGUE_PATH <- file.path(
  DIRS$patterns,
  "tac102_pattern_catalogue.csv"
)


PRIMARY_OUTPUT_ALL <- file.path(
  DIRS$summaries,
  "tac102_all_validated_evidence.csv"
)


PARTNER_OUTPUT_ALL <- file.path(
  DIRS$summaries,
  "tac102_all_paired_partner_evidence.csv"
)


POSITION_SUMMARY_PATH <- file.path(
  DIRS$summaries,
  "tac102_position_summary.csv"
)


ALLELE_SUMMARY_PATH <- file.path(
  DIRS$summaries,
  "tac102_allele_summary.csv"
)


HAPLOTYPE_SUMMARY_PATH <- file.path(
  DIRS$summaries,
  "tac102_single_read_haplotype_summary.csv"
)


PAIRED_SUMMARY_PATH <- file.path(
  DIRS$summaries,
  "tac102_paired_read_summary.csv"
)


PILOT_COMPARISON_PATH <- file.path(
  DIRS$summaries,
  "tac102_fast_vs_legacy_pilot_comparison.csv"
)


REPORT_PATH <- file.path(
  DIRS$summaries,
  "tac102_4B3p_fast_report.txt"
)


# ==============================================================================
# 3. LOGGING
# ==============================================================================

FAST_LOG_FILE <- file.path(
  DIRS$logs,
  paste0(
    "tac102_fast_",
    format(
      Sys.time(),
      "%Y%m%d_%H%M%S"
    ),
    ".log"
  )
)


log_fast <- function(message) {

  line <- paste0(
    "[",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    "] ",
    message
  )

  cat(line, "\n")

  cat(
    line,
    "\n",
    file = FAST_LOG_FILE,
    append = TRUE
  )
}


log_fast(
  "Starting 04p2 TAC102 fast positional validation."
)


# ==============================================================================
# 4. INPUT VALIDATION
# ==============================================================================

if (!file.exists(R1_FASTQ)) {

  stop(
    "R1 FASTQ not found:\n",
    R1_FASTQ
  )
}


if (!file.exists(R2_FASTQ)) {

  stop(
    "R2 FASTQ not found:\n",
    R2_FASTQ
  )
}


if (CHUNK_SIZE < 1000L) {

  stop("CHUNK_SIZE must be >= 1000.\n")
}


if (PREFILTER_MAX_MISMATCH < 0L) {

  stop("PREFILTER_MAX_MISMATCH cannot be negative.\n")
}


if (VALIDATOR_MAX_AA_MISMATCH < 0L) {

  stop("VALIDATOR_MAX_AA_MISMATCH cannot be negative.\n")
}


cat(
  "Run mode: ",
  RUN_MODE,
  "\n",
  "R1 FASTQ: ",
  R1_FASTQ,
  "\n",
  "R2 FASTQ: ",
  R2_FASTQ,
  "\n",
  "Maximum reads / FASTQ: ",
  ifelse(
    is.infinite(MAX_READS),
    "FULL DATASET",
    format(MAX_READS, big.mark = ",")
  ),
  "\n",
  "Chunk size: ",
  format(CHUNK_SIZE, big.mark = ","),
  "\n",
  "Target positions: ",
  paste(TARGET_POSITIONS, collapse = ", "),
  "\n",
  "Anchor width: ",
  ANCHOR_WIDTH,
  " aa\n",
  "Minimum anchor identity: ",
  MIN_ANCHOR_IDENTITY,
  "\n",
  "Prefilter maximum nucleotide mismatches: ",
  PREFILTER_MAX_MISMATCH,
  "\n",
  "Validator maximum AA mismatches: ",
  VALIDATOR_MAX_AA_MISMATCH,
  "\n\n",
  sep = ""
)


# ==============================================================================
# 5. FIND TAC102 PROTEIN REFERENCE
# ==============================================================================

resolve_reference_path <- function(reference_file) {

  candidates <- character()


  if (file.exists(reference_file)) {

    candidates <- c(
      candidates,
      reference_file
    )
  }


  candidates <- c(
    candidates,

    file.path(
      here::here("raw_read_analysis"),
      reference_file
    ),

    file.path(
      here::here("raw_read_analysis", "reference"),
      basename(reference_file)
    )
  )


  if ("reference" %in% names(PATHS)) {

    candidates <- c(
      candidates,
      file.path(
        PATHS$reference,
        basename(reference_file)
      )
    )
  }


  if ("reference_dir" %in% names(PATHS)) {

    candidates <- c(
      candidates,
      file.path(
        PATHS$reference_dir,
        basename(reference_file)
      )
    )
  }


  candidates <- unique(candidates)

  found <- candidates[file.exists(candidates)]


  if (length(found) == 0L) {

    stop(
      "Unable to resolve TAC102 reference FASTA.\n",
      "Metadata path: ",
      reference_file,
      "\n"
    )
  }


  found[[1]]
}


tac_metadata <- TARGETS[
  toupper(TARGETS$gene) == "TAC102",
  ,
  drop = FALSE
]


if (nrow(tac_metadata) != 1L) {

  stop("Expected exactly one TAC102 row in target metadata.\n")
}


TAC_REF_FILE <- resolve_reference_path(
  tac_metadata$reference_file[[1]]
)


TAC_REF_SET <- Biostrings::readAAStringSet(
  TAC_REF_FILE
)


if (length(TAC_REF_SET) != 1L) {

  stop(
    "TAC102 reference FASTA must contain exactly one protein sequence.\n",
    "Found: ",
    length(TAC_REF_SET),
    "\n"
  )
}


tac_ref_chr <- toupper(
  as.character(
    TAC_REF_SET[[1]]
  )
)


# ------------------------------------------------------------------------------
# Remove a terminal stop character if present.
#
# The project coordinate system is the protein sequence without terminal stop.
# ------------------------------------------------------------------------------

if (
  nchar(tac_ref_chr) > 0L &&
  endsWith(tac_ref_chr, "*")
) {

  tac_ref_chr <- substr(
    tac_ref_chr,
    1L,
    nchar(tac_ref_chr) - 1L
  )
}


TAC_REF_LENGTH <- nchar(tac_ref_chr)


cat(
  "TAC102 reference: ",
  TAC_REF_FILE,
  "\n",
  "TAC102 reference length: ",
  TAC_REF_LENGTH,
  " aa\n\n",
  sep = ""
)


# ==============================================================================
# 6. REFERENCE IDENTITY LOCKS
# ==============================================================================

reference_653 <- substr(tac_ref_chr, 653L, 653L)
reference_698 <- substr(tac_ref_chr, 698L, 698L)


if (reference_653 != "K") {

  stop(
    "REFERENCE LOCK FAILED: TAC102 position 653 is not K.\n",
    "Observed: ",
    reference_653,
    "\n"
  )
}


if (reference_698 != "L") {

  stop(
    "REFERENCE LOCK FAILED: TAC102 position 698 is not L.\n",
    "Observed: ",
    reference_698,
    "\n"
  )
}


# ==============================================================================
# 7. ANCHORS
# ==============================================================================

get_anchor <- function(reference, position, flank = 8L) {

  start <- position - flank
  end <- position + flank


  if (start < 1L || end > nchar(reference)) {

    stop("Anchor outside reference bounds.\n")
  }


  substr(reference, start, end)
}


ANCHORS <- setNames(
  lapply(
    TARGET_POSITIONS,
    function(p) {

      get_anchor(
        tac_ref_chr,
        p,
        FLANK_WINDOW
      )
    }
  ),
  as.character(TARGET_POSITIONS)
)


for (p in TARGET_POSITIONS) {

  cat(
    p,
    " anchor: ",
    p - FLANK_WINDOW,
    " - ",
    p + FLANK_WINDOW,
    " | ",
    ANCHORS[[as.character(p)]],
    "\n",
    sep = ""
  )
}


cat("\n")


# ==============================================================================
# 8. INTERNAL-CDS TRANSLATION
# ==============================================================================

translate_codon <- function(codon) {

  codon <- toupper(as.character(codon))


  if (
    length(codon) != 1L ||
    nchar(codon) != 3L ||
    grepl("[^ACGT]", codon)
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


# ==============================================================================
# 9. STANDARD GENETIC CODE / CODON MAP
# ==============================================================================

GENETIC_CODONS <- split(
  names(Biostrings::GENETIC_CODE),
  unname(Biostrings::GENETIC_CODE)
)


GENETIC_CODONS <- GENETIC_CODONS[
  names(GENETIC_CODONS) != "*"
]


BASESET_TO_IUPAC <- c(

  "A" = "A",
  "C" = "C",
  "G" = "G",
  "T" = "T",

  "AG" = "R",
  "CT" = "Y",
  "CG" = "S",
  "AT" = "W",
  "GT" = "K",
  "AC" = "M",

  "ACG" = "V",
  "ACT" = "H",
  "AGT" = "D",
  "CGT" = "B",

  "ACGT" = "N"
)


aa_to_iupac_codon <- function(aa) {

  codons <- GENETIC_CODONS[[aa]]


  if (is.null(codons)) {

    stop(
      "No codon mapping for amino acid: ",
      aa
    )
  }


  positions <- lapply(
    seq_len(3L),
    function(j) {

      sort(
        unique(
          substr(codons, j, j)
        )
      )
    }
  )


  codes <- vapply(
    positions,
    function(x) {

      key <- paste(x, collapse = "")

      if (!key %in% names(BASESET_TO_IUPAC)) {

        stop(
          "Unable to convert nucleotide set to IUPAC code: ",
          key
        )
      }

      BASESET_TO_IUPAC[[key]]
    },
    character(1)
  )


  paste0(codes, collapse = "")
}


# ==============================================================================
# 10. BUILD CODON-DEGENERATE NUCLEOTIDE PATTERN
# ==============================================================================

build_nucleotide_pattern <- function(anchor, target_offset) {

  aa <- strsplit(anchor, split = "")[[1]]


  codons <- vapply(
    aa,
    aa_to_iupac_codon,
    character(1)
  )


  codons[target_offset] <- "NNN"


  paste0(codons, collapse = "")
}


PATTERN_ROWS <- list()


for (p in TARGET_POSITIONS) {

  anchor <- ANCHORS[[as.character(p)]]


  pattern_forward <- build_nucleotide_pattern(
    anchor,
    TARGET_OFFSET_IN_ANCHOR
  )


  pattern_reverse <- as.character(
    Biostrings::reverseComplement(
      Biostrings::DNAString(pattern_forward)
    )
  )


  PATTERN_ROWS[[length(PATTERN_ROWS) + 1L]] <- data.frame(

    target_position = p,

    target_reference_aa = substr(tac_ref_chr, p, p),

    anchor_start = p - FLANK_WINDOW,

    anchor_end = p + FLANK_WINDOW,

    anchor_peptide = anchor,

    target_offset_in_anchor = TARGET_OFFSET_IN_ANCHOR,

    pattern_type = "protein_derived_codon_degenerate",

    target_codon_pattern = "NNN",

    forward_pattern = pattern_forward,

    reverse_complement_pattern = pattern_reverse,

    pattern_length_nt = nchar(pattern_forward),

    prefilter_max_mismatch = PREFILTER_MAX_MISMATCH,

    stringsAsFactors = FALSE
  )
}


PATTERN_CATALOGUE <- do.call(rbind, PATTERN_ROWS)


write.csv(
  PATTERN_CATALOGUE,
  PATTERN_CATALOGUE_PATH,
  row.names = FALSE
)


cat("Reference-derived nucleotide prefilter patterns:\n\n")


for (i in seq_len(nrow(PATTERN_CATALOGUE))) {

  cat(
    "Position ",
    PATTERN_CATALOGUE$target_position[i],
    "\n",
    "  peptide: ",
    PATTERN_CATALOGUE$anchor_peptide[i],
    "\n",
    "  forward: ",
    PATTERN_CATALOGUE$forward_pattern[i],
    "\n",
    "  reverse: ",
    PATTERN_CATALOGUE$reverse_complement_pattern[i],
    "\n\n",
    sep = ""
  )
}


# ==============================================================================
# 11. RUN METADATA / RESUME LOCK
# ==============================================================================

input_info <- function(path) {

  info <- file.info(path)


  data.frame(

    path = normalizePath(path, winslash = "/", mustWork = TRUE),

    size_bytes = as.numeric(info$size),

    modified = as.character(info$mtime),

    stringsAsFactors = FALSE
  )
}


R1_INFO <- input_info(R1_FASTQ)
R2_INFO <- input_info(R2_FASTQ)


PATTERN_SIGNATURE <- paste(
  PATTERN_CATALOGUE$forward_pattern,
  PATTERN_CATALOGUE$reverse_complement_pattern,
  collapse = "|"
)


CURRENT_METADATA <- data.frame(

  run_mode = RUN_MODE,

  max_reads = ifelse(
    is.infinite(MAX_READS),
    "FULL",
    as.character(MAX_READS)
  ),

  chunk_size = as.character(CHUNK_SIZE),

  prefilter_max_mismatch = as.character(PREFILTER_MAX_MISMATCH),

  validator_max_aa_mismatch = as.character(VALIDATOR_MAX_AA_MISMATCH),

  r1_path = R1_INFO$path,
  r1_size_bytes = R1_INFO$size_bytes,
  r1_modified = R1_INFO$modified,

  r2_path = R2_INFO$path,
  r2_size_bytes = R2_INFO$size_bytes,
  r2_modified = R2_INFO$modified,

  tac102_reference_path = normalizePath(
    TAC_REF_FILE,
    winslash = "/",
    mustWork = TRUE
  ),

  tac102_reference_length = TAC_REF_LENGTH,

  tac102_reference_653 = reference_653,
  tac102_reference_698 = reference_698,

  pattern_signature = PATTERN_SIGNATURE,

  stringsAsFactors = FALSE
)


if (file.exists(RUN_METADATA_PATH)) {

  OLD_METADATA <- read.csv(
    RUN_METADATA_PATH,
    stringsAsFactors = FALSE
  )


  compare_cols <- intersect(
    names(CURRENT_METADATA),
    names(OLD_METADATA)
  )


  mismatch <- vapply(
    compare_cols,
    function(col) {

      as.character(CURRENT_METADATA[[col]][1]) !=
        as.character(OLD_METADATA[[col]][1])
    },
    logical(1)
  )


  if (any(mismatch)) {

    stop(
      "\nRUN CONFIGURATION MISMATCH.\n\n",
      "The existing resumable run has different configuration/input metadata.\n",
      "Do not mix chunks from different runs.\n\n",
      "Start a fresh run directory or restore the original configuration.\n\n",
      "Mismatches:\n",
      paste(compare_cols[mismatch], collapse = "\n"),
      "\n"
    )
  }

} else {

  write.csv(
    CURRENT_METADATA,
    RUN_METADATA_PATH,
    row.names = FALSE
  )
}


# ==============================================================================
# 12. ATOMIC FILE WRITING
# ==============================================================================

atomic_write_csv <- function(x, path) {

  directory <- dirname(path)


  if (!dir.exists(directory)) {

    dir.create(
      directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }


  temp_path <- paste0(
    path,
    ".tmp_",
    Sys.getpid(),
    "_",
    as.integer(runif(1, 1, 1e9))
  )


  write.csv(
    x,
    temp_path,
    row.names = FALSE,
    quote = TRUE
  )


  if (!file.exists(temp_path) || file.info(temp_path)$size <= 0) {

    if (file.exists(temp_path)) {

      file.remove(temp_path)
    }

    stop(
      "Atomic write failed before rename:\n",
      temp_path
    )
  }


  if (file.exists(path)) {

    file.remove(path)
  }


  ok <- file.rename(temp_path, path)


  if (!ok || !file.exists(path)) {

    if (file.exists(temp_path)) {

      file.remove(temp_path)
    }

    stop(
      "Atomic rename failed:\n",
      path
    )
  }


  invisible(path)
}


# ==============================================================================
# 13. MANIFEST
# ==============================================================================

MANIFEST_COLUMNS <- c(

  "chunk_id",
  "read_start",
  "read_end",

  "r1_n_reads",
  "r2_n_reads",

  "r1_prefilter_candidates",
  "r2_prefilter_candidates",

  "r1_primary_validated",
  "r2_primary_validated",

  "paired_partner_validated",

  "elapsed_seconds",

  "status",

  "started_at",
  "completed_at",

  "candidate_output",
  "validated_output",
  "paired_output"
)


init_manifest <- function() {

  if (!file.exists(MANIFEST_PATH)) {

    empty <- as.data.frame(
      setNames(
        lapply(
          MANIFEST_COLUMNS,
          function(x) {

            if (
              x %in% c(
                "chunk_id",
                "read_start",
                "read_end",
                "r1_n_reads",
                "r2_n_reads",
                "r1_prefilter_candidates",
                "r2_prefilter_candidates",
                "r1_primary_validated",
                "r2_primary_validated",
                "paired_partner_validated"
              )
            ) {

              integer()

            } else {

              character()
            }
          }
        ),
        MANIFEST_COLUMNS
      )
    )


    write.csv(
      empty,
      MANIFEST_PATH,
      row.names = FALSE
    )
  }
}


init_manifest()


read_manifest <- function() {

  read.csv(
    MANIFEST_PATH,
    stringsAsFactors = FALSE
  )
}


write_manifest <- function(manifest) {

  atomic_write_csv(manifest, MANIFEST_PATH)
}


chunk_complete <- function(chunk_id) {

  manifest <- read_manifest()


  row <- manifest[
    manifest$chunk_id == chunk_id &
      manifest$status == "COMPLETE",
    ,
    drop = FALSE
  ]


  if (nrow(row) != 1L) {

    return(FALSE)
  }


  output_paths <- c(
    row$candidate_output,
    row$validated_output,
    row$paired_output
  )


  all(file.exists(output_paths)) &&
    all(file.info(output_paths)$size > 0)
}


mark_chunk_started <- function(chunk_id, read_start, read_end, started_at) {

  manifest <- read_manifest()


  manifest <- manifest[
    manifest$chunk_id != chunk_id,
    ,
    drop = FALSE
  ]


  new_row <- data.frame(

    chunk_id = chunk_id,

    read_start = read_start,
    read_end = read_end,

    r1_n_reads = NA_integer_,
    r2_n_reads = NA_integer_,

    r1_prefilter_candidates = NA_integer_,
    r2_prefilter_candidates = NA_integer_,

    r1_primary_validated = NA_integer_,
    r2_primary_validated = NA_integer_,

    paired_partner_validated = NA_integer_,

    elapsed_seconds = NA_character_,

    status = "RUNNING",

    started_at = as.character(started_at),
    completed_at = NA_character_,

    candidate_output = NA_character_,
    validated_output = NA_character_,
    paired_output = NA_character_,

    stringsAsFactors = FALSE
  )


  manifest <- rbind(manifest, new_row)


  write_manifest(manifest)
}


mark_chunk_complete <- function(chunk_id, values) {

  manifest <- read_manifest()


  idx <- which(manifest$chunk_id == chunk_id)


  if (length(idx) != 1L) {

    stop(
      "Cannot update manifest for chunk ",
      chunk_id
    )
  }


  for (nm in names(values)) {

    manifest[idx, nm] <- values[[nm]]
  }


  manifest[idx, "status"] <- "COMPLETE"

  manifest[idx, "completed_at"] <- as.character(Sys.time())


  write_manifest(manifest)
}


# ==============================================================================
# 14. READ-ID NORMALIZATION
# ==============================================================================

normalize_pair_id <- function(ids) {

  x <- sub("\\s.*$", "", as.character(ids))
  x <- sub("(/1|/2)$", "", x)
  x <- sub("([._-])[12]$", "", x)

  x
}


# ==============================================================================
# 15. QUALITY ACCESSOR / DECODER
#
# These helpers support:
#
#   - ShortReadQ objects with a direct quality slot
#   - FastqQuality objects that cannot be coerced by as.character()
#   - fixed-length reads as an integer matrix
#   - variable-length reads as a list of integer vectors
# ==============================================================================

safe_quality <- function(reads) {

  if (
    inherits(reads, "ShortReadQ") &&
    "quality" %in% methods::slotNames(reads)
  ) {

    return(
      methods::slot(reads, "quality")
    )
  }


  stop(
    "Unable to extract quality slot from ShortReadQ object.\n"
  )
}


quality_strings <- function(qobj) {

  # --------------------------------------------------------------------------
  # 1. Direct S4 coercion to character
  # --------------------------------------------------------------------------

  out <- tryCatch(
    {
      x <- methods::as(qobj, "character")

      if (is.character(x)) {
        x
      } else {
        NULL
      }
    },
    error = function(e) {
      NULL
    }
  )


  if (!is.null(out)) {

    return(out)
  }


  # --------------------------------------------------------------------------
  # 2. Coerce to BStringSet, then to character
  # --------------------------------------------------------------------------

  out <- tryCatch(
    {
      b <- methods::as(qobj, "BStringSet")

      x <- methods::as(b, "character")

      if (is.character(x)) {
        x
      } else {
        NULL
      }
    },
    error = function(e) {
      NULL
    }
  )


  if (!is.null(out)) {

    return(out)
  }


  # --------------------------------------------------------------------------
  # 3. Try known S4 slots
  # --------------------------------------------------------------------------

  if (isS4(qobj)) {

    sn <- methods::slotNames(qobj)

    candidate_slots <- c(
      "quality",
      ".Data",
      "data",
      "values"
    )


    for (s in candidate_slots) {

      if (s %in% sn) {

        sv <- methods::slot(qobj, s)


        if (is.character(sv)) {

          return(sv)
        }


        out <- tryCatch(
          {
            x <- methods::as(sv, "character")

            if (is.character(x)) {
              x
            } else {
              NULL
            }
          },
          error = function(e) {
            NULL
          }
        )


        if (!is.null(out)) {

          return(out)
        }


        out <- tryCatch(
          {
            b <- methods::as(sv, "BStringSet")

            x <- methods::as(b, "character")

            if (is.character(x)) {
              x
            } else {
              NULL
            }
          },
          error = function(e) {
            NULL
          }
        )


        if (!is.null(out)) {

          return(out)
        }
      }
    }
  }


  stop(
    "Unable to extract ASCII quality strings from FastqQuality object.\n"
  )
}


safe_quality_data <- function(reads) {

  qobj <- safe_quality(reads)

  expected_n <- length(reads)


  # --------------------------------------------------------------------------
  # 1. Try native coercion first, but only accept it if completely NA-free
  # --------------------------------------------------------------------------

  mat <- tryCatch(
    as(qobj, "matrix"),
    error = function(e) {
      NULL
    }
  )


  if (!is.null(mat)) {

    storage.mode(mat) <- "integer"


    if (!anyNA(mat) && nrow(mat) == expected_n) {

      dimnames(mat) <- NULL

      return(
        list(
          type = "matrix",
          data = mat
        )
      )
    }
  }


  # --------------------------------------------------------------------------
  # 2. Extract ASCII quality strings
  # --------------------------------------------------------------------------

  qchar <- quality_strings(qobj)


  if (length(qchar) != expected_n) {

    stop(
      "Quality string extraction returned ",
      length(qchar),
      " strings but expected ",
      expected_n,
      ".\n"
    )
  }


  if (expected_n == 0L) {

    return(
      list(
        type = "list",
        data = list()
      )
    )
  }


  widths <- nchar(qchar)


  offset <- if (inherits(qobj, "SFastqQuality")) {
    64L
  } else {
    33L
  }


  decode_one <- function(z, off) {

    as.integer(charToRaw(z)) - off
  }


  decode_all <- function(off) {

    if (length(unique(widths)) == 1L) {

      read_width <- widths[[1]]


      m <- t(
        vapply(
          qchar,
          function(z) {
            decode_one(z, off)
          },
          integer(read_width)
        )
      )


      dimnames(m) <- NULL


      return(
        list(
          type = "matrix",
          data = m
        )
      )
    }


    l <- lapply(
      qchar,
      function(z) {
        decode_one(z, off)
      }
    )


    list(
      type = "list",
      data = l
    )
  }


  qd <- decode_all(offset)


  flatten_scores <- function(qd_object) {

    if (qd_object$type == "matrix") {

      return(qd_object$data)
    }


    unlist(
      qd_object$data,
      use.names = FALSE
    )
  }


  scores <- flatten_scores(qd)


  if (any(scores < 0L, na.rm = TRUE)) {

    alt_offset <- if (offset == 33L) 64L else 33L

    qd_alt <- decode_all(alt_offset)

    scores_alt <- flatten_scores(qd_alt)


    if (
      sum(scores_alt < 0L, na.rm = TRUE) <
      sum(scores < 0L, na.rm = TRUE)
    ) {

      qd <- qd_alt
      scores <- scores_alt
    }
  }


  if (anyNA(scores)) {

    stop("FASTQ quality decoding produced NA values.\n")
  }


  if (any(scores < 0L)) {

    stop(
      "FASTQ quality decoding produced negative values.\n",
      "The FASTQ quality offset could not be determined reliably.\n"
    )
  }


  qd
}


quality_vector <- function(qd, i) {

  if (qd$type == "matrix") {

    return(qd$data[i, ])
  }

  qd$data[[i]]
}


quality_nrow <- function(qd) {

  if (qd$type == "matrix") {

    return(nrow(qd$data))
  }

  length(qd$data)
}


quality_range <- function(qd) {

  if (qd$type == "matrix") {

    return(range(qd$data))
  }


  scores <- unlist(
    qd$data,
    use.names = FALSE
  )


  if (length(scores) == 0L) {

    return(c(NA_integer_, NA_integer_))
  }


  range(scores)
}


# ==============================================================================
# 16. RIGOROUS READ VALIDATION
# ==============================================================================

validate_read_against_target <- function(
    seq_obj,
    qual_vec,
    read_id,
    pair_id,
    mate,
    read_index,
    target_position,
    prefilter_count_forward = NA_integer_,
    prefilter_count_reverse = NA_integer_,
    evidence_role = "PRIMARY"
) {

  if (!inherits(seq_obj, "DNAString")) {

    seq_obj <- Biostrings::DNAString(as.character(seq_obj))
  }


  read_length <- nchar(as.character(seq_obj))


  if (length(qual_vec) != read_length) {

    stop(
      "Sequence/quality length mismatch for read ",
      read_id,
      "."
    )
  }


  anchor <- ANCHORS[[as.character(target_position)]]

  target_offset <- TARGET_OFFSET_IN_ANCHOR

  results <- list()

  result_counter <- 0L


  for (orientation in c("fwd", "rev")) {

    oriented_seq <- if (orientation == "fwd") {

      seq_obj

    } else {

      Biostrings::reverseComplement(seq_obj)
    }


    oriented_qual <- if (orientation == "fwd") {

      qual_vec

    } else {

      rev(qual_vec)
    }


    for (frame in 1L:3L) {

      remaining <- read_length - frame + 1L


      if (remaining < 3L) {

        next
      }


      usable_nt <- floor(remaining / 3) * 3


      frame_seq <- Biostrings::subseq(
        oriented_seq,
        start = frame,
        width = usable_nt
      )


      aa_seq <- suppressWarnings(
        Biostrings::translate(
          frame_seq,
          no.init.codon = TRUE,
          if.fuzzy.codon = "X"
        )
      )


      matches <- Biostrings::matchPattern(
        Biostrings::AAString(anchor),
        aa_seq,
        max.mismatch = VALIDATOR_MAX_AA_MISMATCH
      )


      if (length(matches) == 0L) {

        next
      }


      match_starts <- start(matches)
      match_ends <- end(matches)


      for (j in seq_along(match_starts)) {

        aa_start <- match_starts[j]
        aa_end <- match_ends[j]


        if (aa_end > nchar(as.character(aa_seq))) {

          next
        }


        anchor_candidate <- substr(
          as.character(aa_seq),
          aa_start,
          aa_end
        )


        if (nchar(anchor_candidate) != ANCHOR_WIDTH) {

          next
        }


        ref_chars <- strsplit(anchor, "")[[1]]
        cand_chars <- strsplit(anchor_candidate, "")[[1]]

        keep_anchor <- seq_len(ANCHOR_WIDTH) != target_offset


        anchor_identity <- mean(
          cand_chars[keep_anchor] == ref_chars[keep_anchor]
        )


        if (anchor_identity < MIN_ANCHOR_IDENTITY) {

          next
        }


        target_aa_from_anchor <- cand_chars[target_offset]


        # ----------------------------------------------------------------
        # Reject ambiguous or stop residues at the target position.
        # ----------------------------------------------------------------

        if (target_aa_from_anchor %in% c("X", "*")) {

          next
        }


        codon_start_oriented <- frame +
          (aa_start - 1L) * 3L +
          (target_offset - 1L) * 3L


        codon_end_oriented <- codon_start_oriented + 2L


        if (
          codon_start_oriented < 1L ||
          codon_end_oriented > read_length
        ) {

          next
        }


        raw_codon <- Biostrings::subseq(
          oriented_seq,
          start = codon_start_oriented,
          end = codon_end_oriented
        )


        raw_codon_qual <- oriented_qual[
          codon_start_oriented:codon_end_oriented
        ]


        if (orientation == "rev") {

          original_codon_start <- read_length - codon_end_oriented + 1L
          original_codon_end <- read_length - codon_start_oriented + 1L

        } else {

          original_codon_start <- codon_start_oriented
          original_codon_end <- codon_end_oriented
        }


        codon_string <- toupper(as.character(raw_codon))


        # ----------------------------------------------------------------
        # Target codon must be unambiguous A/C/G/T.
        # ----------------------------------------------------------------

        if (grepl("[^ACGT]", codon_string)) {

          next
        }


        target_aa <- translate_codon(codon_string)


        # ----------------------------------------------------------------
        # Reject ambiguous or stop translations at the target position.
        # ----------------------------------------------------------------

        if (target_aa %in% c("X", "*")) {

          next
        }


        if (target_aa != target_aa_from_anchor) {

          next
        }


        expected_reference_aa <- substr(
          tac_ref_chr,
          target_position,
          target_position
        )


        pct_reference <- 100 * as.numeric(
          target_aa == expected_reference_aa
        )


        variant_class <- if (target_aa == expected_reference_aa) {

          "REFERENCE_MATCH"

        } else {

          paste0(
            expected_reference_aa,
            target_position,
            target_aa
          )
        }


        result_counter <- result_counter + 1L


        results[[result_counter]] <- data.frame(

          read_id = read_id,
          pair_id = pair_id,
          mate = mate,
          read_index = read_index,
          evidence_role = evidence_role,

          target_position = target_position,
          orientation = orientation,
          frame = frame,

          anchor_reference_start = target_position - FLANK_WINDOW,
          anchor_reference_end = target_position + FLANK_WINDOW,

          anchor_peptide = anchor,
          anchor_candidate = anchor_candidate,
          anchor_identity = anchor_identity,

          target_aa = target_aa,
          target_codon = codon_string,
          expected_reference_aa = expected_reference_aa,

          pct_reference_aa = pct_reference,
          variant_class = variant_class,

          read_codon_start = original_codon_start,
          read_codon_end = original_codon_end,

          codon_min_qual = min(raw_codon_qual),
          codon_mean_qual = mean(raw_codon_qual),

          read_length = read_length,

          prefilter_forward_count = prefilter_count_forward,
          prefilter_reverse_count = prefilter_count_reverse,

          validator_max_aa_mismatch = VALIDATOR_MAX_AA_MISMATCH,

          stringsAsFactors = FALSE
        )
      }
    }
  }


  if (length(results) == 0L) {

    return(NULL)
  }


  do.call(rbind, results)
}


# ==============================================================================
# 17. FAST NUCLEOTIDE PREFILTER
# ==============================================================================

fast_nucleotide_prefilter <- function(dna_set) {

  if (!inherits(dna_set, "DNAStringSet")) {

    dna_set <- Biostrings::DNAStringSet(dna_set)
  }


  candidate_table <- data.frame(

    read_index = seq_along(dna_set),

    forward_count_653 = 0L,
    reverse_count_653 = 0L,

    forward_count_698 = 0L,
    reverse_count_698 = 0L,

    stringsAsFactors = FALSE
  )


  for (i in seq_len(nrow(PATTERN_CATALOGUE))) {

    p <- PATTERN_CATALOGUE[i, , drop = FALSE]

    target <- p$target_position[[1]]


    forward_counts <- Biostrings::vcountPattern(
      p$forward_pattern[[1]],
      dna_set,
      max.mismatch = PREFILTER_MAX_MISMATCH,
      fixed = "subject"
    )


    reverse_counts <- Biostrings::vcountPattern(
      p$reverse_complement_pattern[[1]],
      dna_set,
      max.mismatch = PREFILTER_MAX_MISMATCH,
      fixed = "subject"
    )


    if (target == 653L) {

      candidate_table$forward_count_653 <- as.integer(forward_counts)
      candidate_table$reverse_count_653 <- as.integer(reverse_counts)

    } else {

      candidate_table$forward_count_698 <- as.integer(forward_counts)
      candidate_table$reverse_count_698 <- as.integer(reverse_counts)
    }
  }


  candidate_table$any_653 <- (
    candidate_table$forward_count_653 > 0L |
      candidate_table$reverse_count_653 > 0L
  )


  candidate_table$any_698 <- (
    candidate_table$forward_count_698 > 0L |
      candidate_table$reverse_count_698 > 0L
  )


  candidate_table$any_candidate <- (
    candidate_table$any_653 |
      candidate_table$any_698
  )


  candidate_table
}


# ==============================================================================
# 18. SYNTHETIC READ GENERATOR
# ==============================================================================

REPRESENTATIVE_CODONS <- vapply(
  GENETIC_CODONS,
  function(x) {
    x[[1]]
  },
  character(1)
)


make_synthetic_read <- function(
    target_position,
    target_aa,
    orientation = "fwd",
    frame = 1L,
    read_length = 120L,
    anchor_start_nt = 10L
) {

  anchor <- ANCHORS[[as.character(target_position)]]

  aa <- strsplit(anchor, "")[[1]]

  aa[TARGET_OFFSET_IN_ANCHOR] <- target_aa


  anchor_dna <- paste0(
    vapply(
      aa,
      function(a) {
        REPRESENTATIVE_CODONS[[a]]
      },
      character(1)
    ),
    collapse = ""
  )


  if (
    anchor_start_nt < 1L ||
    anchor_start_nt + nchar(anchor_dna) - 1L > read_length
  ) {

    stop("Synthetic anchor does not fit inside read.")
  }


  if (((anchor_start_nt - frame) %% 3L) != 0L) {

    stop("Synthetic anchor start is not compatible with requested frame.")
  }


  filler_length_left <- anchor_start_nt - 1L

  filler_length_right <- read_length -
    (anchor_start_nt + nchar(anchor_dna) - 1L)


  raw <- paste0(
    paste(rep("C", filler_length_left), collapse = ""),
    anchor_dna,
    paste(rep("G", filler_length_right), collapse = "")
  )


  raw_dna <- Biostrings::DNAString(raw)


  if (orientation == "rev") {

    raw_dna <- Biostrings::reverseComplement(raw_dna)
  }


  list(

    seq = raw_dna,

    quality = rep(35L, read_length),

    expected_target_aa = target_aa,

    expected_target_position = target_position,

    expected_orientation = orientation,

    expected_frame = frame
  )
}


# ==============================================================================
# 19. PRE-FLIGHT TESTS
# ==============================================================================

run_preflight_tests <- function() {

  results <- list()


  results[["01_DNAString_type"]] <- inherits(
    Biostrings::DNAString("ACGT"),
    "DNAString"
  )


  results[["02_AAString_type"]] <- inherits(
    Biostrings::AAString("M"),
    "AAString"
  )


  results[["03_reference_length_951"]] <- TAC_REF_LENGTH == 951L


  results[["04_reference_653_K"]] <- identical(reference_653, "K")


  results[["05_reference_698_L"]] <- identical(reference_698, "L")


  results[["06_anchor_width_17"]] <- all(
    nchar(unlist(ANCHORS)) == 17L
  )


  results[["07_target_offset_center"]] <- TARGET_OFFSET_IN_ANCHOR == 9L


  results[["08_pattern_length_51nt"]] <- all(
    PATTERN_CATALOGUE$pattern_length_nt == 51L
  )


  valid_iupac <- function(x) {

    !grepl("[^ACGTRYSWKMBDHVN]", x)
  }


  results[["09_forward_patterns_valid_IUPAC"]] <- all(
    vapply(
      PATTERN_CATALOGUE$forward_pattern,
      valid_iupac,
      logical(1)
    )
  )


  results[["10_reverse_patterns_valid_IUPAC"]] <- all(
    vapply(
      PATTERN_CATALOGUE$reverse_complement_pattern,
      valid_iupac,
      logical(1)
    )
  )


  translation_test_codons <- c(

    ATG = "M",
    CTG = "L",
    TTG = "L",
    GTG = "V",
    AAA = "K",
    AAG = "K",
    TTT = "F",
    GGG = "G"
  )


  observed <- vapply(
    names(translation_test_codons),
    translate_codon,
    character(1)
  )


  observed <- unname(observed)


  results[["11_internal_codon_translation"]] <- identical(
    observed,
    unname(translation_test_codons)
  )


  syn_fwd <- make_synthetic_read(
    target_position = 653L,
    target_aa = "M",
    orientation = "fwd",
    frame = 1L
  )


  fwd_prefilter <- fast_nucleotide_prefilter(
    Biostrings::DNAStringSet(syn_fwd$seq)
  )


  results[["12_forward_prefilter_653"]] <- isTRUE(
    fwd_prefilter$any_653[1]
  )


  syn_rev <- make_synthetic_read(
    target_position = 653L,
    target_aa = "M",
    orientation = "rev",
    frame = 1L
  )


  rev_prefilter <- fast_nucleotide_prefilter(
    Biostrings::DNAStringSet(syn_rev$seq)
  )


  results[["13_reverse_prefilter_653"]] <- isTRUE(
    rev_prefilter$any_653[1]
  )


  fwd_validation <- validate_read_against_target(
    seq_obj = syn_fwd$seq,
    qual_vec = syn_fwd$quality,
    read_id = "SYNTHETIC_FWD",
    pair_id = "SYNTHETIC_FWD",
    mate = "R1",
    read_index = 1L,
    target_position = 653L
  )


  results[["14_forward_rigorous_validation"]] <- (
    !is.null(fwd_validation) &&
      any(fwd_validation$target_aa == "M")
  )


  rev_validation <- validate_read_against_target(
    seq_obj = syn_rev$seq,
    qual_vec = syn_rev$quality,
    read_id = "SYNTHETIC_REV",
    pair_id = "SYNTHETIC_REV",
    mate = "R1",
    read_index = 1L,
    target_position = 653L
  )


  results[["15_reverse_rigorous_validation"]] <- (
    !is.null(rev_validation) &&
      any(rev_validation$target_aa == "M")
  )


  if (!is.null(rev_validation)) {

    reverse_rows <- rev_validation[
      rev_validation$target_aa == "M",
      ,
      drop = FALSE
    ]

    results[["16_reverse_coordinate_bounds"]] <- all(
      reverse_rows$read_codon_start >= 1L &
        reverse_rows$read_codon_end <= reverse_rows$read_length
    )

  } else {

    results[["16_reverse_coordinate_bounds"]] <- FALSE
  }


  if (!is.null(fwd_validation)) {

    fwd_rows <- fwd_validation[
      fwd_validation$target_aa == "M",
      ,
      drop = FALSE
    ]

    results[["17_forward_coordinate_bounds"]] <- all(
      fwd_rows$read_codon_start >= 1L &
        fwd_rows$read_codon_end <= fwd_rows$read_length
    )

  } else {

    results[["17_forward_coordinate_bounds"]] <- FALSE
  }


  results[["18_exact_phred_coordinates"]] <- (
    !is.null(fwd_validation) &&
      all(fwd_validation$codon_min_qual == 35L)
  )


  results[["19_reverse_quality_handling"]] <- (
    !is.null(rev_validation) &&
      all(rev_validation$codon_min_qual == 35L)
  )


  syn_698 <- make_synthetic_read(
    target_position = 698L,
    target_aa = "P",
    orientation = "fwd",
    frame = 1L
  )


  validation_698 <- validate_read_against_target(
    seq_obj = syn_698$seq,
    qual_vec = syn_698$quality,
    read_id = "SYNTHETIC_698",
    pair_id = "SYNTHETIC_698",
    mate = "R1",
    read_index = 1L,
    target_position = 698L
  )


  results[["20_detects_698_variant"]] <- (
    !is.null(validation_698) &&
      any(validation_698$target_aa == "P")
  )


  results[["21_rejects_unrelated_read"]] <- {

    unrelated <- Biostrings::DNAString(
      paste(rep("ACGT", 30), collapse = "")
    )


    unrelated_result <- fast_nucleotide_prefilter(
      Biostrings::DNAStringSet(unrelated)
    )


    !isTRUE(unrelated_result$any_candidate[1])
  }


  results[["22_no_init_translation"]] <- (
    translate_codon("CTG") == "L" &&
      translate_codon("TTG") == "L"
  )


  results[["23_pair_id_normalization"]] <- (
    identical(normalize_pair_id("READ.123/1"), "READ.123") &&
      identical(normalize_pair_id("READ.123/2"), "READ.123")
  )


  results[["24_length_guard"]] <- {

    short_seq <- Biostrings::DNAString("ACGTACGTAC")

    short_quality <- rep(30L, 10L)


    out <- validate_read_against_target(
      seq_obj = short_seq,
      qual_vec = short_quality,
      read_id = "SHORT",
      pair_id = "SHORT",
      mate = "R1",
      read_index = 1L,
      target_position = 653L
    )


    is.null(out)
  }


  results[["25_reference_653_assertion"]] <- (
    substr(tac_ref_chr, 653L, 653L) == "K"
  )


  results[["26_reference_698_assertion"]] <- (
    substr(tac_ref_chr, 698L, 698L) == "L"
  )


  data.frame(

    test = names(results),

    passed = vapply(results, isTRUE, logical(1)),

    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 20. RUN PREFLIGHT
# ==============================================================================

cat(
  "\n",
  "====================================================================\n",
  " RUNNING 04p2 PRE-FLIGHT TEST SUITE\n",
  "====================================================================\n",
  "\n",
  sep = ""
)


PREFLIGHT <- run_preflight_tests()


print(PREFLIGHT, row.names = FALSE)


write.csv(
  PREFLIGHT,
  PREFLIGHT_PATH,
  row.names = FALSE
)


if (!all(PREFLIGHT$passed)) {

  failed <- PREFLIGHT$test[!PREFLIGHT$passed]


  cat(
    "\nFAILED PRE-FLIGHT TESTS:\n",
    paste(failed, collapse = "\n"),
    "\n\n",
    sep = ""
  )


  stop(
    "\n04p2 FASTQ PROCESSING REFUSED TO START.\n",
    "Failed pre-flight tests:\n",
    paste0("  - ", failed, collapse = "\n"),
    "\n\nFull pre-flight output:\n  ",
    PREFLIGHT_PATH,
    "\n"
  )
}


cat(
  "\n",
  "====================================================================\n",
  " PRE-FLIGHT PASSED: ALL TESTS\n",
  " FASTQ processing is now permitted.\n",
  "====================================================================\n",
  "\n",
  sep = ""
)


log_fast(
  paste0(
    "Pre-flight passed: ",
    nrow(PREFLIGHT),
    "/",
    nrow(PREFLIGHT),
    " tests."
  )
)


# ==============================================================================
# 21. FASTQ PAIRING PREFLIGHT
# ==============================================================================

check_initial_pairing <- function(r1_path, r2_path, n = 1000L) {

  fq1 <- ShortRead::FastqStreamer(r1_path, n = n)
  fq2 <- ShortRead::FastqStreamer(r2_path, n = n)


  on.exit(try(close(fq1), silent = TRUE), add = TRUE)
  on.exit(try(close(fq2), silent = TRUE), add = TRUE)


  x1 <- ShortRead::yield(fq1)
  x2 <- ShortRead::yield(fq2)


  if (length(x1) == 0L || length(x2) == 0L) {

    return(
      list(
        passed = FALSE,
        n_compared = 0L,
        message = "One FASTQ file returned zero reads during pairing check."
      )
    )
  }


  q1 <- safe_quality(x1)

  q1_class <- class(q1)[1]

  q1_data <- safe_quality_data(x1)

  q1_range <- quality_range(q1_data)


  cat(
    "Quality pre-flight:\n",
    "  R1 class: ",
    q1_class,
    "\n",
    "  R1 Phred range: ",
    q1_range[1],
    " - ",
    q1_range[2],
    "\n",
    sep = ""
  )


  stopifnot(
    quality_nrow(q1_data) == length(x1),
    !anyNA(q1_range),
    q1_range[1] >= 0
  )


  ids1 <- normalize_pair_id(ShortRead::id(x1))
  ids2 <- normalize_pair_id(ShortRead::id(x2))


  n_compare <- min(length(ids1), length(ids2))


  same <- ids1[seq_len(n_compare)] == ids2[seq_len(n_compare)]


  list(

    passed = all(same),

    n_compared = n_compare,

    mismatch_count = sum(!same),

    first_mismatch = if (any(!same)) {
      which(!same)[1]
    } else {
      NA_integer_
    },

    message = if (all(same)) {
      "R1/R2 normalized read IDs agree."
    } else {
      "R1/R2 normalized read IDs do not agree."
    }
  )
}


PAIR_PREFLIGHT <- check_initial_pairing(
  R1_FASTQ,
  R2_FASTQ,
  n = 1000L
)


cat(
  "FASTQ pairing pre-flight:\n",
  "  Compared reads: ",
  PAIR_PREFLIGHT$n_compared,
  "\n",
  "  Mismatches: ",
  PAIR_PREFLIGHT$mismatch_count,
  "\n",
  "  Status: ",
  ifelse(PAIR_PREFLIGHT$passed, "PASSED", "FAILED"),
  "\n\n",
  sep = ""
)


if (!PAIR_PREFLIGHT$passed) {

  stop(
    "FASTQ pairing pre-flight failed.\n",
    PAIR_PREFLIGHT$message,
    "\n"
  )
}


# ==============================================================================
# 22. CHUNK OUTPUT SCHEMAS
# ==============================================================================

empty_candidate_table <- function() {

  data.frame(

    mate = character(),
    read_index = integer(),
    read_id = character(),
    pair_id = character(),

    prefilter_653_forward = integer(),
    prefilter_653_reverse = integer(),

    prefilter_698_forward = integer(),
    prefilter_698_reverse = integer(),

    any_653 = logical(),
    any_698 = logical(),
    any_candidate = logical(),

    stringsAsFactors = FALSE
  )
}


empty_validated_table <- function() {

  data.frame(

    read_id = character(),
    pair_id = character(),
    mate = character(),
    read_index = integer(),
    evidence_role = character(),

    target_position = integer(),
    orientation = character(),
    frame = integer(),

    anchor_reference_start = integer(),
    anchor_reference_end = integer(),

    anchor_peptide = character(),
    anchor_candidate = character(),
    anchor_identity = numeric(),

    target_aa = character(),
    target_codon = character(),
    expected_reference_aa = character(),

    pct_reference_aa = numeric(),
    variant_class = character(),

    read_codon_start = integer(),
    read_codon_end = integer(),

    codon_min_qual = integer(),
    codon_mean_qual = numeric(),

    read_length = integer(),

    prefilter_forward_count = integer(),
    prefilter_reverse_count = integer(),

    validator_max_aa_mismatch = integer(),

    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 23. PROCESS ONE MATE'S PREFILTER
# ==============================================================================

prefilter_mate <- function(reads, mate_label, global_indices) {

  dna <- ShortRead::sread(reads)

  prefilter <- fast_nucleotide_prefilter(dna)

  ids <- normalize_pair_id(ShortRead::id(reads))


  data.frame(

    mate = mate_label,

    read_index = as.integer(global_indices),

    read_id = sub(
      "\\s.*$",
      "",
      as.character(ShortRead::id(reads))
    ),

    pair_id = ids,

    prefilter_653_forward = prefilter$forward_count_653,
    prefilter_653_reverse = prefilter$reverse_count_653,

    prefilter_698_forward = prefilter$forward_count_698,
    prefilter_698_reverse = prefilter$reverse_count_698,

    any_653 = prefilter$any_653,
    any_698 = prefilter$any_698,
    any_candidate = prefilter$any_candidate,

    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# 24. VALIDATE PREFILTER CANDIDATES
# ==============================================================================

validate_primary_candidates <- function(reads, candidate_table) {

  candidate_indices <- which(candidate_table$any_candidate)


  if (length(candidate_indices) == 0L) {

    return(empty_validated_table())
  }


  dna <- ShortRead::sread(reads)


  qual_data <- safe_quality_data(reads)


  stopifnot(
    quality_nrow(qual_data) == length(dna)
  )


  results <- list()

  counter <- 0L


  for (i in candidate_indices) {

    seq_obj <- dna[[i]]

    qual_vec <- quality_vector(qual_data, i)


    read_id <- candidate_table$read_id[i]
    pair_id <- candidate_table$pair_id[i]
    mate <- candidate_table$mate[i]
    read_index <- candidate_table$read_index[i]


    target_positions_for_read <- integer()


    if (candidate_table$any_653[i]) {

      target_positions_for_read <- c(target_positions_for_read, 653L)
    }


    if (candidate_table$any_698[i]) {

      target_positions_for_read <- c(target_positions_for_read, 698L)
    }


    for (target_position in target_positions_for_read) {

      v <- validate_read_against_target(

        seq_obj = seq_obj,
        qual_vec = qual_vec,

        read_id = read_id,
        pair_id = pair_id,
        mate = mate,
        read_index = read_index,

        target_position = target_position,

        prefilter_count_forward = if (target_position == 653L) {
          candidate_table$prefilter_653_forward[i]
        } else {
          candidate_table$prefilter_698_forward[i]
        },

        prefilter_count_reverse = if (target_position == 653L) {
          candidate_table$prefilter_653_reverse[i]
        } else {
          candidate_table$prefilter_698_reverse[i]
        },

        evidence_role = "PRIMARY"
      )


      if (!is.null(v)) {

        counter <- counter + 1L

        results[[counter]] <- v
      }
    }
  }


  if (length(results) == 0L) {

    return(empty_validated_table())
  }


  do.call(rbind, results)
}


# ==============================================================================
# 25. VALIDATE PAIRED PARTNERS
# ==============================================================================

validate_paired_partners <- function(
    reads_r1,
    reads_r2,
    candidate_r1,
    candidate_r2
) {

  candidate_pairs <- unique(
    c(
      candidate_r1$pair_id[candidate_r1$any_candidate],
      candidate_r2$pair_id[candidate_r2$any_candidate]
    )
  )


  candidate_pairs <- candidate_pairs[nzchar(candidate_pairs)]


  if (length(candidate_pairs) == 0L) {

    return(empty_validated_table())
  }


  ids_r1 <- normalize_pair_id(ShortRead::id(reads_r1))
  ids_r2 <- normalize_pair_id(ShortRead::id(reads_r2))


  q_r1_data <- safe_quality_data(reads_r1)

  stopifnot(
    quality_nrow(q_r1_data) == length(reads_r1)
  )


  q_r2_data <- safe_quality_data(reads_r2)

  stopifnot(
    quality_nrow(q_r2_data) == length(reads_r2)
  )


  s_r1 <- ShortRead::sread(reads_r1)
  s_r2 <- ShortRead::sread(reads_r2)


  results <- list()

  counter <- 0L


  r1_candidate_indices <- which(candidate_r1$any_candidate)


  for (i in r1_candidate_indices) {

    pair_id <- candidate_r1$pair_id[i]

    partner_idx <- match(pair_id, ids_r2)


    if (is.na(partner_idx)) {

      next
    }


    for (target_position in TARGET_POSITIONS) {

      v <- validate_read_against_target(

        seq_obj = s_r2[[partner_idx]],

        qual_vec = quality_vector(q_r2_data, partner_idx),

        read_id = sub(
          "\\s.*$",
          "",
          as.character(ShortRead::id(reads_r2)[[partner_idx]])
        ),

        pair_id = pair_id,

        mate = "R2",

        read_index = candidate_r1$read_index[i],

        target_position = target_position,

        prefilter_count_forward = NA_integer_,
        prefilter_count_reverse = NA_integer_,

        evidence_role = "PAIRED_PARTNER"
      )


      if (!is.null(v)) {

        counter <- counter + 1L

        results[[counter]] <- v
      }
    }
  }


  r2_candidate_indices <- which(candidate_r2$any_candidate)


  for (i in r2_candidate_indices) {

    pair_id <- candidate_r2$pair_id[i]

    partner_idx <- match(pair_id, ids_r1)


    if (is.na(partner_idx)) {

      next
    }


    for (target_position in TARGET_POSITIONS) {

      v <- validate_read_against_target(

        seq_obj = s_r1[[partner_idx]],

        qual_vec = quality_vector(q_r1_data, partner_idx),

        read_id = sub(
          "\\s.*$",
          "",
          as.character(ShortRead::id(reads_r1)[[partner_idx]])
        ),

        pair_id = pair_id,

        mate = "R1",

        read_index = candidate_r2$read_index[i],

        target_position = target_position,

        prefilter_count_forward = NA_integer_,
        prefilter_count_reverse = NA_integer_,

        evidence_role = "PAIRED_PARTNER"
      )


      if (!is.null(v)) {

        counter <- counter + 1L

        results[[counter]] <- v
      }
    }
  }


  if (length(results) == 0L) {

    return(empty_validated_table())
  }


  do.call(rbind, results)
}


# ==============================================================================
# 26. PAIR ALIGNMENT CHECK
# ==============================================================================

assert_pair_alignment <- function(reads_r1, reads_r2) {

  n1 <- length(reads_r1)
  n2 <- length(reads_r2)


  if (n1 != n2) {

    stop(
      "R1/R2 chunk length mismatch: ",
      n1,
      " vs ",
      n2,
      ".\n"
    )
  }


  ids1 <- normalize_pair_id(ShortRead::id(reads_r1))
  ids2 <- normalize_pair_id(ShortRead::id(reads_r2))


  mismatch <- which(ids1 != ids2)


  if (length(mismatch) > 0L) {

    i <- mismatch[1]

    stop(
      "\nR1/R2 pair alignment FAILED inside processing chunk.\n",
      "Local read: ",
      i,
      "\n",
      "R1: ",
      ids1[i],
      "\n",
      "R2: ",
      ids2[i],
      "\n"
    )
  }


  TRUE
}


# ==============================================================================
# 27. PROCESS FASTQ DATASET
# ==============================================================================

process_dataset <- function() {

  fq1 <- ShortRead::FastqStreamer(R1_FASTQ, n = CHUNK_SIZE)
  fq2 <- ShortRead::FastqStreamer(R2_FASTQ, n = CHUNK_SIZE)


  on.exit(try(close(fq1), silent = TRUE), add = TRUE)
  on.exit(try(close(fq2), silent = TRUE), add = TRUE)


  chunk_id <- 0L
  global_read_start <- 1L

  total_r1 <- 0L
  total_r2 <- 0L

  total_primary <- 0L
  total_partner <- 0L


  overall_start <- Sys.time()


  repeat {

    if (RUN_MODE == "PILOT" && total_r1 >= MAX_READS) {

      break
    }


    reads_r1 <- ShortRead::yield(fq1)
    reads_r2 <- ShortRead::yield(fq2)


    n1 <- length(reads_r1)
    n2 <- length(reads_r2)


    if (n1 == 0L && n2 == 0L) {

      break
    }


    if (n1 == 0L || n2 == 0L) {

      stop("R1/R2 FASTQ lengths diverged during processing.\n")
    }


    if (
      RUN_MODE == "PILOT" &&
      (total_r1 + n1 > MAX_READS)
    ) {

      keep_n <- MAX_READS - total_r1


      if (keep_n <= 0L) {

        break
      }


      reads_r1 <- reads_r1[seq_len(keep_n)]
      reads_r2 <- reads_r2[seq_len(keep_n)]

      n1 <- keep_n
      n2 <- keep_n
    }


    chunk_id <- chunk_id + 1L

    read_start <- global_read_start
    read_end <- global_read_start + n1 - 1L


    cat(
      "\n------------------------------------------------------------\n",
      "Chunk ",
      chunk_id,
      " | reads ",
      format(read_start, big.mark = ","),
      " - ",
      format(read_end, big.mark = ","),
      "\n",
      "------------------------------------------------------------\n",
      sep = ""
    )


    if (chunk_complete(chunk_id)) {

      cat("Status: COMPLETE in manifest; skipping analysis.\n")


      total_r1 <- total_r1 + n1
      total_r2 <- total_r2 + n2

      global_read_start <- read_end + 1L

      next
    }


    started_at <- Sys.time()


    mark_chunk_started(
      chunk_id,
      read_start,
      read_end,
      started_at
    )


    assert_pair_alignment(reads_r1, reads_r2)


    global_indices <- read_start:read_end


    pre_r1 <- prefilter_mate(reads_r1, "R1", global_indices)
    pre_r2 <- prefilter_mate(reads_r2, "R2", global_indices)


    candidate_table <- rbind(pre_r1, pre_r2)


    candidate_output <- file.path(
      DIRS$candidates,
      sprintf("chunk_%06d.csv", chunk_id)
    )


    atomic_write_csv(candidate_table, candidate_output)


    cat(
      "R1 nucleotide candidates: ",
      sum(pre_r1$any_candidate),
      "\n",
      "R2 nucleotide candidates: ",
      sum(pre_r2$any_candidate),
      "\n",
      sep = ""
    )


    validated_r1 <- validate_primary_candidates(reads_r1, pre_r1)
    validated_r2 <- validate_primary_candidates(reads_r2, pre_r2)


    primary_validated <- rbind(validated_r1, validated_r2)


    validated_output <- file.path(
      DIRS$validated,
      sprintf("chunk_%06d.csv", chunk_id)
    )


    atomic_write_csv(primary_validated, validated_output)


    paired_validated <- validate_paired_partners(
      reads_r1,
      reads_r2,
      pre_r1,
      pre_r2
    )


    paired_output <- file.path(
      DIRS$paired,
      sprintf("chunk_%06d.csv", chunk_id)
    )


    atomic_write_csv(paired_validated, paired_output)


    elapsed <- as.numeric(
      difftime(Sys.time(), started_at, units = "secs")
    )


    mark_chunk_complete(
      chunk_id,
      list(

        r1_n_reads = as.integer(n1),
        r2_n_reads = as.integer(n2),

        r1_prefilter_candidates = as.integer(sum(pre_r1$any_candidate)),
        r2_prefilter_candidates = as.integer(sum(pre_r2$any_candidate)),

        r1_primary_validated = as.integer(nrow(validated_r1)),
        r2_primary_validated = as.integer(nrow(validated_r2)),

        paired_partner_validated = as.integer(nrow(paired_validated)),

        elapsed_seconds = as.character(elapsed),

        candidate_output = normalizePath(
          candidate_output,
          winslash = "/",
          mustWork = TRUE
        ),

        validated_output = normalizePath(
          validated_output,
          winslash = "/",
          mustWork = TRUE
        ),

        paired_output = normalizePath(
          paired_output,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )


    total_r1 <- total_r1 + n1
    total_r2 <- total_r2 + n2


    total_primary <- total_primary + nrow(primary_validated)
    total_partner <- total_partner + nrow(paired_validated)


    global_read_start <- read_end + 1L


    cat(
      "Chunk complete.\n",
      "  elapsed: ",
      sprintf("%.2f", elapsed),
      " sec\n",
      "  R1 validated: ",
      nrow(validated_r1),
      "\n",
      "  R2 validated: ",
      nrow(validated_r2),
      "\n",
      "  paired-partner validations: ",
      nrow(paired_validated),
      "\n",
      sep = ""
    )


    if (RUN_MODE == "PILOT" && total_r1 >= MAX_READS) {

      break
    }
  }


  overall_elapsed <- as.numeric(
    difftime(Sys.time(), overall_start, units = "secs")
  )


  list(

    chunks = chunk_id,

    r1_reads = total_r1,
    r2_reads = total_r2,

    primary_validated = total_primary,
    partner_validated = total_partner,

    elapsed_seconds = overall_elapsed,

    reads_per_second = if (overall_elapsed > 0) {

      (total_r1 + total_r2) / overall_elapsed

    } else {

      NA_real_
    }
  )
}


# ==============================================================================
# 28. RUN DATASET
# ==============================================================================

cat(
  "\n",
  "====================================================================\n",
  " STARTING FAST NUCLEOTIDE + RIGOROUS POSITIONAL SCAN\n",
  "====================================================================\n",
  "\n",
  sep = ""
)


RUN_STATS <- process_dataset()


cat(
  "\n",
  "====================================================================\n",
  " FAST SCAN COMPLETE\n",
  "====================================================================\n",
  "\n",
  "Chunks processed: ",
  RUN_STATS$chunks,
  "\n",
  "R1 reads: ",
  format(RUN_STATS$r1_reads, big.mark = ","),
  "\n",
  "R2 reads: ",
  format(RUN_STATS$r2_reads, big.mark = ","),
  "\n",
  "Primary validated records: ",
  RUN_STATS$primary_validated,
  "\n",
  "Paired-partner validated records: ",
  RUN_STATS$partner_validated,
  "\n",
  "Elapsed: ",
  sprintf("%.2f", RUN_STATS$elapsed_seconds),
  " sec\n",
  "Combined read rate: ",
  sprintf("%.2f", RUN_STATS$reads_per_second),
  " reads/sec\n",
  "\n",
  sep = ""
)


# ==============================================================================
# 29. READ ALL COMPLETED VALIDATION CHUNKS
# ==============================================================================

read_chunk_files <- function(directory) {

  files <- list.files(
    directory,
    pattern = "^chunk_[0-9]+\\.csv$",
    full.names = TRUE
  )


  if (length(files) == 0L) {

    return(NULL)
  }


  files <- files[order(files)]


  tables <- lapply(
    files,
    function(f) {

      read.csv(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  )


  if (length(tables) == 0L) {

    return(NULL)
  }


  do.call(rbind, tables)
}


ALL_VALIDATED <- read_chunk_files(DIRS$validated)


if (is.null(ALL_VALIDATED)) {

  ALL_VALIDATED <- empty_validated_table()
}


ALL_PARTNER_VALIDATED <- read_chunk_files(DIRS$paired)


if (is.null(ALL_PARTNER_VALIDATED)) {

  ALL_PARTNER_VALIDATED <- empty_validated_table()
}


atomic_write_csv(ALL_VALIDATED, PRIMARY_OUTPUT_ALL)
atomic_write_csv(ALL_PARTNER_VALIDATED, PARTNER_OUTPUT_ALL)


# ==============================================================================
# 30. DEDUPLICATE READ-LEVEL EVIDENCE
# ==============================================================================

deduplicate_evidence <- function(x) {

  if (nrow(x) == 0L) {

    return(x)
  }


  x$role_rank <- ifelse(x$evidence_role == "PRIMARY", 1L, 2L)


  x <- x[
    order(
      x$read_id,
      x$mate,
      x$target_position,
      x$orientation,
      x$frame,
      x$role_rank,
      -x$anchor_identity,
      -x$codon_min_qual
    ),
    ,
    drop = FALSE
  ]


  key <- paste(
    x$read_id,
    x$mate,
    x$target_position,
    x$orientation,
    x$frame,
    sep = "|"
  )


  keep <- !duplicated(key)


  x <- x[keep, , drop = FALSE]

  x$role_rank <- NULL

  rownames(x) <- NULL


  x
}


ALL_EVIDENCE <- deduplicate_evidence(
  rbind(ALL_VALIDATED, ALL_PARTNER_VALIDATED)
)


FINAL_EVIDENCE_PATH <- file.path(
  DIRS$summaries,
  "tac102_final_deduplicated_evidence.csv"
)


atomic_write_csv(ALL_EVIDENCE, FINAL_EVIDENCE_PATH)


# ==============================================================================
# 31. POSITION SUMMARY
# ==============================================================================

make_position_summary <- function(evidence) {

  out <- list()


  for (p in TARGET_POSITIONS) {

    x <- evidence[
      evidence$target_position == p,
      ,
      drop = FALSE
    ]


    if (nrow(x) == 0L) {

      out[[length(out) + 1L]] <- data.frame(

        target_position = p,

        total_hits = 0L,
        unique_reads = 0L,

        M = 0L,
        L = 0L,
        P = 0L,
        other = 0L,

        reference_aa = substr(tac_ref_chr, p, p),

        reference_matches = 0L,

        pct_reference_aa = NA_real_,

        mean_anchor_identity = NA_real_,
        median_anchor_identity = NA_real_,

        mean_codon_min_qual = NA_real_,
        median_codon_min_qual = NA_real_,

        stringsAsFactors = FALSE
      )

      next
    }


    allele_counts <- table(x$target_aa)


    get_count <- function(aa) {

      if (aa %in% names(allele_counts)) {

        as.integer(allele_counts[[aa]])

      } else {

        0L
      }
    }


    reference_aa <- substr(tac_ref_chr, p, p)

    reference_matches <- sum(x$target_aa == reference_aa)


    out[[length(out) + 1L]] <- data.frame(

      target_position = p,

      total_hits = nrow(x),

      unique_reads = length(unique(x$read_id)),

      M = get_count("M"),
      L = get_count("L"),
      P = get_count("P"),

      other = sum(!x$target_aa %in% c("M", "L", "P")),

      reference_aa = reference_aa,

      reference_matches = reference_matches,

      pct_reference_aa = 100 * reference_matches / nrow(x),

      mean_anchor_identity = mean(x$anchor_identity, na.rm = TRUE),
      median_anchor_identity = median(x$anchor_identity, na.rm = TRUE),

      mean_codon_min_qual = mean(x$codon_min_qual, na.rm = TRUE),
      median_codon_min_qual = median(x$codon_min_qual, na.rm = TRUE),

      stringsAsFactors = FALSE
    )
  }


  do.call(rbind, out)
}


POSITION_SUMMARY <- make_position_summary(ALL_EVIDENCE)

atomic_write_csv(POSITION_SUMMARY, POSITION_SUMMARY_PATH)


# ==============================================================================
# 32. ALLELE SUMMARY
# ==============================================================================

if (nrow(ALL_EVIDENCE) > 0L) {

  ALLELE_SUMMARY <- as.data.frame(
    table(
      ALL_EVIDENCE$mate,
      ALL_EVIDENCE$target_position,
      ALL_EVIDENCE$target_aa
    ),
    stringsAsFactors = FALSE
  )


  names(ALLELE_SUMMARY) <- c(
    "mate",
    "target_position",
    "target_aa",
    "hit_count"
  )


  ALLELE_SUMMARY$percent <- NA_real_


  for (p in TARGET_POSITIONS) {

    total <- sum(ALL_EVIDENCE$target_position == p)

    idx <- ALLELE_SUMMARY$target_position == p


    if (total > 0L) {

      ALLELE_SUMMARY$percent[idx] <-
        100 * ALLELE_SUMMARY$hit_count[idx] / total
    }
  }

} else {

  ALLELE_SUMMARY <- data.frame(

    mate = character(),
    target_position = integer(),
    target_aa = character(),
    hit_count = integer(),
    percent = numeric(),

    stringsAsFactors = FALSE
  )
}


atomic_write_csv(ALLELE_SUMMARY, ALLELE_SUMMARY_PATH)


# ==============================================================================
# 33. SINGLE-READ HAPLOTYPE SUMMARY
# ==============================================================================

if (nrow(ALL_EVIDENCE) > 0L) {

  read_position <- unique(
    ALL_EVIDENCE[
      ,
      c(
        "read_id",
        "mate",
        "pair_id",
        "target_position",
        "target_aa",
        "target_codon"
      ),
      drop = FALSE
    ]
  )


  hap_key <- paste(
    read_position$read_id,
    read_position$mate,
    sep = "|"
  )


  split_hap <- split(read_position, hap_key)


  hap_records <- list()

  hcounter <- 0L


  for (z in split_hap) {

    if (all(TARGET_POSITIONS %in% z$target_position)) {

      a653 <- z$target_aa[z$target_position == 653L][1]
      a698 <- z$target_aa[z$target_position == 698L][1]

      c653 <- z$target_codon[z$target_position == 653L][1]
      c698 <- z$target_codon[z$target_position == 698L][1]


      hcounter <- hcounter + 1L


      hap_records[[hcounter]] <- data.frame(

        read_id = z$read_id[1],
        mate = z$mate[1],
        pair_id = z$pair_id[1],

        haplotype_aa = paste0(a653, "653-", a698, "698"),
        haplotype_codon = paste(c653, c698, sep = "-"),

        stringsAsFactors = FALSE
      )
    }
  }


  if (length(hap_records) > 0L) {

    HAPLOTYPES <- do.call(rbind, hap_records)


    HAP_SUMMARY <- aggregate(
      read_id ~ haplotype_aa + haplotype_codon,
      data = HAPLOTYPES,
      FUN = length
    )


    names(HAP_SUMMARY)[names(HAP_SUMMARY) == "read_id"] <- "read_count"


    HAP_SUMMARY$percent <- 100 *
      HAP_SUMMARY$read_count /
      nrow(HAPLOTYPES)

  } else {

    HAP_SUMMARY <- data.frame(

      haplotype_aa = character(),
      haplotype_codon = character(),
      read_count = integer(),
      percent = numeric(),

      stringsAsFactors = FALSE
    )
  }

} else {

  HAP_SUMMARY <- data.frame(

    haplotype_aa = character(),
    haplotype_codon = character(),
    read_count = integer(),
    percent = numeric(),

    stringsAsFactors = FALSE
  )
}


atomic_write_csv(HAP_SUMMARY, HAPLOTYPE_SUMMARY_PATH)


# ==============================================================================
# 34. PAIRED-READ SUMMARY
# ==============================================================================

make_paired_summary <- function(evidence) {

  if (nrow(evidence) == 0L) {

    return(
      data.frame(

        pair_id = character(),
        r1_653 = character(),
        r2_653 = character(),
        r1_698 = character(),
        r2_698 = character(),
        paired_653_status = character(),
        paired_698_status = character(),

        stringsAsFactors = FALSE
      )
    )
  }


  pair_ids <- unique(evidence$pair_id)

  output <- list()


  for (pair_id in pair_ids) {

    x <- evidence[
      evidence$pair_id == pair_id,
      ,
      drop = FALSE
    ]


    get_allele <- function(mate, position) {

      y <- x[
        x$mate == mate &
          x$target_position == position,
        ,
        drop = FALSE
      ]


      if (nrow(y) == 0L) {

        return(NA_character_)
      }


      paste(unique(y$target_aa), collapse = "/")
    }


    r1_653 <- get_allele("R1", 653L)
    r2_653 <- get_allele("R2", 653L)

    r1_698 <- get_allele("R1", 698L)
    r2_698 <- get_allele("R2", 698L)


    status_function <- function(a, b) {

      if (is.na(a) && is.na(b)) {

        return("NO_POSITIONAL_EVIDENCE")
      }


      if (!is.na(a) && !is.na(b)) {

        if (identical(a, b)) {

          return("BOTH_MATES_CONCORDANT")

        } else {

          return("BOTH_MATES_DISCORDANT")
        }
      }


      if (!is.na(a) && is.na(b)) {

        return("R1_ONLY_POSITIONAL_EVIDENCE")
      }


      "R2_ONLY_POSITIONAL_EVIDENCE"
    }


    output[[length(output) + 1L]] <- data.frame(

      pair_id = pair_id,

      r1_653 = r1_653,
      r2_653 = r2_653,

      r1_698 = r1_698,
      r2_698 = r2_698,

      paired_653_status = status_function(r1_653, r2_653),
      paired_698_status = status_function(r1_698, r2_698),

      stringsAsFactors = FALSE
    )
  }


  do.call(rbind, output)
}


PAIRED_SUMMARY <- make_paired_summary(ALL_EVIDENCE)

atomic_write_csv(PAIRED_SUMMARY, PAIRED_SUMMARY_PATH)


# ==============================================================================
# 35. LEGACY PILOT COMPARISON
# ==============================================================================

LEGACY_EVIDENCE_PATH <- file.path(
  PATHS$reports,
  "tac102_4B3p_all_read_positional_evidence.csv"
)


compare_with_legacy <- function(new_evidence, legacy_path) {

  if (!file.exists(legacy_path)) {

    return(
      data.frame(
        status = "LEGACY_FILE_NOT_FOUND",
        stringsAsFactors = FALSE
      )
    )
  }


  old <- read.csv(
    legacy_path,
    stringsAsFactors = FALSE
  )


  if (nrow(old) == 0L) {

    return(
      data.frame(
        status = "LEGACY_FILE_EMPTY",
        stringsAsFactors = FALSE
      )
    )
  }


  old_key <- unique(
    paste(old$read_id, old$mate, old$target_position, sep = "|")
  )


  new_key <- unique(
    paste(
      new_evidence$read_id,
      new_evidence$mate,
      new_evidence$target_position,
      sep = "|"
    )
  )


  overlap <- intersect(old_key, new_key)

  missing_from_new <- setdiff(old_key, new_key)

  new_only <- setdiff(new_key, old_key)


  data.frame(

    status = "COMPARED",

    legacy_unique_read_position_records = length(old_key),

    fast_unique_read_position_records = length(new_key),

    overlap_records = length(overlap),

    legacy_records_not_recovered_by_fast = length(missing_from_new),

    fast_records_not_in_legacy = length(new_only),

    legacy_recovery_percent = if (length(old_key) > 0) {

      100 * length(overlap) / length(old_key)

    } else {

      NA_real_
    },

    stringsAsFactors = FALSE
  )
}


PILOT_COMPARISON <- compare_with_legacy(
  ALL_EVIDENCE,
  LEGACY_EVIDENCE_PATH
)


atomic_write_csv(PILOT_COMPARISON, PILOT_COMPARISON_PATH)


# ==============================================================================
# 36. FINAL REPORT
# ==============================================================================

report_lines <- c(

  "====================================================================",
  "TAC102 FAST ALL-READ POSITIONAL VALIDATION",
  "Step 04p2",
  "====================================================================",
  "",
  paste0("Date: ", Sys.time()),
  paste0("Run mode: ", RUN_MODE),
  paste0(
    "Maximum reads per FASTQ: ",
    ifelse(
      is.infinite(MAX_READS),
      "FULL DATASET",
      format(MAX_READS, big.mark = ",")
    )
  ),
  paste0("Chunk size: ", format(CHUNK_SIZE, big.mark = ",")),
  "",
  "REFERENCE",
  paste0("  TAC102 reference length: ", TAC_REF_LENGTH, " aa"),
  paste0("  Position 653 reference: ", reference_653),
  paste0("  Position 698 reference: ", reference_698),
  "",
  "TARGETS",
  paste0("  Position 653 anchor: ", ANCHORS[["653"]]),
  paste0("  Position 698 anchor: ", ANCHORS[["698"]]),
  "",
  "PREFILTER",
  "  Type: protein-derived codon-degenerate nucleotide pattern",
  "  Target codon represented as NNN",
  paste0("  Maximum nucleotide mismatches: ", PREFILTER_MAX_MISMATCH),
  "  Both forward and reverse-complement patterns searched",
  "",
  "VALIDATOR",
  paste0("  Minimum anchor identity: ", MIN_ANCHOR_IDENTITY),
  paste0(
    "  Maximum AA mismatches during anchor localization: ",
    VALIDATOR_MAX_AA_MISMATCH
  ),
  "  Six frames evaluated only for nucleotide-prefilter candidates",
  "  Internal codons translated with no.init.codon=TRUE",
  "  Ambiguous target codons rejected",
  "",
  "PRE-FLIGHT",
  paste0(
    "  Synthetic tests passed: ",
    sum(PREFLIGHT$passed),
    "/",
    nrow(PREFLIGHT)
  ),
  paste0(
    "  FASTQ pair preflight: ",
    ifelse(PAIR_PREFLIGHT$passed, "PASSED", "FAILED")
  ),
  "",
  "PROCESSING",
  paste0("  Chunks: ", RUN_STATS$chunks),
  paste0("  R1 reads: ", format(RUN_STATS$r1_reads, big.mark = ",")),
  paste0("  R2 reads: ", format(RUN_STATS$r2_reads, big.mark = ",")),
  paste0("  Primary validated records: ", RUN_STATS$primary_validated),
  paste0("  Paired-partner validated records: ", RUN_STATS$partner_validated),
  paste0("  Elapsed seconds: ", sprintf("%.2f", RUN_STATS$elapsed_seconds)),
  paste0(
    "  Combined read rate: ",
    sprintf("%.2f", RUN_STATS$reads_per_second),
    " reads/sec"
  ),
  "",
  "SCIENTIFIC INTERPRETATION",
  "  The script identifies reads that can be independently anchored",
  "  to TAC102 positions 653 and 698.",
  "",
  "  The current reference is:",
  "      K653",
  "      L698",
  "",
  "  A read carrying M at position 653 is therefore classified as:",
  "      K653M",
  "",
  "  A read carrying L at position 698 is classified as:",
  "      L698 reference match",
  "",
  "  This analysis does NOT establish genomic fixation, population",
  "  allele frequency, or phenotype by itself.",
  "",
  "  Paired-read support should be evaluated separately from primary",
  "  positional observations.",
  "",
  "OUTPUTS",
  paste0("  Pattern catalogue: ", PATTERN_CATALOGUE_PATH),
  paste0("  Preflight: ", PREFLIGHT_PATH),
  paste0("  Chunk manifest: ", MANIFEST_PATH),
  paste0("  Candidate chunks: ", DIRS$candidates),
  paste0("  Validated chunks: ", DIRS$validated),
  paste0("  Paired chunks: ", DIRS$paired),
  paste0("  Final evidence: ", FINAL_EVIDENCE_PATH),
  paste0("  Position summary: ", POSITION_SUMMARY_PATH),
  paste0("  Allele summary: ", ALLELE_SUMMARY_PATH),
  paste0("  Single-read haplotypes: ", HAPLOTYPE_SUMMARY_PATH),
  paste0("  Paired-read summary: ", PAIRED_SUMMARY_PATH),
  paste0("  Legacy pilot comparison: ", PILOT_COMPARISON_PATH),
  paste0("  Report: ", REPORT_PATH),
  "",
  "===================================================================="
)


writeLines(report_lines, REPORT_PATH)


# ==============================================================================
# 37. FINAL LOGGING
# ==============================================================================

log_fast(
  paste0(
    "04p2 complete. R1=",
    RUN_STATS$r1_reads,
    "; R2=",
    RUN_STATS$r2_reads,
    "; primary=",
    RUN_STATS$primary_validated,
    "; paired=",
    RUN_STATS$partner_validated
  )
)


cat(
  "\n",
  "====================================================================\n",
  " 04p2 COMPLETE\n",
  "====================================================================\n",
  "\n",
  "Final evidence:\n",
  "  ",
  FINAL_EVIDENCE_PATH,
  "\n\n",
  "Position summary:\n",
  "  ",
  POSITION_SUMMARY_PATH,
  "\n\n",
  "Paired-read summary:\n",
  "  ",
  PAIRED_SUMMARY_PATH,
  "\n\n",
  "Report:\n",
  "  ",
  REPORT_PATH,
  "\n\n",
  sep = ""
)