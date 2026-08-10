# ==============================================================================
# Step 4B.3c
# TAC102 Orientation / Coordinate Diagnostic
#
# PURPOSE
# -------
# Independently reconstruct the nucleotide sequence corresponding to each
# candidate's TAC102 gap coordinates.
#
# This script DOES NOT reuse the 4B.3b extraction functions.
#
# For every candidate we test:
#   1. raw segment translated as-is
#   2. reverse-complemented raw segment translated
#
# Each is compared against:
#   A. forward TAC102 reference
#   B. reversed-reference control
#
# The purpose is to determine whether apparent TAC102 correspondence is
# biologically oriented correctly or could be produced by coordinate/reversal
# handling errors.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. LOAD SETUP
# ------------------------------------------------------------------------------

if (!exists("PATHS") || !exists("CONFIG")) {
    source(
        here::here(
            "raw_read_analysis",
            "scripts",
            "00_setup.R"
        )
    )
}


# ------------------------------------------------------------------------------
# 1. INPUTS
# ------------------------------------------------------------------------------

tails_file <- file.path(
    PATHS$reports,
    "tac102_41pairs_unaligned_tails.csv"
)

reference_file <- file.path(
    PATHS$intermediate,
    "reference_sequences.rds"
)

r1_file <- file.path(
    PATHS$raw_recovery_derivative,
    "tac102_41pair_R1.fasta"
)

r2_file <- file.path(
    PATHS$raw_recovery_derivative,
    "tac102_41pair_R2.fasta"
)


tails <- read.csv(
    tails_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

reference_sequences <- readRDS(reference_file)

tac_ref_chr <- as.character(
    reference_sequences[["TAC102"]]
)

r1_seqs <- Biostrings::readDNAStringSet(r1_file)
r2_seqs <- Biostrings::readDNAStringSet(r2_file)

# Normalize FASTQ names to pair IDs
names(r1_seqs) <- sub("\\s.*$", "", names(r1_seqs))
names(r2_seqs) <- sub("\\s.*$", "", names(r2_seqs))


# ------------------------------------------------------------------------------
# 2. BASIC CHECKS
# ------------------------------------------------------------------------------

cat("\n==============================================================\n")
cat("TAC102 ORIENTATION / COORDINATE DIAGNOSTIC\n")
cat("Step 4B.3c\n")
cat("==============================================================\n\n")

cat("Candidate rows:", nrow(tails), "\n")
cat("TAC102 reference length:", nchar(tac_ref_chr), "aa\n")
cat("R1 recovered reads:", length(r1_seqs), "\n")
cat("R2 recovered reads:", length(r2_seqs), "\n\n")


# ------------------------------------------------------------------------------
# 3. SAFE CHARACTER EXTRACTION
# ------------------------------------------------------------------------------

safe_chr <- function(x) {

    if (length(x) == 0 || is.na(x[1])) {
        return(NA_character_)
    }

    as.character(x[1])
}


safe_int <- function(x) {

    if (length(x) == 0 || is.na(x[1])) {
        return(NA_integer_)
    }

    as.integer(x[1])
}


# ------------------------------------------------------------------------------
# 4. TRANSLATION FUNCTION
# ------------------------------------------------------------------------------

translate_dna_safe <- function(dna) {

    if (is.null(dna)) {
        return(NA_character_)
    }

    dna <- as.character(dna)

    if (!nzchar(dna)) {
        return(NA_character_)
    }

    usable_length <- nchar(dna) -
        (nchar(dna) %% 3L)

    if (usable_length < 3L) {
        return(NA_character_)
    }

    dna <- substr(
        dna,
        1L,
        usable_length
    )

    result <- tryCatch(
        Biostrings::translate(
            Biostrings::DNAString(dna),
            if.fuzzy.codon = "X"
        ),
        error = function(e) NULL
    )

    if (is.null(result)) {
        return(NA_character_)
    }

    as.character(result)
}


# ------------------------------------------------------------------------------
# 5. IDENTITY FUNCTION
# ------------------------------------------------------------------------------

protein_identity <- function(
    candidate,
    reference
) {

    if (
        is.na(candidate) ||
        is.na(reference) ||
        !nzchar(candidate) ||
        !nzchar(reference)
    ) {
        return(NA_real_)
    }

    candidate_chars <- strsplit(
        candidate,
        "",
        fixed = TRUE
    )[[1]]

    reference_chars <- strsplit(
        reference,
        "",
        fixed = TRUE
    )[[1]]

    n <- min(
        length(candidate_chars),
        length(reference_chars)
    )

    if (n == 0L) {
        return(NA_real_)
    }

    mean(
        candidate_chars[seq_len(n)] ==
            reference_chars[seq_len(n)]
    ) * 100
}


# ------------------------------------------------------------------------------
# 6. REFERENCE REVERSE CONTROL
#
# IMPORTANT:
# This is deliberately only a diagnostic control.
# It is NOT a biologically valid alternative TAC102 orientation.
# ------------------------------------------------------------------------------

reverse_reference <- function(x) {

    if (
        is.na(x) ||
        !nzchar(x)
    ) {
        return(NA_character_)
    }

    paste(
        rev(strsplit(x, "", fixed = TRUE)[[1]]),
        collapse = ""
    )
}


# ------------------------------------------------------------------------------
# 7. ONE-ROW EMPTY RESULT
#
# This is deliberately explicit so that every error still produces exactly
# one row with exactly the same columns as a successful result.
# ------------------------------------------------------------------------------

empty_result <- function(
    row,
    note
) {

    data.frame(

        pair_id = safe_chr(row$pair_id),

        mate = safe_chr(row$mate),

        evidence_class = safe_chr(row$evidence_class),

        hsp_strand_reported = safe_chr(
            row$hsp_strand
        ),

        hsp_sstart = safe_int(
            row$hsp_sstart
        ),

        hsp_send = safe_int(
            row$hsp_send
        ),

        hsp_qstart = safe_int(
            row$hsp_qstart
        ),

        hsp_qend = safe_int(
            row$hsp_qend
        ),

        gap_position_start = safe_int(
            row$gap_position_start
        ),

        gap_position_end = safe_int(
            row$gap_position_end
        ),

        reconstructed_forward = NA,

        nt_lo = NA_integer_,

        nt_hi = NA_integer_,

        nt_window_length = NA_integer_,

        raw_segment = NA_character_,

        translated_as_is = NA_character_,

        translated_reverse_complement = NA_character_,

        reference_forward = NA_character_,

        reference_reversed_control = NA_character_,

        id_as_is_vs_reference = NA_real_,

        id_reverse_complement_vs_reference = NA_real_,

        id_as_is_vs_reversed_control = NA_real_,

        id_reverse_complement_vs_reversed_control = NA_real_,

        best_biological_orientation = NA_character_,

        candidate_reference_length_match = NA,

        candidate_reference_length_difference = NA_integer_,

        note = as.character(note),

        stringsAsFactors = FALSE,

        check.names = FALSE
    )
}


# ------------------------------------------------------------------------------
# 8. DIAGNOSTIC FOR ONE CANDIDATE
# ------------------------------------------------------------------------------

diagnose_row <- function(row) {

    tryCatch({

        # ----------------------------------------------------------------------
        # Coordinates
        # ----------------------------------------------------------------------

        s_raw <- safe_int(row$hsp_sstart)
        e_raw <- safe_int(row$hsp_send)

        q_start <- safe_int(row$hsp_qstart)
        q_end <- safe_int(row$hsp_qend)

        p_start <- safe_int(row$gap_position_start)
        p_end <- safe_int(row$gap_position_end)

        if (
            any(
                is.na(
                    c(
                        s_raw,
                        e_raw,
                        q_start,
                        q_end,
                        p_start,
                        p_end
                    )
                )
            )
        ) {

            return(
                empty_result(
                    row,
                    "insufficient coordinates"
                )
            )
        }


        # ----------------------------------------------------------------------
        # Determine orientation directly from raw subject coordinates
        # ----------------------------------------------------------------------

        forward <- s_raw <= e_raw

        strand <- if (forward) {
            "forward"
        } else {
            "reverse"
        }


        # ----------------------------------------------------------------------
        # Recover the original read
        # ----------------------------------------------------------------------

        pair_id <- safe_chr(row$pair_id)

        mate <- safe_chr(row$mate)

        if (mate == "R1") {

            read_seq <- r1_seqs[[pair_id]]

        } else if (mate == "R2") {

            read_seq <- r2_seqs[[pair_id]]

        } else {

            return(
                empty_result(
                    row,
                    paste(
                        "unknown mate:",
                        mate
                    )
                )
            )
        }


        if (is.null(read_seq)) {

            return(
                empty_result(
                    row,
                    "read not found"
                )
            )
        }


        # ----------------------------------------------------------------------
        # IMPORTANT:
        # Use length(), not width(), for DNAString.
        # ----------------------------------------------------------------------

        read_len <- length(read_seq)


        # ----------------------------------------------------------------------
        # Reconstruct nucleotide window
        #
        # Protein position p_start corresponds to q_start in the HSP.
        #
        # Forward:
        #   sstart + (p - qstart)*3
        #
        # Reverse:
        #   send + (qstart - p)*3
        #
        # The raw nucleotide interval is then normalized to low/high genomic
        # coordinates in the read.
        # ----------------------------------------------------------------------

        if (forward) {

            nt_first <- s_raw +
                (p_start - q_start) * 3L

            nt_last <- s_raw +
                (p_end - q_start) * 3L + 2L

        } else {

            nt_first <- e_raw -
                (p_start - q_start) * 3L

            nt_last <- e_raw -
                (p_end - q_start) * 3L - 2L
        }


        nt_lo <- min(
            nt_first,
            nt_last
        )

        nt_hi <- max(
            nt_first,
            nt_last
        )


        # ----------------------------------------------------------------------
        # Bounds check
        # ----------------------------------------------------------------------

        if (
            nt_lo < 1L ||
            nt_hi > read_len ||
            nt_lo > nt_hi
        ) {

            return(
                empty_result(
                    row,
                    paste0(
                        "reconstructed nucleotide window ",
                        nt_lo,
                        ":",
                        nt_hi,
                        " outside read length ",
                        read_len
                    )
                )
            )
        }


        # ----------------------------------------------------------------------
        # Extract raw nucleotide segment
        # ----------------------------------------------------------------------

        raw_segment <- Biostrings::subseq(
            read_seq,
            start = nt_lo,
            end = nt_hi
        )


        raw_segment_chr <- as.character(
            raw_segment
        )


        nt_window_length <- nchar(
            raw_segment_chr
        )


        # ----------------------------------------------------------------------
        # Translate both orientations
        # ----------------------------------------------------------------------

        translated_as_is <- translate_dna_safe(
            raw_segment_chr
        )


        translated_reverse_complement <- translate_dna_safe(
            as.character(
                Biostrings::reverseComplement(
                    raw_segment
                )
            )
        )


        # ----------------------------------------------------------------------
        # Reference slice
        # ----------------------------------------------------------------------

        reference_forward <- substr(
            tac_ref_chr,
            p_start,
            p_end
        )


        reference_reversed_control <- reverse_reference(
            reference_forward
        )


        # ----------------------------------------------------------------------
        # Identity calculations
        #
        # We calculate identity over the overlapping length rather than
        # requiring exact length equality. Length differences are recorded
        # separately.
        # ----------------------------------------------------------------------

        id_as_is_vs_reference <- protein_identity(
            translated_as_is,
            reference_forward
        )


        id_reverse_complement_vs_reference <- protein_identity(
            translated_reverse_complement,
            reference_forward
        )


        id_as_is_vs_reversed_control <- protein_identity(
            translated_as_is,
            reference_reversed_control
        )


        id_reverse_complement_vs_reversed_control <- protein_identity(
            translated_reverse_complement,
            reference_reversed_control
        )


        # ----------------------------------------------------------------------
        # Determine best biologically meaningful orientation
        #
        # Only comparison against the FORWARD reference is biologically
        # meaningful.
        # ----------------------------------------------------------------------

        biological_scores <- c(
            as_is = id_as_is_vs_reference,
            reverse_complement =
                id_reverse_complement_vs_reference
        )


        if (all(is.na(biological_scores))) {

            best_orientation <- NA_character_

        } else {

            best_orientation <- names(
                which.max(
                    biological_scores
                )
            )
        }


        # ----------------------------------------------------------------------
        # Length diagnostics
        # ----------------------------------------------------------------------

        reference_length <- nchar(
            reference_forward
        )

        candidate_lengths <- c(
            as_is =
                ifelse(
                    is.na(translated_as_is),
                    NA_integer_,
                    nchar(translated_as_is)
                ),

            reverse_complement =
                ifelse(
                    is.na(translated_reverse_complement),
                    NA_integer_,
                    nchar(translated_reverse_complement)
                )
        )


        best_candidate_length <- if (
            all(is.na(biological_scores))
        ) {
            NA_integer_
        } else {

            candidate_lengths[
                names(
                    which.max(
                        biological_scores
                    )
                )
            ]
        }


        length_difference <- if (
            is.na(best_candidate_length)
        ) {
            NA_integer_
        } else {

            best_candidate_length -
                reference_length
        }


        length_match <- if (
            is.na(length_difference)
        ) {
            NA
        } else {

            length_difference == 0L
        }


        # ----------------------------------------------------------------------
        # Final result
        # ----------------------------------------------------------------------

        data.frame(

            pair_id = pair_id,

            mate = mate,

            evidence_class =
                safe_chr(row$evidence_class),

            hsp_strand_reported =
                safe_chr(row$hsp_strand),

            hsp_sstart = s_raw,

            hsp_send = e_raw,

            hsp_qstart = q_start,

            hsp_qend = q_end,

            gap_position_start = p_start,

            gap_position_end = p_end,

            reconstructed_forward =
                strand == "forward",

            nt_lo = nt_lo,

            nt_hi = nt_hi,

            nt_window_length =
                nt_window_length,

            raw_segment =
                raw_segment_chr,

            translated_as_is =
                translated_as_is,

            translated_reverse_complement =
                translated_reverse_complement,

            reference_forward =
                reference_forward,

            reference_reversed_control =
                reference_reversed_control,

            id_as_is_vs_reference =
                id_as_is_vs_reference,

            id_reverse_complement_vs_reference =
                id_reverse_complement_vs_reference,

            id_as_is_vs_reversed_control =
                id_as_is_vs_reversed_control,

            id_reverse_complement_vs_reversed_control =
                id_reverse_complement_vs_reversed_control,

            best_biological_orientation =
                best_orientation,

            candidate_reference_length_match =
                length_match,

            candidate_reference_length_difference =
                length_difference,

            note =
                paste0(
                    "Independent reconstruction; reported HSP strand = ",
                    strand
                ),

            stringsAsFactors = FALSE,

            check.names = FALSE
        )

    }, error = function(e) {

        empty_result(
            row,
            paste(
                "Diagnostic error:",
                conditionMessage(e)
            )
        )
    })
}


# ------------------------------------------------------------------------------
# 9. RUN ALL 82 CANDIDATES
# ------------------------------------------------------------------------------

diag_list <- lapply(
    seq_len(nrow(tails)),
    function(i) {

        diagnose_row(
            tails[i, , drop = FALSE]
        )
    }
)


# ------------------------------------------------------------------------------
# 10. SAFELY COMBINE RESULTS
# ------------------------------------------------------------------------------

expected_columns <- c(

    "pair_id",
    "mate",
    "evidence_class",
    "hsp_strand_reported",
    "hsp_sstart",
    "hsp_send",
    "hsp_qstart",
    "hsp_qend",
    "gap_position_start",
    "gap_position_end",
    "reconstructed_forward",
    "nt_lo",
    "nt_hi",
    "nt_window_length",
    "raw_segment",
    "translated_as_is",
    "translated_reverse_complement",
    "reference_forward",
    "reference_reversed_control",
    "id_as_is_vs_reference",
    "id_reverse_complement_vs_reference",
    "id_as_is_vs_reversed_control",
    "id_reverse_complement_vs_reversed_control",
    "best_biological_orientation",
    "candidate_reference_length_match",
    "candidate_reference_length_difference",
    "note"
)


# Force every element to have the same columns and order.
diag_list <- lapply(
    diag_list,
    function(x) {

        missing_cols <- setdiff(
            expected_columns,
            names(x)
        )

        for (mc in missing_cols) {
            x[[mc]] <- NA
        }

        x <- x[, expected_columns, drop = FALSE]

        # Enforce exactly one row.
        if (nrow(x) != 1L) {
            stop(
                "Diagnostic result did not return exactly one row."
            )
        }

        x
    }
)


diag_results <- do.call(
    rbind,
    diag_list
)

rownames(diag_results) <- NULL


# ------------------------------------------------------------------------------
# 11. SAVE FULL DIAGNOSTIC TABLE
# ------------------------------------------------------------------------------

output_file <- file.path(
    PATHS$reports,
    "tac102_4B3c_orientation_diagnostic.csv"
)

write.csv(
    diag_results,
    output_file,
    row.names = FALSE,
    na = ""
)


# ------------------------------------------------------------------------------
# 12. SUMMARY
# ------------------------------------------------------------------------------

cat("\n==============================================================\n")
cat("DIAGNOSTIC SUMMARY\n")
cat("==============================================================\n\n")

cat(
    "Candidate rows analysed:",
    nrow(diag_results),
    "\n\n"
)


cat("Reported HSP strand:\n")
print(
    table(
        diag_results$hsp_strand_reported,
        useNA = "ifany"
    )
)


cat("\nDiagnostic notes:\n")
print(
    sort(
        table(diag_results$note),
        decreasing = TRUE
    )
)


cat("\nBest biologically meaningful orientation:\n")
print(
    table(
        diag_results$best_biological_orientation,
        useNA = "ifany"
    )
)


cat("\nReference length match:\n")
print(
    table(
        diag_results$candidate_reference_length_match,
        useNA = "ifany"
    )
)


# ------------------------------------------------------------------------------
# 13. REVERSAL CONTROL
# ------------------------------------------------------------------------------

reversal_excess <- with(
    diag_results,
    pmax(
        id_as_is_vs_reversed_control -
            id_as_is_vs_reference,
        id_reverse_complement_vs_reversed_control -
            id_reverse_complement_vs_reference,
        na.rm = TRUE
    )
)

# If both comparisons are NA, pmax(..., na.rm=TRUE) gives -Inf.
# Convert those to NA.
reversal_excess[is.infinite(reversal_excess)] <- NA_real_


diag_results$reversal_control_excess <-
    reversal_excess


cat("\nRows with >=20 percentage-point reversal-control excess:\n")

print(
    sum(
        diag_results$reversal_control_excess >= 20,
        na.rm = TRUE
    )
)


# ------------------------------------------------------------------------------
# 14. SAVE UPDATED TABLE WITH REVERSAL EXCESS
# ------------------------------------------------------------------------------

write.csv(
    diag_results,
    output_file,
    row.names = FALSE,
    na = ""
)


# ------------------------------------------------------------------------------
# 15. WRITE HUMAN-READABLE SUMMARY
# ------------------------------------------------------------------------------

summary_file <- file.path(
    PATHS$reports,
    "tac102_4B3c_orientation_summary.txt"
)

sink(summary_file)

cat("==============================================================\n")
cat("TAC102 ORIENTATION / COORDINATE DIAGNOSTIC\n")
cat("Step 4B.3c\n")
cat("==============================================================\n\n")

cat(
    "Candidate rows analysed:",
    nrow(diag_results),
    "\n"
)

cat(
    "TAC102 reference length:",
    nchar(tac_ref_chr),
    "aa\n\n"
)

cat("Reported HSP strand:\n")
print(
    table(
        diag_results$hsp_strand_reported,
        useNA = "ifany"
    )
)

cat("\nBest biologically meaningful orientation:\n")
print(
    table(
        diag_results$best_biological_orientation,
        useNA = "ifany"
    )
)

cat("\nReference length match:\n")
print(
    table(
        diag_results$candidate_reference_length_match,
        useNA = "ifany"
    )
)

cat("\nReversal-control diagnostics:\n")

cat(
    "Rows with >=20 percentage-point reversal-control excess:",
    sum(
        diag_results$reversal_control_excess >= 20,
        na.rm = TRUE
    ),
    "\n"
)

cat("\nDiagnostic notes:\n")
print(
    sort(
        table(diag_results$note),
        decreasing = TRUE
    )
)

cat("\nInterpretation rule:\n")
cat(
    "The biologically meaningful comparison is candidate translation\n",
    "versus the forward TAC102 reference.\n",
    "The reversed-reference comparison is diagnostic only.\n",
    "A strong match to the reversed reference is suspicious for\n",
    "sequence-orientation/reversal handling errors.\n"
)

sink()


# ------------------------------------------------------------------------------
# 16. FINAL CONSOLE MESSAGE
# ------------------------------------------------------------------------------

cat("\n==============================================================\n")
cat("4B.3c COMPLETE\n")
cat("==============================================================\n")

cat(
    "Diagnostic CSV:\n",
    output_file,
    "\n\n"
)

cat(
    "Summary:\n",
    summary_file,
    "\n"
)