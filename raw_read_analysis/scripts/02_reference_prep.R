# ============================================================================== 
# Step 2: Reference Integrity Verification & Preparation
# Project: Phase II Validation of Candidate Genes from T. equiperdum IVM-t1 Raw Reads
# File: scripts/02_reference_prep.R
# ============================================================================== 

cat("======================================================\n")
cat(" Phase II: Reference Integrity Verification          \n")
cat("======================================================\n\n")

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("TARGETS")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

log_info("Starting Step 2: Reference Integrity Verification...", LOG_FILE)

if (is.null(PATHS$intermediate)) {
  PATHS$intermediate <- here::here("raw_read_analysis", "results", "intermediate")
}
if (!dir.exists(PATHS$intermediate)) {
  dir.create(PATHS$intermediate, showWarnings = FALSE, recursive = TRUE)
}

VALID_STOP_CODONS <- c("TAA", "TAG", "TGA")
VALID_START_CODONS <- c("ATG", "GTG", "TTG")

# The original workflow used Biostrings::translate() without a genetic.code
# argument, i.e. the standard genetic code. No translation table is changed
# here. Confirm the appropriate code before validating any nucleotide reference.
if (is.null(CONFIG$translation_genetic_code)) {
  log_warn(
    "CONFIG does not define translation_genetic_code; nucleotide validation uses Biostrings' existing default (standard genetic code). Confirm this scientific assumption before use.",
    LOG_FILE
  )
}

empty_result <- function(gene, file, reference_type = NA_character_,
                         header = NA_character_, checksum = NA_character_, reason) {
  list(
    gene = gene, file = file, header = header, checksum = checksum,
    reference_type = reference_type, seq_length_bp = NA_integer_,
    expected_bp = NA_integer_, cds_length_match = NA,
    start_codon = NA_character_, start_codon_valid = NA,
    stop_codon = NA_character_, stop_codon_valid = NA,
    cds_divisible_by_3 = NA, internal_stops = NA_integer_,
    ambiguous_bases = NA_integer_, aa_length_including_terminal_stop = NA_integer_,
    aa_length_excluding_terminal_stop = NA_integer_, expected_aa_len = NA_integer_,
    aa_length_match = NA, terminal_stop_included = NA, gc_percent = NA_real_,
    header_valid = FALSE, header_contains_target_gene = FALSE, header_unique = FALSE, translation_success = FALSE,
    status = "FAIL", reason = reason, translation = NULL
  )
}

append_reason <- function(existing_reason, new_reason) {
  if (identical(existing_reason, "OK")) return(new_reason)
  paste(existing_reason, new_reason, sep = "; ")
}

# Hash only the normalized sequence characters, not FASTA layout or header, so
# the checksum is reproducible after harmless line-wrapping changes.
sequence_checksum <- function(sequence) {
  checksum_file <- tempfile("reference_sequence_", fileext = ".txt")
  on.exit(unlink(checksum_file), add = TRUE)
  writeLines(as.character(sequence), checksum_file, useBytes = TRUE)
  unname(tools::md5sum(checksum_file))
}

# Metadata remains authoritative. Content inspection is used only when metadata
# is absent or unrecognised, preventing a protein reference from being read as DNA.
normalize_reference_type <- function(reference_type) {
  type <- tolower(trimws(as.character(reference_type)))
  if (type %in% c("nucleotide", "dna", "cds")) return("nucleotide")
  if (type %in% c("protein", "amino_acid", "amino acid")) return("protein")
  NA_character_
}

detect_seq_type_from_content <- function(fasta_path) {
  raw_lines <- readLines(fasta_path, warn = FALSE)
  seq_lines <- raw_lines[!grepl("^>", raw_lines)]
  seq_chars <- toupper(gsub("[[:space:]]", "", paste(seq_lines, collapse = "")))
  if (!nzchar(seq_chars)) stop("FASTA contains no sequence characters.")
  nucleotide_codes <- c("A", "C", "G", "T", "U", "N", "R", "Y", "W", "S", "K", "M", "B", "D", "H", "V", "-")
  chars <- unique(strsplit(seq_chars, "", fixed = TRUE)[[1]])
  if (all(chars %in% nucleotide_codes)) "nucleotide" else "protein"
}

validate_nucleotide_ref <- function(gene, dna_set, expected_cds_len, expected_aa_len, log_file) {
  sequence <- dna_set[[1]]
  sequence_chr <- as.character(sequence)
  cds_len <- length(sequence)
  failures <- character()

  cds_length_match <- !is.na(expected_cds_len) && expected_cds_len > 0 && cds_len == expected_cds_len
  if (!cds_length_match) failures <- c(failures, sprintf("CDS length %d bp does not match expected %s bp", cds_len, expected_cds_len))

  cds_divisible_by_3 <- cds_len %% 3 == 0
  if (!cds_divisible_by_3) failures <- c(failures, "CDS length is not divisible by 3")

  start_codon <- if (cds_len >= 3) substr(sequence_chr, 1, 3) else NA_character_
  start_codon_valid <- !is.na(start_codon) && start_codon %in% VALID_START_CODONS
  if (!start_codon_valid) failures <- c(failures, sprintf("invalid start codon: %s", start_codon))

  stop_codon <- if (cds_len >= 3) substr(sequence_chr, cds_len - 2, cds_len) else NA_character_
  stop_codon_valid <- !is.na(stop_codon) && stop_codon %in% VALID_STOP_CODONS
  if (!stop_codon_valid) failures <- c(failures, sprintf("invalid stop codon: %s", stop_codon))

  ambiguous_bases <- sum(!(strsplit(sequence_chr, "", fixed = TRUE)[[1]] %in% c("A", "C", "G", "T")))
  if (ambiguous_bases > 0) failures <- c(failures, sprintf("ambiguous nucleotides: %d", ambiguous_bases))

  translated <- tryCatch(
    Biostrings::translate(sequence, if.fuzzy.codon = "X"),
    error = function(e) {
      log_error(sprintf("  [FAIL] %s: Translation failed - %s", gene, conditionMessage(e)), log_file)
      NULL
    }
  )
  translation_success <- !is.null(translated)
  aa_chr <- if (translation_success) as.character(translated) else NA_character_
  terminal_stop_included <- translation_success && nchar(aa_chr) > 0 && substr(aa_chr, nchar(aa_chr), nchar(aa_chr)) == "*"
  aa_len_including <- if (translation_success) nchar(aa_chr) else NA_integer_
  aa_len_excluding <- if (translation_success) aa_len_including - as.integer(terminal_stop_included) else NA_integer_
  internal_stops <- if (translation_success) {
    body <- if (terminal_stop_included) substr(aa_chr, 1, nchar(aa_chr) - 1) else aa_chr
    sum(strsplit(body, "", fixed = TRUE)[[1]] == "*")
  } else NA_integer_

  if (!translation_success) failures <- c(failures, "translation failed")
  if (translation_success && internal_stops > 0) failures <- c(failures, sprintf("internal stop codons: %d", internal_stops))
  aa_length_match <- translation_success && !is.na(expected_aa_len) && expected_aa_len > 0 && aa_len_excluding == expected_aa_len
  if (!aa_length_match) failures <- c(failures, sprintf("translated protein length %s aa (excluding terminal stop) does not match expected %s aa", aa_len_excluding, expected_aa_len))

  gc_percent <- round(sum(strsplit(sequence_chr, "", fixed = TRUE)[[1]] %in% c("G", "C")) / cds_len * 100, 2)
  list(
    reference_type = "nucleotide", seq_length_bp = cds_len, expected_bp = expected_cds_len,
    cds_length_match = cds_length_match, start_codon = start_codon,
    start_codon_valid = start_codon_valid, stop_codon = stop_codon,
    stop_codon_valid = stop_codon_valid, cds_divisible_by_3 = cds_divisible_by_3,
    internal_stops = internal_stops, ambiguous_bases = ambiguous_bases,
    aa_length_including_terminal_stop = aa_len_including,
    aa_length_excluding_terminal_stop = aa_len_excluding,
    expected_aa_len = expected_aa_len, aa_length_match = aa_length_match,
    terminal_stop_included = terminal_stop_included, gc_percent = gc_percent,
    translation_success = translation_success,
    status = if (length(failures) == 0) "PASS" else "FAIL",
    reason = if (length(failures) == 0) "OK" else paste(failures, collapse = "; "),
    translation = translated
  )
}

validate_protein_ref <- function(aa_set, expected_aa_len) {
  sequence <- aa_set[[1]]
  aa_chr <- as.character(sequence)
  aa_len_including <- length(sequence)
  terminal_stop_included <- endsWith(aa_chr, "*")
  aa_len_excluding <- aa_len_including - as.integer(terminal_stop_included)
  body <- if (terminal_stop_included) substr(aa_chr, 1, nchar(aa_chr) - 1) else aa_chr
  internal_stops <- sum(strsplit(body, "", fixed = TRUE)[[1]] == "*")
  invalid_residues <- setdiff(unique(strsplit(body, "", fixed = TRUE)[[1]]), Biostrings::AA_STANDARD)
  invalid_count <- sum(strsplit(body, "", fixed = TRUE)[[1]] %in% invalid_residues)
  aa_length_match <- !is.na(expected_aa_len) && expected_aa_len > 0 && aa_len_excluding == expected_aa_len
  failures <- character()
  if (!aa_length_match) failures <- c(failures, sprintf("protein length %d aa (excluding terminal stop) does not match expected %s aa", aa_len_excluding, expected_aa_len))
  if (internal_stops > 0) failures <- c(failures, sprintf("internal stop codons: %d", internal_stops))
  if (invalid_count > 0) failures <- c(failures, sprintf("invalid/ambiguous amino-acid residues: %d", invalid_count))

  list(
    reference_type = "protein", seq_length_bp = NA_integer_, expected_bp = NA_integer_,
    cds_length_match = NA, start_codon = NA_character_, start_codon_valid = NA,
    stop_codon = NA_character_, stop_codon_valid = NA, cds_divisible_by_3 = NA,
    internal_stops = internal_stops, ambiguous_bases = invalid_count,
    aa_length_including_terminal_stop = aa_len_including,
    aa_length_excluding_terminal_stop = aa_len_excluding,
    expected_aa_len = expected_aa_len, aa_length_match = aa_length_match,
    terminal_stop_included = terminal_stop_included, gc_percent = NA_real_,
    translation_success = TRUE,
    status = if (length(failures) == 0) "PASS" else "FAIL",
    reason = if (length(failures) == 0) "OK" else paste(failures, collapse = "; "),
    translation = sequence
  )
}

validation_log <- list()
translated_proteins <- list()
overall_status <- "PASS"
required_cols <- c("gene", "reference_file", "reference_type", "expected_aa_len", "CDS_expected_bp_length", "retrieval_strategy", "notes")
missing_cols <- setdiff(required_cols, colnames(TARGETS))
if (length(missing_cols) > 0 || anyDuplicated(TARGETS$gene) || anyNA(TARGETS$gene) || any(!nzchar(TARGETS$gene))) {
  overall_status <- "FAIL"
  schema_reason <- c(
    if (length(missing_cols) > 0) sprintf("missing columns: %s", paste(missing_cols, collapse = ", ")),
    if ("gene" %in% colnames(TARGETS) && anyDuplicated(TARGETS$gene)) "duplicate gene identifiers",
    if ("gene" %in% colnames(TARGETS) && (anyNA(TARGETS$gene) || any(!nzchar(TARGETS$gene)))) "missing gene identifiers"
  )
  log_error(sprintf("TARGETS schema validation failed: %s", paste(schema_reason, collapse = "; ")), LOG_FILE)
  stop("Execution halted: target metadata schema mismatch")
}

n_genes <- nrow(TARGETS)
for (i in seq_len(n_genes)) {
  gene <- TARGETS$gene[i]
  ref_file <- TARGETS$reference_file[i]
  metadata_type <- normalize_reference_type(TARGETS$reference_type[i])
  log_info(sprintf("\n--- Validating: %s ---", gene), LOG_FILE)

  if (!file.exists(ref_file)) {
    log_error(sprintf("  [FAIL] File not found: %s", ref_file), LOG_FILE)
    validation_log[[gene]] <- empty_result(gene, ref_file, metadata_type, reason = "FASTA file not found")
    overall_status <- "FAIL"
    next
  }

  ref_type <- metadata_type
  if (is.na(ref_type)) {
    ref_type <- tryCatch(detect_seq_type_from_content(ref_file), error = function(e) NA_character_)
    if (is.na(ref_type)) {
      validation_log[[gene]] <- empty_result(gene, ref_file, reason = "sequence-type detection failed")
      log_error(sprintf("  [FAIL] Unable to determine reference type for %s", gene), LOG_FILE)
      overall_status <- "FAIL"
      next
    }
    log_warn(sprintf("  [WARN] %s has an unrecognised metadata reference_type; content fallback classified it as %s", gene, ref_type), LOG_FILE)
  }

  seq_set <- tryCatch(
    if (ref_type == "nucleotide") Biostrings::readDNAStringSet(ref_file) else Biostrings::readAAStringSet(ref_file),
    error = function(e) {
      log_error(sprintf("  [FAIL] Cannot parse %s FASTA for %s: %s", ref_type, gene, conditionMessage(e)), LOG_FILE)
      NULL
    }
  )
  if (is.null(seq_set) || length(seq_set) == 0) {
    validation_log[[gene]] <- empty_result(gene, ref_file, ref_type, reason = "FASTA parsing failed or contains no records")
    overall_status <- "FAIL"
    next
  }
  if (length(seq_set) != 1) {
    log_error(sprintf("  [FAIL] %s contains %d FASTA records; exactly one is required.", gene, length(seq_set)), LOG_FILE)
    validation_log[[gene]] <- empty_result(gene, ref_file, ref_type, reason = sprintf("expected exactly one FASTA record; found %d", length(seq_set)))
    overall_status <- "FAIL"
    next
  }

  sequence <- seq_set[[1]]
  header <- names(seq_set)[1]
  checksum <- sequence_checksum(sequence)
  # Source FASTA records can use a stable locus identifier instead of the local
  # target alias. A valid header is therefore non-empty and identifies either
  # the target alias or a source `gene=` field; alias presence is recorded
  # separately rather than treated as a false scientific failure.
  header_contains_target_gene <- !is.na(header) && grepl(gene, header, ignore.case = TRUE, fixed = TRUE)
  header_valid <- !is.na(header) && nzchar(trimws(header)) &&
    (header_contains_target_gene || grepl("gene=", header, fixed = TRUE))
  if (!header_valid) {
    log_error(sprintf("  [FAIL] Header is empty or lacks a target alias/source gene identifier: %s", header), LOG_FILE)
  } else if (!header_contains_target_gene) {
    log_info(sprintf("  [INFO] Header uses a source gene identifier rather than project alias '%s'.", gene), LOG_FILE)
  }

  result <- if (ref_type == "nucleotide") {
    validate_nucleotide_ref(gene, seq_set, TARGETS$CDS_expected_bp_length[i], TARGETS$expected_aa_len[i], LOG_FILE)
  } else {
    validate_protein_ref(seq_set, TARGETS$expected_aa_len[i])
  }
  result$gene <- gene
  result$file <- ref_file
  result$header <- header
  result$checksum <- checksum
  result$header_valid <- header_valid
  result$header_contains_target_gene <- header_contains_target_gene
  result$header_unique <- NA
  if (!header_valid) {
    result$status <- "FAIL"
    result$reason <- append_reason(result$reason, "header is empty or lacks a target/source gene identifier")
  }
  validation_log[[gene]] <- result
  if (result$status == "PASS") translated_proteins[[gene]] <- result$translation
  if (result$status == "FAIL") overall_status <- "FAIL"
}

# Verify header uniqueness only after every readable, single-record FASTA has
# been collected. It is a cross-reference validation, not an assumption.
headers <- vapply(validation_log, function(x) if (is.null(x$header)) NA_character_ else as.character(x$header), character(1))
duplicate_headers <- unique(headers[!is.na(headers) & duplicated(headers)])
for (gene in names(validation_log)) {
  result <- validation_log[[gene]]
  result$header_unique <- !is.na(result$header) && !(result$header %in% duplicate_headers)
  if (!result$header_unique) {
    result$status <- "FAIL"
    result$reason <- append_reason(result$reason, "FASTA header is not unique across targets")
    translated_proteins[[gene]] <- NULL
    overall_status <- "FAIL"
  }
  validation_log[[gene]] <- result
}

summary_rows <- lapply(validation_log, function(x) {
  data.frame(
    gene = x$gene, reference_file = basename(x$file), header = x$header,
    checksum = x$checksum, reference_type = x$reference_type,
    observed_bp = x$seq_length_bp, expected_bp = x$expected_bp,
    cds_length_match = x$cds_length_match, cds_divisible_by_3 = x$cds_divisible_by_3,
    start_codon = x$start_codon, start_codon_valid = x$start_codon_valid,
    stop_codon = x$stop_codon, stop_codon_valid = x$stop_codon_valid,
    internal_stops = x$internal_stops, ambiguous_bases = x$ambiguous_bases,
    aa_length_including_terminal_stop = x$aa_length_including_terminal_stop,
    aa_length_excluding_terminal_stop = x$aa_length_excluding_terminal_stop,
    expected_aa_len = x$expected_aa_len, aa_length_match = x$aa_length_match,
    terminal_stop_included = x$terminal_stop_included, gc_percent = x$gc_percent,
    header_valid = x$header_valid, header_contains_target_gene = x$header_contains_target_gene,
    header_unique = x$header_unique,
    translation_success = x$translation_success, status = x$status, reason = x$reason,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[match(TARGETS$gene, summary_df$gene), , drop = FALSE]
write.csv(summary_df, file.path(PATHS$reports, "reference_summary.csv"), row.names = FALSE)
log_info("Saved reference_summary.csv to reports directory.", LOG_FILE)

failed_genes <- summary_df$gene[summary_df$status == "FAIL"]
all_headers_unique <- all(summary_df$header_unique %in% TRUE)
report_lines <- c(
  "==================================", "REFERENCE VALIDATION REPORT", "==================================", "",
  sprintf("Date: %s", Sys.Date()), sprintf("Pipeline: %s", PROJECT$name),
  sprintf("Organism: %s %s", PROJECT$species, PROJECT$strain),
  sprintf("Reference file directory: %s", PATHS$ref), sprintf("Total genes: %d", n_genes),
  sprintf("Overall status: %s", overall_status),
  sprintf("Header uniqueness explicitly checked: %s", if (all_headers_unique) "PASS" else "FAIL"), ""
)
if (overall_status == "PASS") {
  report_lines <- c(report_lines, "All target references passed the checks applicable to their declared sequence type.",
                    "No Rsubread index is constructed in this validation step.")
} else {
  report_lines <- c(report_lines, sprintf("Failed genes: %s", paste(failed_genes, collapse = ", ")),
                    "Correct the reported reference integrity failures before mapping.")
}
report_lines <- c(report_lines, "", "==================================", "DETAILED STATISTICS", "==================================")
for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, ]
  report_lines <- c(report_lines,
    sprintf("%s: %s", row$gene, row$status),
    sprintf("  File: %s", row$reference_file),
    sprintf("  Header: %s", row$header),
    sprintf("  Sequence MD5: %s", row$checksum),
    sprintf("  Header valid/contains target alias/unique: %s / %s / %s", row$header_valid, row$header_contains_target_gene, row$header_unique),
    sprintf("  CDS length: %s bp (expected %s; match %s; divisible by 3 %s)", row$observed_bp, row$expected_bp, row$cds_length_match, row$cds_divisible_by_3),
    sprintf("  Start/stop: %s (%s) / %s (%s)", row$start_codon, row$start_codon_valid, row$stop_codon, row$stop_codon_valid),
    sprintf("  Protein length: %s aa including terminal stop; %s aa excluding terminal stop (expected %s; match %s)", row$aa_length_including_terminal_stop, row$aa_length_excluding_terminal_stop, row$expected_aa_len, row$aa_length_match),
    sprintf("  Internal stops: %s; ambiguous residues/bases: %s; translation success: %s", row$internal_stops, row$ambiguous_bases, row$translation_success),
    sprintf("  Validation reason: %s", row$reason), ""
  )
}
writeLines(report_lines, file.path(PATHS$reports, "reference_validation_report.txt"))
log_info("Saved reference_validation_report.txt to reports directory.", LOG_FILE)

if (overall_status == "PASS") {
  validated_aas <- Biostrings::AAStringSet(vapply(translated_proteins, as.character, character(1)))
  names(validated_aas) <- names(translated_proteins)
  saveRDS(validated_aas, file.path(PATHS$intermediate, "reference_sequences.rds"))
  log_info(sprintf("Validated protein translations saved to: %s", file.path(PATHS$intermediate, "reference_sequences.rds")), LOG_FILE)
} else {
  log_warn("Skipping RDS export because reference integrity validation failed.", LOG_FILE)
}

log_system_memory(LOG_FILE)
log_info(sprintf("Step 2 Reference Integrity Verification complete: %s.", overall_status), LOG_FILE)
file.copy(LOG_FILE, file.path(PATHS$logs, "latest_execution.log"), overwrite = TRUE)

if (overall_status == "PASS") {
  cat(sprintf("\n[SUCCESS] Reference validation PASSED for all %d genes.\n", n_genes))
} else {
  cat(sprintf("\n[ERROR] Reference validation FAILED for: %s\n", paste(failed_genes, collapse = ", ")))
  stop("Step 2 reference integrity validation failed; see reports/reference_validation_report.txt")
}
