# ==============================================================================
# Step 4B.3c-1: BLAST-HSP Ground-Truth Orientation/Coordinate Validation
#
# Decisive comparison: BLAST sseq <-> independent translation of FASTQ
# nucleotides extracted using BLAST's sstart/send for that HSP.
# ==============================================================================

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("CONFIG")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

GT_DIR <- file.path(PATHS$reports, "tac102_4B3c1_ground_truth")
if (!dir.exists(GT_DIR)) dir.create(GT_DIR, recursive = TRUE, showWarnings = FALSE)

log_info("Starting Step 4B.3c-1: BLAST-HSP ground-truth validation.", LOG_FILE)

# ------------------------------------------------------------------------------
# 1. Selection: mandatory 3 + deterministic quantile spread, stratified by strand
# ------------------------------------------------------------------------------
tails <- read.csv(file.path(PATHS$reports, "tac102_41pairs_unaligned_tails.csv"), stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# Explicitly select the originally flagged HSP/mate for mandatory validation.
# This avoids ambiguity when both mates have non-NA gap_position_start.
# ------------------------------------------------------------------------------

mandatory_mate_map <- c(
  "SRR7910035.1964432" = "R1",
  "SRR7910035.24197645" = "R1",
  "SRR7910035.3831778" = "R2"
)

mandatory_rows <- do.call(rbind, lapply(names(mandatory_mate_map), function(pid) {

  target_mate <- mandatory_mate_map[[pid]]

  d <- tails[
    tails$pair_id == pid &
    tails$mate == target_mate &
    !is.na(tails$hsp_strand),
  ]

  if (nrow(d) == 0) {
    stop(sprintf(
      "Explicitly required mandatory row not found: %s / %s",
      pid,
      target_mate
    ))
  }

  d[1, ]
}))

log_info(
  sprintf(
    "Mandatory HSPs explicitly selected: %s",
    paste(
      sprintf(
        "%s(%s,%s)",
        mandatory_rows$pair_id,
        mandatory_rows$mate,
        mandatory_rows$hsp_strand
      ),
      collapse = "; "
    )
  ),
  LOG_FILE
)

cat("\n=== MANDATORY HSP SELECTION ===\n")
print(
  mandatory_rows[, c(
    "pair_id",
    "mate",
    "hsp_strand",
    "pident_aa",
    "gap_position_start"
  )],
  row.names = FALSE
)

target_forward <- 4L; target_reverse <- 4L
need_fwd <- target_forward - sum(mandatory_rows$hsp_strand == "forward")
need_rev <- target_reverse - sum(mandatory_rows$hsp_strand == "reverse")

pick_by_quantile <- function(pool, n) {
  pool <- pool[!is.na(pool$pident_aa), ]
  pool <- pool[order(pool$pident_aa), ]
  if (nrow(pool) == 0 || n <= 0) return(pool[0, ])
  if (nrow(pool) <= n) return(pool)
  idx <- unique(round(seq(1, nrow(pool), length.out = n)))
  pool[idx, ]
}

pool_fwd <- tails[tails$hsp_strand == "forward" & !tails$pair_id %in% names(mandatory_mate_map), ]
pool_rev <- tails[tails$hsp_strand == "reverse" & !tails$pair_id %in% names(mandatory_mate_map), ]
extra_fwd <- pick_by_quantile(pool_fwd, need_fwd)
extra_rev <- pick_by_quantile(pool_rev, need_rev)

selected <- rbind(mandatory_rows, extra_fwd, extra_rev)
selected <- selected[!duplicated(paste(selected$pair_id, selected$mate)), ]

write.csv(selected[, c("pair_id","mate","hsp_strand","pident_aa")],
          file.path(GT_DIR, "selected_hsps_audit_trail.csv"), row.names = FALSE)

cat("\n=== DETERMINISTIC SELECTION AUDIT TRAIL ===\n")
print(selected[, c("pair_id","mate","hsp_strand","pident_aa")])
cat("===========================================\n\n")

log_info(sprintf("Selected %d HSPs (mandatory: %d, forward-fill: %d, reverse-fill: %d).",
                 nrow(selected), nrow(mandatory_rows), nrow(extra_fwd), nrow(extra_rev)), LOG_FILE)

# ------------------------------------------------------------------------------
# 2. Targeted tblastn re-query
# ------------------------------------------------------------------------------
GT_OUTFMT <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
tac_query <- file.path(PATHS$raw_recovery_derivative, "queries", "TAC102.fasta")

run_gt_tblastn <- function(subject_fasta, out_tsv) {
  args <- c("-query", tac_query, "-subject", subject_fasta,
            "-outfmt", shQuote(GT_OUTFMT), "-evalue", "1e-10", "-out", out_tsv)
  res <- system2("tblastn", args = args, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(res, "status")
  if (!is.null(exit_status) && exit_status != 0) stop(paste("tblastn failed:", paste(res, collapse="\n")))
  out_tsv
}

r1_gt <- run_gt_tblastn(file.path(PATHS$raw_recovery_derivative, "tac102_41pair_R1.fasta"), file.path(GT_DIR, "tac102_ground_truth_R1.tsv"))
r2_gt <- run_gt_tblastn(file.path(PATHS$raw_recovery_derivative, "tac102_41pair_R2.fasta"), file.path(GT_DIR, "tac102_ground_truth_R2.tsv"))

gt_cols <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qseq","sseq")
gt_r1 <- read.delim(r1_gt, header = FALSE, col.names = gt_cols, stringsAsFactors = FALSE); gt_r1$mate <- "R1"
gt_r2 <- read.delim(r2_gt, header = FALSE, col.names = gt_cols, stringsAsFactors = FALSE); gt_r2$mate <- "R2"
gt_all <- rbind(gt_r1, gt_r2)
log_info(sprintf("Ground-truth tblastn: %d HSPs (R1:%d, R2:%d).", nrow(gt_all), nrow(gt_r1), nrow(gt_r2)), LOG_FILE)

# ------------------------------------------------------------------------------
# 3. Four explicit tests per selected HSP (with modulus-3 safeguard)
# ------------------------------------------------------------------------------
r1_seqs <- Biostrings::readDNAStringSet(file.path(PATHS$raw_recovery_derivative, "tac102_41pair_R1.fasta"))
r2_seqs <- Biostrings::readDNAStringSet(file.path(PATHS$raw_recovery_derivative, "tac102_41pair_R2.fasta"))
names(r1_seqs) <- sub("\\s.*$", "", names(r1_seqs))
names(r2_seqs) <- sub("\\s.*$", "", names(r2_seqs))

pct_identity_or_na <- function(a, b) {
  if (is.na(a) || is.na(b) || nchar(a) == 0 || nchar(b) == 0) return(NA_real_)
  if (nchar(a) != nchar(b)) return(NA_real_)
  mean(strsplit(a, "")[[1]] == strsplit(b, "")[[1]]) * 100
}

make_empty_row <- function(pair_id, mate, status, note, forward = NA, sstart = NA_integer_, send = NA_integer_,
                           qstart = NA_integer_, qend = NA_integer_, pident = NA_real_, hsp_len = NA_integer_,
                           bitscore = NA_real_, raw_nt = "", raw_nt_length = 0L, blast_sseq = "", blast_qseq = "") {
  b_sseq_un <- gsub("-", "", blast_sseq)
  b_qseq_un <- gsub("-", "", blast_qseq)
  data.frame(
    pair_id = pair_id, mate = mate,
    hsp_strand = if (is.na(forward)) NA_character_ else if (forward) "forward" else "reverse",
    sstart = sstart, send = send, qstart = qstart, qend = qend,
    reported_pident = pident, hsp_length = hsp_len, bitscore = bitscore,
    raw_nt = raw_nt, raw_nt_length = raw_nt_length,
    our_translation = NA_character_, our_translation_length = NA_integer_,
    blast_sseq = blast_sseq, blast_sseq_ungapped = b_sseq_un, blast_sseq_length = nchar(b_sseq_un),
    blast_qseq = blast_qseq, blast_qseq_ungapped = b_qseq_un, blast_qseq_length = nchar(b_qseq_un),
    reconstruction_matches_blast_sseq_pct = NA_real_,
    reconstruction_length_difference = NA_integer_,
    blast_sseq_vs_qseq_identity_pct = NA_real_,
    coordinate_orientation_status = status,
    note = note, stringsAsFactors = FALSE
  )
}

validate_hsp <- function(pair_id, mate) {
  hsp <- gt_all[gt_all$sseqid == pair_id & gt_all$mate == mate, ]
  if (nrow(hsp) == 0) {
    return(make_empty_row(pair_id, mate, "NO_HSP", "no ground-truth HSP found"))
  }
  hsp <- hsp[which.max(hsp$bitscore), ]

  read_seq <- if (mate == "R1") r1_seqs[[pair_id]] else r2_seqs[[pair_id]]
  if (is.null(read_seq)) {
    return(make_empty_row(pair_id, mate, "ERROR", "read sequence not found in fasta",
                          forward = (hsp$sstart < hsp$send), sstart = hsp$sstart, send = hsp$send,
                          qstart = hsp$qstart, qend = hsp$qend, pident = hsp$pident,
                          hsp_len = hsp$length, bitscore = hsp$bitscore,
                          blast_sseq = hsp$sseq, blast_qseq = hsp$qseq))
  }

  forward <- hsp$sstart < hsp$send
  s <- min(hsp$sstart, hsp$send)
  e <- max(hsp$sstart, hsp$send)
  raw_nt_obj <- Biostrings::subseq(read_seq, start = s, end = e)
  raw_nt <- as.character(raw_nt_obj)
  raw_nt_length <- nchar(raw_nt)

  # Safeguard: Verify raw nucleotide interval length is a multiple of 3
  if (raw_nt_length %% 3L != 0L) {
    return(make_empty_row(
      pair_id, mate, "NON_MULTIPLE_OF_3",
      "BLAST sstart:send interval is not divisible by 3; translation skipped",
      forward = forward, sstart = hsp$sstart, send = hsp$send,
      qstart = hsp$qstart, qend = hsp$qend, pident = hsp$pident,
      hsp_len = hsp$length, bitscore = hsp$bitscore,
      raw_nt = raw_nt, raw_nt_length = raw_nt_length,
      blast_sseq = hsp$sseq, blast_qseq = hsp$qseq
    ))
  }

  our_translation <- tryCatch({
    seq_to_translate <- if (forward) raw_nt_obj else Biostrings::reverseComplement(raw_nt_obj)
    as.character(Biostrings::translate(seq_to_translate, if.fuzzy.codon = "X", no.init.codon = TRUE))
  }, error = function(e) NA_character_)

  blast_sseq_ungapped <- gsub("-", "", hsp$sseq)
  blast_qseq_ungapped <- gsub("-", "", hsp$qseq)

  len_match_3 <- nchar(our_translation) == nchar(blast_sseq_ungapped)
  pident_3 <- pct_identity_or_na(our_translation, blast_sseq_ungapped)
  pident_4 <- pct_identity_or_na(blast_sseq_ungapped, blast_qseq_ungapped)

  status <- if (is.na(our_translation)) {
    "ERROR"
  } else if (!len_match_3) {
    "LENGTH_MISMATCH"
  } else if (!is.na(pident_3) && pident_3 == 100) {
    "PASS"
  } else {
    "FAIL"
  }

  data.frame(
    pair_id = pair_id, mate = mate,
    hsp_strand = if (forward) "forward" else "reverse",
    sstart = hsp$sstart, send = hsp$send, qstart = hsp$qstart, qend = hsp$qend,
    reported_pident = hsp$pident, hsp_length = hsp$length, bitscore = hsp$bitscore,
    raw_nt = raw_nt, raw_nt_length = raw_nt_length,
    our_translation = our_translation, our_translation_length = nchar(our_translation),
    blast_sseq = hsp$sseq, blast_sseq_ungapped = blast_sseq_ungapped, blast_sseq_length = nchar(blast_sseq_ungapped),
    blast_qseq = hsp$qseq, blast_qseq_ungapped = blast_qseq_ungapped, blast_qseq_length = nchar(blast_qseq_ungapped),
    reconstruction_matches_blast_sseq_pct = pident_3,
    reconstruction_length_difference = nchar(our_translation) - nchar(blast_sseq_ungapped),
    blast_sseq_vs_qseq_identity_pct = pident_4,
    coordinate_orientation_status = status,
    note = "", stringsAsFactors = FALSE
  )
}

gt_results <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
  tryCatch(
    validate_hsp(selected$pair_id[i], selected$mate[i]),
    error = function(e) make_empty_row(selected$pair_id[i], selected$mate[i], "ERROR",
                                       paste("uncaught error:", conditionMessage(e)))
  )
}))

write.csv(gt_results, file.path(GT_DIR, "tac102_ground_truth_validation_results.csv"), row.names = FALSE)

cat("\n=== GROUND-TRUTH VALIDATION SUMMARY ===\n")
print(table(strand = gt_results$hsp_strand, status = gt_results$coordinate_orientation_status, useNA = "ifany"))
cat("=======================================\n")

log_info(sprintf("Step 4B.3c-1 completed. Status distribution: %s",
                 paste(capture.output(table(gt_results$hsp_strand, gt_results$coordinate_orientation_status)), collapse=" | ")), LOG_FILE)

file.copy(LOG_FILE, file.path(PATHS$logs, "latest_execution.log"), overwrite = TRUE)
cat("\nFull results saved to:", file.path(GT_DIR, "tac102_ground_truth_validation_results.csv"), "\n")