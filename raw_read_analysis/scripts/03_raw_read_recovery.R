# ==============================================================================
# Step 3: Raw Read Recovery — Translated-Read Protein Homology Search
# Project: Phase II Validation of Candidate Genes from T. equiperdum IVM-t1 Raw Reads
# File: scripts/03_raw_read_recovery.R
# ==============================================================================

cat("======================================================\n")
cat(" Phase II: Candidate Gene Recovery from Raw Reads      \n")
cat(" Step 3: Raw Read Recovery (03_raw_read_recovery.R)    \n")
cat("======================================================\n\n")

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("TARGETS") || !exists("CONFIG")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

log_info("Starting Step 3: Raw Read Recovery...", LOG_FILE)

RR <- CONFIG$raw_recovery
if (is.null(RR)) stop("Execution halted: CONFIG$raw_recovery not defined. Update 00_setup.R first.")

ref_seqs_path <- file.path(PATHS$intermediate, "reference_sequences.rds")
if (!file.exists(ref_seqs_path)) {
  log_error(sprintf("Missing Step 2 output: %s", ref_seqs_path), LOG_FILE)
  stop("Execution halted: run 02_reference_prep.R first.")
}
reference_sequences <- readRDS(ref_seqs_path)
log_info(sprintf("Loaded %d validated protein reference(s) from Step 2.", length(reference_sequences)), LOG_FILE)

# ------------------------------------------------------------------------------
# 3.1 Derivative FASTA Construction
# ------------------------------------------------------------------------------
fastq_to_fasta_derivative <- function(fastq_path, fasta_out, log_file = NULL) {
  if (file.exists(fasta_out)) {
    log_info(sprintf("Derivative FASTA already exists, reusing: %s", basename(fasta_out)), log_file)
    return(invisible(fasta_out))
  }
  log_info(sprintf("Streaming %s -> derivative FASTA %s ...", basename(fastq_path), basename(fasta_out)), log_file)
  streamer <- ShortRead::FastqStreamer(fastq_path, n = 2e6)
  on.exit(close(streamer))
  first_chunk <- TRUE
  n_written <- 0
  while (length(fq <- ShortRead::yield(streamer))) {
    seqs <- ShortRead::sread(fq)
    ids  <- as.character(ShortRead::id(fq))
    names(seqs) <- ids
    Biostrings::writeXStringSet(seqs, filepath = fasta_out, append = !first_chunk, format = "fasta")
    first_chunk <- FALSE
    n_written <- n_written + length(fq)
  }
  log_info(sprintf("Wrote %s reads to %s.", format(n_written, big.mark = ","), basename(fasta_out)), log_file)
  invisible(fasta_out)
}

fasta_r1 <- file.path(PATHS$raw_recovery_derivative, "SRR7910035_1.fasta")
fasta_r2 <- file.path(PATHS$raw_recovery_derivative, "SRR7910035_2.fasta")
fastq_to_fasta_derivative(PATHS$fastq_r1, fasta_r1, LOG_FILE)
fastq_to_fasta_derivative(PATHS$fastq_r2, fasta_r2, LOG_FILE)

# ------------------------------------------------------------------------------
# 3.1b Read-ID Convention Detection
# ------------------------------------------------------------------------------
detect_pair_id_pattern <- function(fasta_r1_path, fasta_r2_path, log_file, n_sample = 500) {
  r1_ids <- names(Biostrings::readDNAStringSet(fasta_r1_path, nrec = n_sample))
  r2_ids <- names(Biostrings::readDNAStringSet(fasta_r2_path, nrec = n_sample))

  log_info(sprintf("Sample R1 read ID: %s", r1_ids[1]), log_file)
  log_info(sprintf("Sample R2 read ID: %s", r2_ids[1]), log_file)

  candidates <- list(
    first_token   = function(x) sub("\\s.*$", "", x),
    dot_suffix    = function(x) sub("\\.[12]$", "", x),
    slash_suffix  = function(x) sub("/[12]$", "", x),
    space_suffix  = function(x) sub("\\s+[12](:.*)?$", "", x),
    identical_ids = function(x) x
  )

  for (pattern_name in names(candidates)) {
    strip_fn <- candidates[[pattern_name]]
    r1_stripped <- strip_fn(r1_ids)
    r2_stripped <- strip_fn(r2_ids)
    match_rate <- mean(r1_stripped == r2_stripped)
    if (match_rate == 1) {
      log_info(sprintf("Detected read-ID pairing convention: '%s' (100%% match on %d-record sample).",
                       pattern_name, n_sample), log_file)
      return(list(pattern = pattern_name, strip_fn = strip_fn, confirmed = TRUE))
    }
  }

  log_warn(sprintf(
    "Could not confirm a consistent R1/R2 ID pairing convention from a %d-record sample. Mate concordance will be reported as NOT_DETERMINED rather than assumed.",
    n_sample), log_file)
  list(pattern = "undetermined", strip_fn = NULL, confirmed = FALSE)
}

pair_id_info <- detect_pair_id_pattern(fasta_r1, fasta_r2, LOG_FILE)

# ------------------------------------------------------------------------------
# 3.2 BLAST Database Construction
# ------------------------------------------------------------------------------
build_blastdb <- function(fasta_path, db_prefix, log_file = NULL) {
  db_check <- paste0(db_prefix, ".nsq")
  db_check_alt <- paste0(db_prefix, ".00.nsq")
  if (file.exists(db_check) || file.exists(db_check_alt)) {
    log_info(sprintf("BLAST db already exists, reusing: %s", basename(db_prefix)), log_file)
    return(invisible(db_prefix))
  }
  log_info(sprintf("Building BLAST nucleotide db: %s", basename(db_prefix)), log_file)
  args <- c("-in", fasta_path, "-dbtype", "nucl", "-out", db_prefix, "-parse_seqids")
  res <- system2("makeblastdb", args = args, stdout = TRUE, stderr = TRUE)
  writeLines(res, file.path(PATHS$logs, sprintf("makeblastdb_%s.log", basename(db_prefix))))
  if (!file.exists(db_check) && !file.exists(db_check_alt)) {
    log_error(sprintf("makeblastdb did not produce expected output for %s", db_prefix), log_file)
    stop("Execution halted: BLAST database creation failed.")
  }
  invisible(db_prefix)
}

db_r1 <- file.path(PATHS$blastdb, "ivm_t1_R1_blastdb")
db_r2 <- file.path(PATHS$blastdb, "ivm_t1_R2_blastdb")
build_blastdb(fasta_r1, db_r1, LOG_FILE)
build_blastdb(fasta_r2, db_r2, LOG_FILE)

# ------------------------------------------------------------------------------
# 3.3 tBLASTn Execution (Layer A)
# ------------------------------------------------------------------------------
colnames_tblastn <- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                      "qstart","qend","sstart","send","evalue","bitscore")

run_tblastn <- function(gene, query_fasta_path, db_prefix, mate_label, out_path, log_file = NULL) {
  if (file.exists(out_path) && file.info(out_path)$size > 0) {
    log_info(sprintf("%s (%s): Layer A output already exists, reusing.", gene, mate_label), log_file)
    return(invisible(out_path))
  }
  log_info(sprintf("%s (%s): running tblastn against %s ...", gene, mate_label, basename(db_prefix)), log_file)
  
  args <- c(
    "-query", query_fasta_path,
    "-db", db_prefix,
    "-out", out_path,
    "-outfmt", shQuote(RR$blast_outfmt),
    "-evalue", format(RR$blast_evalue, scientific = TRUE),
    "-max_target_seqs", format(RR$blast_max_target_seqs, scientific = FALSE),
    "-num_threads", RR$blast_threads
  )
  
  res <- system2("tblastn", args = args, stdout = TRUE, stderr = TRUE)
  exit_status <- attr(res, "status")
  writeLines(res, file.path(PATHS$logs, sprintf("tblastn_%s_%s.log", gene, mate_label)))
  
  if (!is.null(exit_status) && exit_status != 0) {
    log_error(sprintf("%s (%s): tblastn exited with status %d. See log: tblastn_%s_%s.log",
                       gene, mate_label, exit_status, gene, mate_label), log_file)
    stop(sprintf("Execution halted: tblastn failed for %s (%s). Check logs/tblastn_%s_%s.log",
                 gene, mate_label, gene, mate_label))
  }
  
  if (!file.exists(out_path)) {
    log_error(sprintf("%s (%s): tblastn reported success but no output file was created.", gene, mate_label), log_file)
    stop(sprintf("Execution halted: expected tblastn output missing for %s (%s).", gene, mate_label))
  }
  
  n_hits <- length(readLines(out_path))
  log_info(sprintf("%s (%s): tblastn complete, %d raw HSP row(s) (Layer A).", gene, mate_label, n_hits), log_file)
  invisible(out_path)
}

gene_names <- names(reference_sequences)
query_dir <- file.path(PATHS$raw_recovery_derivative, "queries")
if (!dir.exists(query_dir)) dir.create(query_dir, recursive = TRUE, showWarnings = FALSE)

expected_aa_len <- Biostrings::width(reference_sequences)
names(expected_aa_len) <- names(reference_sequences)

log_info(sprintf("expected_aa_len built for %d gene(s): %s",
                  length(expected_aa_len),
                  paste(sprintf("%s=%d", names(expected_aa_len), expected_aa_len), collapse = ", ")),
          LOG_FILE)

for (gene in gene_names) {
  query_path <- file.path(query_dir, sprintf("%s.fasta", gene))
  if (!file.exists(query_path)) {
    aa_seq <- reference_sequences[gene]
    Biostrings::writeXStringSet(aa_seq, filepath = query_path, format = "fasta")
  }

  out_r1 <- file.path(PATHS$tblastn_raw, sprintf("%s_R1.tsv", gene))
  out_r2 <- file.path(PATHS$tblastn_raw, sprintf("%s_R2.tsv", gene))
  run_tblastn(gene, query_path, db_r1, "R1", out_r1, LOG_FILE)
  run_tblastn(gene, query_path, db_r2, "R2", out_r2, LOG_FILE)
}

# ------------------------------------------------------------------------------
# 3.4 Layer A' — Consolidated & Normalized HSP Derivative
# ------------------------------------------------------------------------------
merge_intervals <- function(starts, ends) {
  ord <- order(starts)
  starts <- starts[ord]; ends <- ends[ord]
  merged_start <- starts[1]; merged_end <- ends[1]
  out_start <- integer(0); out_end <- integer(0)
  if (length(starts) > 1) {
    for (i in 2:length(starts)) {
      if (starts[i] <= merged_end + 1) {
        merged_end <- max(merged_end, ends[i])
      } else {
        out_start <- c(out_start, merged_start); out_end <- c(out_end, merged_end)
        merged_start <- starts[i]; merged_end <- ends[i]
      }
    }
  }
  out_start <- c(out_start, merged_start); out_end <- c(out_end, merged_end)
  list(start = out_start, end = out_end, covered_len = sum(out_end - out_start + 1))
}

read_hits <- function(path, gene, mate_label) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(matrix(ncol = length(colnames_tblastn) + 2, nrow = 0,
                             dimnames = list(NULL, c(colnames_tblastn, "gene", "mate")))))
  }
  df <- read.delim(path, header = FALSE, col.names = colnames_tblastn, stringsAsFactors = FALSE)
  df$gene <- gene
  df$mate <- mate_label
  df
}

all_hits <- list()
for (gene in gene_names) {
  h1 <- read_hits(file.path(PATHS$tblastn_raw, sprintf("%s_R1.tsv", gene)), gene, "R1")
  h2 <- read_hits(file.path(PATHS$tblastn_raw, sprintf("%s_R2.tsv", gene)), gene, "R2")
  all_hits[[gene]] <- rbind(h1, h2)
}
all_hits_df <- do.call(rbind, all_hits)
rownames(all_hits_df) <- NULL

if (nrow(all_hits_df) > 0) {
  all_hits_df$qstart_raw <- all_hits_df$qstart
  all_hits_df$qend_raw   <- all_hits_df$qend
  all_hits_df$sstart_raw <- all_hits_df$sstart
  all_hits_df$send_raw   <- all_hits_df$send

  all_hits_df$qstart <- pmin(all_hits_df$qstart_raw, all_hits_df$qend_raw)
  all_hits_df$qend   <- pmax(all_hits_df$qstart_raw, all_hits_df$qend_raw)
  all_hits_df$sstart <- pmin(all_hits_df$sstart_raw, all_hits_df$send_raw)
  all_hits_df$send   <- pmax(all_hits_df$sstart_raw, all_hits_df$send_raw)
}

if (!is.null(RR$min_alignment_length_aa) && nrow(all_hits_df) > 0) {
  all_hits_df <- all_hits_df[all_hits_df$length >= RR$min_alignment_length_aa, ]
}
if (!is.null(RR$min_pident) && nrow(all_hits_df) > 0) {
  all_hits_df <- all_hits_df[all_hits_df$pident >= RR$min_pident, ]
}

write.csv(all_hits_df, file.path(PATHS$reports, "raw_read_recovery_all_hsps.csv"), row.names = FALSE)
log_info(sprintf("Layer A' consolidated HSP rows: %d", nrow(all_hits_df)), LOG_FILE)

# ------------------------------------------------------------------------------
# 3.5 Layer B — Per-Read Aggregation
# ------------------------------------------------------------------------------
derive_pair_id <- function(read_ids, strip_fn) {
  if (is.null(strip_fn)) return(rep(NA_character_, length(read_ids)))
  strip_fn(read_ids)
}

if (nrow(all_hits_df) > 0) {
  split_keys <- interaction(all_hits_df$gene, all_hits_df$mate, all_hits_df$sseqid, drop = TRUE)
  grouped <- split(all_hits_df, split_keys)
  per_read_rows <- lapply(grouped, function(g) {
    q_merged <- merge_intervals(g$qstart, g$qend)
    s_merged <- merge_intervals(g$sstart, g$send)
    best <- g[which.max(g$bitscore), ]
    data.frame(
      gene                = g$gene[1],
      mate                = g$mate[1],
      read_id             = g$sseqid[1],
      pair_id             = derive_pair_id(g$sseqid[1], pair_id_info$strip_fn),
      n_hsps              = nrow(g),
      qstart_union        = paste(q_merged$start, collapse = ";"),
      qend_union          = paste(q_merged$end, collapse = ";"),
      aa_covered_by_read  = q_merged$covered_len,
      sstart_union        = paste(s_merged$start, collapse = ";"),
      send_union          = paste(s_merged$end, collapse = ";"),
      best_pident         = best$pident,
      best_evalue         = best$evalue,
      best_bitscore       = best$bitscore,
      stringsAsFactors = FALSE
    )
  })
  per_read_df <- do.call(rbind, per_read_rows)
  rownames(per_read_df) <- NULL
} else {
  per_read_df <- data.frame(
    gene=character(), mate=character(), read_id=character(), pair_id=character(),
    n_hsps=integer(), qstart_union=character(), qend_union=character(),
    aa_covered_by_read=integer(), sstart_union=character(), send_union=character(),
    best_pident=numeric(), best_evalue=numeric(), best_bitscore=numeric()
  )
}

write.csv(per_read_df, file.path(PATHS$reports, "raw_read_recovery_per_read.csv"), row.names = FALSE)
log_info(sprintf("Saved Layer B (deduplicated reads): %d rows.", nrow(per_read_df)), LOG_FILE)

# ------------------------------------------------------------------------------
# 3.6 Layer C — Per-Gene Summary
# ------------------------------------------------------------------------------
summarize_gene <- function(gene, protein_len, per_read_df, depth_thresholds, pair_id_confirmed) {
  gdf <- per_read_df[per_read_df$gene == gene, ]

  r1_reads <- unique(gdf$read_id[gdf$mate == "R1"])
  r2_reads <- unique(gdf$read_id[gdf$mate == "R2"])
  n_r1 <- length(r1_reads); n_r2 <- length(r2_reads)
  n_total_reads <- n_r1 + n_r2

  support_depth_vec <- numeric(protein_len)
  if (nrow(gdf) > 0) {
    delta <- numeric(protein_len + 2)
    for (i in seq_len(nrow(gdf))) {
      starts <- as.integer(strsplit(gdf$qstart_union[i], ";")[[1]])
      ends   <- as.integer(strsplit(gdf$qend_union[i], ";")[[1]])
      for (k in seq_along(starts)) {
        s <- max(1, starts[k]); e <- min(protein_len, ends[k])
        if (s <= e) { delta[s] <- delta[s] + 1; delta[e + 1] <- delta[e + 1] - 1 }
      }
    }
    support_depth_vec <- cumsum(delta[1:protein_len])
  }

  coverage_pct <- if (protein_len > 0) mean(support_depth_vec >= 1) * 100 else NA_real_
  depth_frac_at <- setNames(
    sapply(depth_thresholds, function(t) mean(support_depth_vec >= t) * 100),
    sprintf("pct_protein_support_depth_ge_%d", depth_thresholds)
  )

  if (pair_id_confirmed && nrow(gdf) > 0) {
    pairs_r1 <- unique(gdf$pair_id[gdf$mate == "R1"])
    pairs_r2 <- unique(gdf$pair_id[gdf$mate == "R2"])
    pairs_union <- union(pairs_r1, pairs_r2)
    pairs_both  <- intersect(pairs_r1, pairs_r2)
    mate_concordance_pct <- if (length(pairs_union) > 0) length(pairs_both) / length(pairs_union) * 100 else NA_real_
    n_pairs_union <- length(pairs_union)
    n_pairs_both  <- length(pairs_both)
    concordance_status <- sprintf("%.2f", mate_concordance_pct)
  } else {
    mate_concordance_pct <- NA_real_
    n_pairs_union <- NA_integer_
    n_pairs_both  <- NA_integer_
    concordance_status <- "NOT_DETERMINED"
  }

  pident_vec <- gdf$best_pident
  evalue_vec <- gdf$best_evalue

  status <- if (n_total_reads == 0) {
    "NO_SIGNIFICANT_SUPPORT"
  } else if (coverage_pct >= 100) {
    "DETECTED_FULL_LENGTH_COVERAGE"
  } else {
    "DETECTED_PARTIAL_COVERAGE"
  }

  data.frame(
    gene = gene,
    protein_length_aa = protein_len,
    unique_supporting_reads_R1 = n_r1,
    unique_supporting_reads_R2 = n_r2,
    unique_supporting_reads_total = n_total_reads,
    protein_query_coverage_pct = round(coverage_pct, 2),
    mean_support_depth = round(mean(support_depth_vec), 3),
    median_support_depth = median(support_depth_vec),
    min_support_depth = min(support_depth_vec),
    max_support_depth = max(support_depth_vec),
    t(depth_frac_at),
    mean_pident = if (length(pident_vec) > 0) round(mean(pident_vec), 2) else NA,
    median_pident = if (length(pident_vec) > 0) round(median(pident_vec), 2) else NA,
    min_pident = if (length(pident_vec) > 0) round(min(pident_vec), 2) else NA,
    max_pident = if (length(pident_vec) > 0) round(max(pident_vec), 2) else NA,
    min_evalue = if (length(evalue_vec) > 0) min(evalue_vec) else NA,
    median_evalue = if (length(evalue_vec) > 0) median(evalue_vec) else NA,
    max_evalue = if (length(evalue_vec) > 0) max(evalue_vec) else NA,
    n_read_pairs_with_any_support = n_pairs_union,
    n_read_pairs_with_both_mates_supporting = n_pairs_both,
    mate_concordance_pct = mate_concordance_pct,
    mate_concordance_status = concordance_status,
    observation_status = status,
    stringsAsFactors = FALSE
  )
}

summary_rows <- lapply(gene_names, function(g) {
  if (!g %in% names(expected_aa_len)) {
    log_error(sprintf(
      "Gene '%s' not found in expected_aa_len names. Available names: %s",
      g, paste(names(expected_aa_len), collapse = ", ")
    ), LOG_FILE)
    stop(sprintf("Execution halted: no expected_aa_len entry for gene '%s'.", g))
  }
  protein_len <- unname(expected_aa_len[g])
  summarize_gene(g, protein_len, per_read_df, RR$depth_report_thresholds, pair_id_info$confirmed)
})
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(PATHS$reports, "raw_read_recovery_summary.csv"), row.names = FALSE)
log_info("Saved raw_read_recovery_summary.csv (Layer C).", LOG_FILE)

# ------------------------------------------------------------------------------
# 3.7 Human-readable report
# ------------------------------------------------------------------------------
report_lines <- c(
  "==================================",
  "RAW READ RECOVERY REPORT (Step 3)",
  "==================================",
  "",
  sprintf("Date: %s", Sys.Date()),
  "Strategy: translated-read protein homology (tblastn), protein query vs. raw reads",
  sprintf("Evalue threshold: %s | threads: %d | max_target_seqs: %s",
        format(RR$blast_evalue, scientific = TRUE), RR$blast_threads,
        format(RR$blast_max_target_seqs, scientific = FALSE)),
  sprintf("Additional filters applied: min_alignment_length_aa=%s, min_pident=%s",
          ifelse(is.null(RR$min_alignment_length_aa), "none", RR$min_alignment_length_aa),
          ifelse(is.null(RR$min_pident), "none", RR$min_pident)),
  sprintf("Read-ID pairing convention detected: %s (confirmed: %s)", pair_id_info$pattern, pair_id_info$confirmed),
  "",
  "Evidence layers:",
  "  Layer A  : results/intermediate/tblastn_raw/<gene>_<mate>.tsv  (untouched BLAST output)",
  "  Layer A' : raw_read_recovery_all_hsps.csv  (annotated/normalized DERIVATIVE of Layer A — not itself untouched output; *_raw columns preserve original BLAST values)",
  "  Layer B  : raw_read_recovery_per_read.csv  (deduplicated read-level evidence, query AND subject coordinates)",
  "  Layer C  : raw_read_recovery_summary.csv   (gene-level empirical summary)",
  "",
  "IMPORTANT: 'support_depth' below is TBLASTN/PROTEIN-HOMOLOGY SUPPORT DEPTH —",
  "the number of reads producing significant protein-level homology at a given",
  "protein position — NOT conventional nucleotide sequencing depth/coverage.",
  "",
  "NOTE: observation_status is a purely descriptive coverage-based label",
  "(NO_SIGNIFICANT_SUPPORT / DETECTED_PARTIAL_COVERAGE / DETECTED_FULL_LENGTH_COVERAGE).",
  "It is NOT a biological recovery/quality classification. Interpretation belongs",
  "to the PI role, pending inspection of empirical distributions across all six genes.",
  "",
  "=================================="
)
for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, ]
  report_lines <- c(report_lines,
    sprintf("%s: %s", row$gene, row$observation_status),
    sprintf("  Protein length: %d aa", row$protein_length_aa),
    sprintf("  Supporting reads: R1=%d, R2=%d, total=%d",
            row$unique_supporting_reads_R1, row$unique_supporting_reads_R2, row$unique_supporting_reads_total),
    sprintf("  Query coverage: %.2f%%", row$protein_query_coverage_pct),
    sprintf("  Support depth (mean/median/min/max): %.2f / %g / %g / %g",
            row$mean_support_depth, row$median_support_depth, row$min_support_depth, row$max_support_depth),
    sprintf("  Identity (mean/median/min/max): %s / %s / %s / %s",
            row$mean_pident, row$median_pident, row$min_pident, row$max_pident),
    sprintf("  E-value (min/median/max): %s / %s / %s", row$min_evalue, row$median_evalue, row$max_evalue),
    sprintf("  Read pairs with support (any / both mates): %s / %s (concordance: %s)",
            row$n_read_pairs_with_any_support, row$n_read_pairs_with_both_mates_supporting, row$mate_concordance_status),
    ""
  )
}
writeLines(report_lines, file.path(PATHS$reports, "raw_read_recovery_report.txt"))
log_info("Saved raw_read_recovery_report.txt.", LOG_FILE)

# ------------------------------------------------------------------------------
# 3.8 Session log copy & exit
# ------------------------------------------------------------------------------
log_system_memory(LOG_FILE)
log_info("Step 3 Raw Read Recovery completed successfully.", LOG_FILE)
file.copy(LOG_FILE, file.path(PATHS$logs, "latest_execution.log"), overwrite = TRUE)

cat(sprintf("\n[SUCCESS] Step 3 complete. Evidence generated for %d gene(s).\n", length(gene_names)))