# ==============================================================================
# Step 1: FASTQ Quality Assessment & Yield Diagnostics
# Project: Phase II Validation of Candidate Genes from T. equiperdum IVM-t1 Raw Reads
# File: scripts/01_fastq_qc.R
# ==============================================================================

cat("======================================================\n")
cat(" Phase II: Candidate Gene Recovery from Raw Reads      \n")
cat(" Step 1: FASTQ Quality Assessment (01_fastq_qc.R)      \n")
cat("======================================================\n\n")

# Load setup environment
if (!exists("PATHS") || !exists("LOG_FILE")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

# Load gridExtra for multi-plot layouts
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra", quiet = TRUE)
library(gridExtra)

# Define consistent color palette for all QC plots
QC_COLORS <- list(
  length  = "steelblue",
  gc      = "forestgreen",
  quality = "tomato",
  cycle   = "darkorange",
  n_content = "purple"
)

log_info("Starting Step 1: FASTQ Quality Assessment...", LOG_FILE)

# 2. Fast Streaming Utilities ------------------------------------------------

#' Stream-count TOTAL reads in the entire file (Memory Safe)
count_total_reads <- function(fastq_path, log_file = NULL) {
  log_info(sprintf("Counting absolute total reads in %s...", basename(fastq_path)), log_file)
  streamer <- ShortRead::FastqStreamer(fastq_path, n = 1e6)
  on.exit(close(streamer))
  total <- 0
  while (length(fq <- ShortRead::yield(streamer))) {
    total <- total + length(fq)
  }
  return(total)
}

#' Sample reads to extract deep quality metrics
sample_fastq_qc <- function(fastq_path, sample_size, log_file = NULL) {
  log_info(sprintf("Sampling %s reads from %s for deep QC metrics...", 
                   format(sample_size, big.mark = ","), basename(fastq_path)), log_file)
  sampler <- ShortRead::FastqSampler(fastq_path, n = sample_size)
  on.exit(close(sampler))
  fq <- ShortRead::yield(sampler)

  if (length(fq) == 0) stop("Failed to read FASTQ data.")

  seqs  <- ShortRead::sread(fq)
  # The quality() generic is exported by Biostrings. It dispatches the
  # appropriate method for the ShortReadQ object returned by FastqSampler.
  quals <- as(Biostrings::quality(fq), "matrix")
  
  read_lens  <- width(seqs)
  
  # Per-read metrics
  gc_counts <- rowSums(Biostrings::letterFrequency(seqs, letters = c("G", "C"), as.prob = FALSE)) / read_lens * 100
  n_counts  <- rowSums(Biostrings::letterFrequency(seqs, letters = "N", as.prob = FALSE)) / read_lens * 100
  mean_read_q <- rowMeans(quals, na.rm = TRUE)
  
  # Base composition calculated as percentages
  raw_comp <- colSums(Biostrings::letterFrequency(seqs, letters = c("A", "C", "G", "T", "N"), as.prob = FALSE))
  base_comp_pct <- round(raw_comp / sum(raw_comp) * 100, 2)
  
  # Per-cycle quality
  per_cycle_q <- colMeans(quals, na.rm = TRUE)

  list(
    sampled_count    = length(fq),
    mean_length      = mean(read_lens),
    median_length    = median(read_lens),
    sd_length        = sd(read_lens),
    min_length       = min(read_lens),
    max_length       = max(read_lens),
    unique_lengths   = paste(unique(read_lens), collapse = ", "), 
    mean_gc          = mean(gc_counts),
    median_gc        = median(gc_counts),
    sd_gc            = sd(gc_counts),
    mean_n_pct       = mean(n_counts),
    q20_pct          = mean(quals >= 20, na.rm = TRUE) * 100,
    q30_pct          = mean(quals >= 30, na.rm = TRUE) * 100,
    mean_phred       = mean(mean_read_q),
    per_cycle_q      = per_cycle_q,
    read_lengths     = read_lens,     
    gc_counts        = gc_counts,        
    mean_read_q      = mean_read_q,      
    n_counts         = n_counts,         
    base_comp_pct    = base_comp_pct    
  )
}

# 3. Process Paired FASTQ Reads ------------------------------------------------
log_info("Counting total reads...", LOG_FILE)
total_r1_reads <- count_total_reads(PATHS$fastq_r1, LOG_FILE)
total_r2_reads <- count_total_reads(PATHS$fastq_r2, LOG_FILE)

if (total_r1_reads != total_r2_reads) {
    log_warn(sprintf("Read count imbalance detected! R1: %d, R2: %d", total_r1_reads, total_r2_reads), LOG_FILE)
} else {
  log_info(sprintf("Paired-end symmetry confirmed: %s read pairs.", format(total_r1_reads, big.mark = ",")), LOG_FILE)
}

# FIX 3: Bound sample size by actual reads in file to prevent errors on small datasets
actual_sample_r1 <- min(CONFIG$qc_sample_size, total_r1_reads)
actual_sample_r2 <- min(CONFIG$qc_sample_size, total_r2_reads)

log_info("Sampling reads for quality distributions...", LOG_FILE)
qc_r1 <- sample_fastq_qc(PATHS$fastq_r1, sample_size = actual_sample_r1, log_file = LOG_FILE)
qc_r2 <- sample_fastq_qc(PATHS$fastq_r2, sample_size = actual_sample_r2, log_file = LOG_FILE)

# 4. Assemble Summary Metrics & Coverage Estimates -----------------------------
r1_bp <- total_r1_reads * qc_r1$mean_length
r2_bp <- total_r2_reads * qc_r2$mean_length
total_bp <- r1_bp + r2_bp

estimated_coverage <- round(total_bp / CONFIG$genome_size_bp, 2)

# FIX 1: Removed invalid sprintf() wrapper
log_info("Coverage Estimation Assumptions:", LOG_FILE)
log_info(sprintf("  Genome size assumed : %.2f Mb", CONFIG$genome_size_bp / 1e6), LOG_FILE)
log_info(sprintf("  Total bases yielded : %.2f Mb", total_bp / 1e6), LOG_FILE)
log_info(sprintf("  Estimated coverage  : %.2fx", estimated_coverage), LOG_FILE)

qc_summary_df <- tibble::tibble(
  read_pair           = c("R1", "R2"),
  file_name           = c(basename(PATHS$fastq_r1), basename(PATHS$fastq_r2)),
  total_reads         = c(total_r1_reads, total_r2_reads),
  unique_read_lengths = c(qc_r1$unique_lengths, qc_r2$unique_lengths),
  mean_read_len_bp    = c(round(qc_r1$mean_length, 2), round(qc_r2$mean_length, 2)),
  sd_read_len_bp      = c(round(qc_r1$sd_length, 2), round(qc_r2$sd_length, 2)),
  mean_gc_pct         = c(round(qc_r1$mean_gc, 2), round(qc_r2$mean_gc, 2)),
  median_gc_pct       = c(round(qc_r1$median_gc, 2), round(qc_r2$median_gc, 2)),
  mean_n_pct          = c(round(qc_r1$mean_n_pct, 2), round(qc_r2$mean_n_pct, 2)),
  q20_bases_pct       = c(round(qc_r1$q20_pct, 2), round(qc_r2$q20_pct, 2)),
  q30_bases_pct       = c(round(qc_r1$q30_pct, 2), round(qc_r2$q30_pct, 2)),
  mean_phred_score    = c(round(qc_r1$mean_phred, 2), round(qc_r2$mean_phred, 2)),
  throughput_mb       = c(round(r1_bp / 1e6, 2), round(r2_bp / 1e6, 2)), 
  throughput_gb       = c(round(r1_bp / 1e9, 3), round(r2_bp / 1e9, 3)),
  est_paired_coverage = estimated_coverage 
)

# 5. Generate 6-Panel Diagnostic Plots ----------------------------------------
log_info("Generating diagnostic QC plots (PDF & PNG)...", LOG_FILE)

generate_qc_plots <- function(qc_data, mate_name) {
  p_len <- ggplot2::ggplot(data.frame(length = qc_data$read_lengths), ggplot2::aes(x = length)) +
    ggplot2::geom_histogram(binwidth = 1, fill = QC_COLORS$length, color = "black", alpha = 0.7) +
    ggplot2::labs(title = paste("Read Length -", mate_name), x = "Read Length (bp)", y = "Count") +
    ggplot2::theme_minimal()
    
  p_gc <- ggplot2::ggplot(data.frame(gc = qc_data$gc_counts), ggplot2::aes(x = gc)) +
    ggplot2::geom_histogram(binwidth = 1, fill = QC_COLORS$gc, color = "black", alpha = 0.7) +
    ggplot2::labs(title = paste("GC Content -", mate_name), x = "GC (%)", y = "Count") +
    ggplot2::theme_minimal()

  p_qdist <- ggplot2::ggplot(data.frame(quality = qc_data$mean_read_q), ggplot2::aes(x = quality)) +
    ggplot2::geom_histogram(binwidth = 0.5, fill = QC_COLORS$quality, color = "black", alpha = 0.7) +
    ggplot2::labs(title = paste("Per-Read Mean Quality -", mate_name), x = "Mean Phred Score", y = "Count") +
    ggplot2::theme_minimal()

  p_n <- ggplot2::ggplot(data.frame(n = qc_data$n_counts), ggplot2::aes(x = n)) +
    ggplot2::geom_histogram(binwidth = 0.1, fill = QC_COLORS$n_content, color = "black", alpha = 0.7) +
    ggplot2::coord_cartesian(xlim = c(0, 5)) + 
    ggplot2::labs(title = paste("N Content -", mate_name), x = "N (%)", y = "Count") +
    ggplot2::theme_minimal()

  pos_df <- data.frame(position = 1:length(qc_data$per_cycle_q), mean_quality = qc_data$per_cycle_q)
  p_cycle <- ggplot2::ggplot(pos_df, ggplot2::aes(x = position, y = mean_quality)) +
    ggplot2::geom_line(color = QC_COLORS$cycle, linewidth = 1) +
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "orange") +
    ggplot2::geom_hline(yintercept = 30, linetype = "dashed", color = "green") +
    ggplot2::labs(title = paste("Per-Position Quality -", mate_name), x = "Position (bp)", y = "Mean Phred") +
    ggplot2::theme_minimal()

  comp_df <- data.frame(
    Base = c("A", "C", "G", "T", "N"),
    Percent = qc_data$base_comp_pct
  )
  # FIX 5: Fixed Y-axis to 0-40% for better run-to-run comparison
  p_comp <- ggplot2::ggplot(comp_df, ggplot2::aes(x = Base, y = Percent, fill = Base)) +
    ggplot2::geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::coord_cartesian(ylim = c(0, 40)) +
    ggplot2::labs(title = paste("Base Composition -", mate_name), x = "Nucleotide", y = "Percent (%)") +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none")

  # Return gtable object cleanly
  g <- gridExtra::arrangeGrob(p_len, p_gc, p_qdist, p_n, p_cycle, p_comp, ncol = 3, 
               top = grid::textGrob(paste("SRR7910035 QC Overview:", mate_name), 
                                     gp = grid::gpar(fontsize = 16, fontface = "bold")))
  return(g)
}

# Generate plot objects
g_r1 <- generate_qc_plots(qc_r1, "R1 (Forward)")
g_r2 <- generate_qc_plots(qc_r2, "R2 (Reverse)")

# FIX 2: Save 2-page PDF (Page 1 = R1, Page 2 = R2)
plot_file_pdf <- file.path(PATHS$reports, "fastq_qc_6panel_plots.pdf")
grDevices::pdf(plot_file_pdf, width = 14, height = 10)
grid::grid.draw(g_r1)
grid::grid.newpage()
grid::grid.draw(g_r2)
grDevices::dev.off()

# FIX 2: Save TWO separate PNG files to prevent page overwriting
grDevices::png(file.path(PATHS$reports, "fastq_qc_R1_6panel.png"), width = 14, height = 10, units = "in", res = 300)
grid::grid.draw(g_r1)
grDevices::dev.off()

grDevices::png(file.path(PATHS$reports, "fastq_qc_R2_6panel.png"), width = 14, height = 10, units = "in", res = 300)
grid::grid.draw(g_r2)
grDevices::dev.off()

log_info("Saved 2-page 6-panel QC plots (PDF) and separate mate PNGs to reports directory.", LOG_FILE)

# 7. Read-Length Assessment ----------------------------------------------------
#' Summarise read lengths using trimming-aware QC thresholds.
#'
#' Modern Illumina FASTQ files commonly contain variable read lengths after
#' adapter and quality trimming. This is generally expected and Rsubread aligns
#' variable-length reads. Warnings are therefore reserved for short reads that
#' can reduce mapping specificity or indicate an over-trimmed library.
#'
#' @param read_lengths Numeric vector of read lengths in base pairs.
#' @param read_pair_name Label used in log messages (for example, "R1").
#' @param log_file Optional pipeline log file.
#' @return A named list of reusable read-length statistics and QC flags.
assess_read_lengths <- function(read_lengths, read_pair_name, log_file = NULL) {
  read_lengths <- as.numeric(read_lengths)
  if (length(read_lengths) == 0 || anyNA(read_lengths) || any(read_lengths <= 0)) {
    stop(sprintf("Invalid read lengths supplied for %s; inspect the FASTQ input.", read_pair_name))
  }

  min_length <- min(read_lengths)
  max_length <- max(read_lengths)
  expected_length <- CONFIG$expected_read_length
  if (length(expected_length) != 1 || is.na(expected_length) || expected_length <= 0 ||
      CONFIG$warn_min_read_length <= 0 || CONFIG$warn_short50_fraction < 0 ||
      CONFIG$warn_mean_length_fraction <= 0 || CONFIG$warn_mean_length_fraction > 1 ||
      CONFIG$report_read_length_fraction <= 0 || CONFIG$report_read_length_fraction > 1) {
    stop("Invalid read-length QC settings in CONFIG.")
  }
  read_length_percentiles <- stats::quantile(
    read_lengths, probs = c(0.05, 0.25, 0.75, 0.95), names = FALSE
  )

  stats <- list(
    min_length_bp = min_length,
    max_length_bp = max_length,
    mean_length_bp = mean(read_lengths),
    median_length_bp = median(read_lengths),
    sd_length_bp = sd(read_lengths),
    p05_length_bp = unname(read_length_percentiles[1]),
    p25_length_bp = unname(read_length_percentiles[2]),
    p75_length_bp = unname(read_length_percentiles[3]),
    p95_length_bp = unname(read_length_percentiles[4]),
    pct_at_max_length = mean(read_lengths == max_length) * 100,
    pct_shorter_than_100bp = mean(read_lengths < 100) * 100,
    pct_shorter_than_50bp = mean(read_lengths < 50) * 100,
    pct_at_or_above_report_fraction = mean(read_lengths >= expected_length * CONFIG$report_read_length_fraction) * 100,
    report_read_length_fraction = CONFIG$report_read_length_fraction,
    expected_length_bp = expected_length,
    expected_sequencing_chemistry = CONFIG$expected_sequencing_chemistry,
    mean_length_reduction_pct = (1 - mean(read_lengths) / expected_length) * 100,
    variable_lengths = min_length != max_length
  )

  log_info(sprintf(
    "%s read-length statistics: min %d bp; max %d bp; mean %.2f bp; median %.2f bp; SD %.2f bp.",
    read_pair_name, stats$min_length_bp, stats$max_length_bp,
    stats$mean_length_bp, stats$median_length_bp, stats$sd_length_bp
  ), log_file)
  log_info(sprintf(
    "%s read-length distribution: %.2f%% at maximum length; %.2f%% <100 bp; %.2f%% <50 bp.",
    read_pair_name, stats$pct_at_max_length, stats$pct_shorter_than_100bp,
    stats$pct_shorter_than_50bp
  ), log_file)
  log_info(sprintf(
    "%s expected sequencing: %s; expected length %d bp; observed percentiles (P05/P25/P75/P95): %.2f/%.2f/%.2f/%.2f bp.",
    read_pair_name, stats$expected_sequencing_chemistry, stats$expected_length_bp,
    stats$p05_length_bp, stats$p25_length_bp, stats$p75_length_bp, stats$p95_length_bp
  ), log_file)

  if (stats$variable_lengths) {
    log_info(sprintf(
      "%s: Variable read lengths detected (%d-%d bp). This pattern is consistent with adapter/quality trimming and is generally compatible with Rsubread alignment.",
      read_pair_name, stats$min_length_bp, stats$max_length_bp
    ), log_file)
  } else {
    log_info(sprintf("%s: All sampled reads are %d bp.", read_pair_name, stats$max_length_bp), log_file)
  }

  qc_flags <- character()

  # Reads below the configured threshold often lack enough unique sequence context for reliable
  # genomic placement, even though their presence alone does not imply corruption.
  if (stats$min_length_bp < CONFIG$warn_min_read_length) {
    qc_flags <- c(qc_flags, sprintf(
      "minimum read length (%d bp) is below the %d bp mapping-specificity threshold",
      stats$min_length_bp, CONFIG$warn_min_read_length
    ))
  }

  # The configured fraction of reads below 50 bp identifies libraries that are
  # commonly over-trimmed and can lose mapping sensitivity; this is more
  # informative than length variation.
  if (stats$pct_shorter_than_50bp > CONFIG$warn_short50_fraction * 100) {
    qc_flags <- c(qc_flags, sprintf(
      "%.2f%% of reads are shorter than 50 bp (threshold %.2f%%)",
      stats$pct_shorter_than_50bp, CONFIG$warn_short50_fraction * 100
    ))
  }

  # The retained mean-length fraction is configurable for platforms with 100,
  # 150, 250, or 301 bp chemistry. A lower value suggests substantial trimming
  # that warrants review before mapping.
  if (stats$mean_length_bp < expected_length * CONFIG$warn_mean_length_fraction) {
    qc_flags <- c(qc_flags, sprintf(
      "mean read length (%.2f bp) is below %.0f%% of the expected %d bp",
      stats$mean_length_bp, CONFIG$warn_mean_length_fraction * 100, expected_length
    ))
  }

  stats$qc_status <- if (length(qc_flags) == 0) "PASS" else "WARNING"
  stats$qc_interpretation <- if (length(qc_flags) == 0) {
    "Variable read lengths are consistent with adapter/quality trimming; no evidence of excessive trimming."
  } else {
    paste("Review recommended:", paste(qc_flags, collapse = "; "), ".")
  }
  log_info(sprintf("%s read-length assessment: %s. %s", read_pair_name,
                   stats$qc_status, stats$qc_interpretation), log_file)
  if (length(qc_flags) > 0) {
    log_warn(sprintf("%s read-length QC warning: %s", read_pair_name,
                     paste(qc_flags, collapse = "; ")), log_file)
  }

  stats
}

# Store the reusable assessment in the QC objects so it is retained in the RDS
# output and available to downstream scripts without parsing text log messages.
qc_r1$read_length_assessment <- assess_read_lengths(qc_r1$read_lengths, "R1", log_file = LOG_FILE)
qc_r2$read_length_assessment <- assess_read_lengths(qc_r2$read_lengths, "R2", log_file = LOG_FILE)

# Add the structured assessment to the machine-readable summary. The complete
# assessment lists are also retained in the RDS object below for downstream use.
qc_summary_df <- qc_summary_df |>
  dplyr::mutate(
    min_read_len_bp = c(qc_r1$read_length_assessment$min_length_bp, qc_r2$read_length_assessment$min_length_bp),
    max_read_len_bp = c(qc_r1$read_length_assessment$max_length_bp, qc_r2$read_length_assessment$max_length_bp),
    median_read_len_bp = c(qc_r1$read_length_assessment$median_length_bp, qc_r2$read_length_assessment$median_length_bp),
    p05_read_len_bp = c(qc_r1$read_length_assessment$p05_length_bp, qc_r2$read_length_assessment$p05_length_bp),
    p25_read_len_bp = c(qc_r1$read_length_assessment$p25_length_bp, qc_r2$read_length_assessment$p25_length_bp),
    p75_read_len_bp = c(qc_r1$read_length_assessment$p75_length_bp, qc_r2$read_length_assessment$p75_length_bp),
    p95_read_len_bp = c(qc_r1$read_length_assessment$p95_length_bp, qc_r2$read_length_assessment$p95_length_bp),
    pct_at_max_read_len = c(qc_r1$read_length_assessment$pct_at_max_length, qc_r2$read_length_assessment$pct_at_max_length),
    pct_reads_lt_100bp = c(qc_r1$read_length_assessment$pct_shorter_than_100bp, qc_r2$read_length_assessment$pct_shorter_than_100bp),
    pct_reads_lt_50bp = c(qc_r1$read_length_assessment$pct_shorter_than_50bp, qc_r2$read_length_assessment$pct_shorter_than_50bp),
    pct_reads_ge_report_length = c(qc_r1$read_length_assessment$pct_at_or_above_report_fraction, qc_r2$read_length_assessment$pct_at_or_above_report_fraction),
    expected_read_len_bp = c(qc_r1$read_length_assessment$expected_length_bp, qc_r2$read_length_assessment$expected_length_bp),
    expected_sequencing_chemistry = c(qc_r1$read_length_assessment$expected_sequencing_chemistry, qc_r2$read_length_assessment$expected_sequencing_chemistry),
    mean_length_reduction_pct = c(qc_r1$read_length_assessment$mean_length_reduction_pct, qc_r2$read_length_assessment$mean_length_reduction_pct),
    read_length_qc_status = c(qc_r1$read_length_assessment$qc_status, qc_r2$read_length_assessment$qc_status),
    read_length_qc_interpretation = c(qc_r1$read_length_assessment$qc_interpretation, qc_r2$read_length_assessment$qc_interpretation)
  )
write.csv(qc_summary_df, file.path(PATHS$reports, "fastq_qc_summary.csv"), row.names = FALSE)

# Save the final QC objects only after read-length assessments are attached.
saveRDS(list(R1 = qc_r1, R2 = qc_r2, summary = qc_summary_df),
        file = file.path(PATHS$reports, "fastq_qc_raw_metrics.rds"))

# 8. Generate Text Summary Report ---------------------------------------------
report_txt <- file.path(PATHS$reports, "fastq_qc_report.txt")
format_read_length_report <- function(assessment, read_pair_name) {
  c(
    sprintf("%s read-length assessment: %s", read_pair_name, assessment$qc_status),
    sprintf("Expected sequencing: %s", assessment$expected_sequencing_chemistry),
    sprintf("Expected length: %d bp; observed mean: %.2f bp; observed median: %.2f bp.",
            assessment$expected_length_bp, assessment$mean_length_bp, assessment$median_length_bp),
    sprintf("Range/SD (bp): min %d; max %d; SD %.2f; P05/P25/P75/P95: %.2f/%.2f/%.2f/%.2f.",
            assessment$min_length_bp, assessment$max_length_bp, assessment$sd_length_bp,
            assessment$p05_length_bp, assessment$p25_length_bp,
            assessment$p75_length_bp, assessment$p95_length_bp),
    sprintf("Distribution: %.2f%% at max; %.2f%% >=%.0f%% expected length; %.2f%% <100 bp; %.2f%% <50 bp.",
            assessment$pct_at_max_length,
            assessment$pct_at_or_above_report_fraction,
            assessment$report_read_length_fraction * 100,
            assessment$pct_shorter_than_100bp, assessment$pct_shorter_than_50bp),
    assessment$qc_interpretation
  )
}
report_lines <- c(
  "=====================================",
  "FASTQ QUALITY REPORT",
  "=====================================",
  "",
  sprintf("Read pairs: %s", format(total_r1_reads, big.mark = ",")),
  "Read-length assessment:",
  format_read_length_report(qc_r1$read_length_assessment, "R1"),
  format_read_length_report(qc_r2$read_length_assessment, "R2"),
  "",
  sprintf("Mean Quality (R1/R2): %.1f / %.1f", qc_r1$mean_phred, qc_r2$mean_phred),
  sprintf("Q30 Bases (R1/R2): %.1f%% / %.1f%%", qc_r1$q30_pct, qc_r2$q30_pct),
  sprintf("Mean GC (R1/R2): %.1f%% / %.1f%%", qc_r1$mean_gc, qc_r2$mean_gc),
  sprintf("Mean N%% (R1/R2): %.2f%% / %.2f%%", qc_r1$mean_n_pct, qc_r2$mean_n_pct),
  "",
  sprintf("Throughput: %.2f Gb (%.2f Mb)", total_bp / 1e9, total_bp / 1e6),
  sprintf("Genome assumption: %.2f Mb", CONFIG$genome_size_bp / 1e6),
  sprintf("Estimated coverage: %.2fx", estimated_coverage),
  "",
  "====================================="
)
writeLines(report_lines, report_txt)
log_info("Saved human-readable text report to reports/fastq_qc_report.txt", LOG_FILE)

# 9. Session Log Copy ------------------------------------------------------
log_system_memory(LOG_FILE)
log_info("Step 1 FASTQ Quality Assessment completed successfully.", LOG_FILE)
file.copy(LOG_FILE, file.path(PATHS$logs, "latest_execution.log"), overwrite = TRUE)

cat("\n[SUCCESS] Step 1 finished. Exact read counts, coverage estimates, base composition, and 6-panel QC plots generated.\n")
