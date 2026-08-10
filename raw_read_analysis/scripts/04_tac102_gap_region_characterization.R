# ==============================================================================
# Step 4B.1: TAC102 Reference Region Characterization (residues 649-732)
# Project: Genomic & Transcriptomic Limits on kDNA Retention in T. equiperdum
#
# Objective:
#   Characterize the TREU927 TAC102 reference sequence at the raw-read
#   zero-support interval (649-732 aa) identified in Step 4A:
#     1. Extract the exact subsequence
#     2. Amino-acid composition
#     3. Low-complexity / repetitiveness screen
#     4. Flag for manual domain/conserved-region lookup (no automated DB call
#        yet - TriTrypDB/InterPro cross-reference is a separate, deliberate
#        step, not bundled in here)
#
# This step does NOT touch raw reads. It is independent of Step 4B's
# nucleotide-evidence analysis by design (per PI: parallel, not sequential,
# to avoid biasing composition/domain interpretation with nucleotide results).
#
# Input:
#   raw_read_analysis/results/intermediate/reference_sequences.rds
#   raw_read_analysis/reports/tac102_uncovered_intervals.csv
#
# Output:
#   tac102_gap_region_sequence.fasta
#   tac102_gap_region_composition.csv
#   tac102_gap_region_report.txt
# ==============================================================================

cat("======================================================\n")
cat(" Step 4B.1: TAC102 Reference Region Characterization   \n")
cat("======================================================\n\n")

if (!exists("PATHS") || !exists("LOG_FILE") || !exists("CONFIG")) {
  source(here::here("raw_read_analysis", "scripts", "00_setup.R"))
}

log_info("Starting Step 4B.1: TAC102 gap-region reference characterization.", LOG_FILE)

# ------------------------------------------------------------------------------
# 1. Load inputs
# ------------------------------------------------------------------------------
ref_seqs_path <- file.path(PATHS$intermediate, "reference_sequences.rds")
if (!file.exists(ref_seqs_path)) {
  log_error(sprintf("Missing Step 2 output: %s", ref_seqs_path), LOG_FILE)
  stop("Run 02_reference_prep.R before Step 4B.1.")
}
reference_sequences <- readRDS(ref_seqs_path)

gap_intervals_path <- file.path(PATHS$reports, "tac102_uncovered_intervals.csv")
if (!file.exists(gap_intervals_path)) {
  log_error(sprintf("Missing Step 4A output: %s", gap_intervals_path), LOG_FILE)
  stop("Run 04_tac102_gap_characterization.R before Step 4B.1.")
}
gap_df <- read.csv(gap_intervals_path, stringsAsFactors = FALSE)

if (nrow(gap_df) != 1) {
  log_error(
    sprintf("Expected exactly 1 gap interval for TAC102, found %d. Aborting.", nrow(gap_df)),
    LOG_FILE
  )
  stop("Execution halted: unexpected number of gap intervals.")
}

gap_start <- gap_df$start_aa[1]
gap_end   <- gap_df$end_aa[1]
gap_len   <- gap_df$length_aa[1]

log_info(
  sprintf("Loaded gap interval from Step 4A: %d-%d aa (%d aa).", gap_start, gap_end, gap_len),
  LOG_FILE
)

if (!"TAC102" %in% names(reference_sequences)) {
  log_error("TAC102 not found in reference_sequences.rds.", LOG_FILE)
  stop("Execution halted: TAC102 reference missing.")
}
tac102_full <- reference_sequences[["TAC102"]]
full_len <- length(tac102_full)

if (gap_end > full_len) {
  log_error(
    sprintf("Gap end (%d) exceeds reference length (%d).", gap_end, full_len),
    LOG_FILE
  )
  stop("Execution halted: gap interval inconsistent with reference length.")
}

# ------------------------------------------------------------------------------
# 2. Extract the gap-region subsequence
# ------------------------------------------------------------------------------
gap_region <- Biostrings::subseq(tac102_full, start = gap_start, end = gap_end)
gap_region_chr <- as.character(gap_region)

log_info(
  sprintf("Extracted TAC102 residues %d-%d (%d aa).", gap_start, gap_end, nchar(gap_region_chr)),
  LOG_FILE
)

# Save as FASTA for downstream inspection
fasta_out <- file.path(PATHS$reports, "tac102_gap_region_sequence.fasta")
gap_aaset <- Biostrings::AAStringSet(gap_region_chr)
names(gap_aaset) <- sprintf("TAC102_%d-%d", gap_start, gap_end)
Biostrings::writeXStringSet(gap_aaset, filepath = fasta_out, format = "fasta")
log_info(sprintf("Saved gap-region FASTA: %s", basename(fasta_out)), LOG_FILE)

# ------------------------------------------------------------------------------
# 3. Amino-acid composition (gap region vs. whole-protein baseline)
# ------------------------------------------------------------------------------
full_chr <- as.character(tac102_full)

aa_freq <- function(seq_chr) {
  chars <- base::strsplit(as.character(seq_chr), "")[[1]]
  tab <- table(factor(chars, levels = Biostrings::AA_STANDARD))
  round(as.numeric(tab) / length(chars) * 100, 2)
}

gap_comp  <- aa_freq(gap_region_chr)
full_comp <- aa_freq(full_chr)

composition_df <- data.frame(
  residue          = Biostrings::AA_STANDARD,
  gap_region_pct   = gap_comp,
  full_protein_pct = full_comp,
  enrichment_ratio = round(ifelse(full_comp == 0, NA, gap_comp / full_comp), 2),
  stringsAsFactors = FALSE
)
composition_df <- composition_df[order(-composition_df$enrichment_ratio), ]

comp_out <- file.path(PATHS$reports, "tac102_gap_region_composition.csv")
write.csv(composition_df, comp_out, row.names = FALSE)
log_info(sprintf("Saved amino-acid composition table: %s", basename(comp_out)), LOG_FILE)

# ------------------------------------------------------------------------------
# 4. Low-complexity / repetitiveness screen
# ------------------------------------------------------------------------------
# 4a. Single-residue dominance threshold
dominant_residues <- composition_df[
  !is.na(composition_df$enrichment_ratio) &
  composition_df$enrichment_ratio >= 2 &
  composition_df$gap_region_pct >= 15,
]

# 4b. Shannon entropy (base::strsplit protected)
shannon_entropy <- function(seq_chr) {
  chars <- base::strsplit(as.character(seq_chr), "")[[1]]
  p <- table(chars) / length(chars)
  -sum(p * log2(p))
}
gap_entropy  <- round(shannon_entropy(gap_region_chr), 3)
full_entropy <- round(shannon_entropy(full_chr), 3)
max_entropy  <- round(log2(20), 3)

# 4c. Simple short tandem-repeat scan
find_simple_repeats <- function(seq_chr, unit_size, min_repeats = 3) {
  n <- nchar(seq_chr)
  hits <- character(0)
  i <- 1
  while (i <= n - unit_size * min_repeats + 1) {
    unit <- substr(seq_chr, i, i + unit_size - 1)
    repeat_count <- 1
    j <- i + unit_size
    while (j <= n - unit_size + 1 && substr(seq_chr, j, j + unit_size - 1) == unit) {
      repeat_count <- repeat_count + 1
      j <- j + unit_size
    }
    if (repeat_count >= min_repeats) {
      hits <- c(hits, sprintf("%s x%d at gap-relative pos %d-%d", unit, repeat_count, i, j - 1))
      i <- j
    } else {
      i <- i + 1
    }
  }
  hits
}
repeat_hits_2 <- find_simple_repeats(gap_region_chr, unit_size = 2)
repeat_hits_3 <- find_simple_repeats(gap_region_chr, unit_size = 3)

# ------------------------------------------------------------------------------
# 5. Report generation
# ------------------------------------------------------------------------------
report <- c(
  "==============================================",
  "TAC102 GAP-REGION REFERENCE CHARACTERIZATION",
  "==============================================",
  "",
  sprintf("Date: %s", Sys.Date()),
  sprintf("Gene: TAC102 (reference: TREU927, 951 aa)"),
  sprintf("Gap region (from Step 4A): residues %d-%d (%d aa)", gap_start, gap_end, gap_len),
  "",
  "IMPORTANT SCOPE NOTE:",
  "This step characterizes the REFERENCE sequence at the gap coordinates only.",
  "It does NOT establish why tBLASTn failed to detect homology there - that",
  "requires the separate nucleotide-evidence analysis (Step 4B.3).",
  "",
  "-- Sequence --",
  gap_region_chr,
  "",
  "-- Compositional bias screen --",
  sprintf("Shannon entropy, gap region:    %.3f bits (max possible: %.3f)", gap_entropy, max_entropy),
  sprintf("Shannon entropy, whole protein: %.3f bits (max possible: %.3f)", full_entropy, max_entropy),
  sprintf("Entropy reduction in gap region: %.3f bits (%.1f%% of whole-protein entropy)",
          full_entropy - gap_entropy, (gap_entropy / full_entropy) * 100),
  "",
  if (nrow(dominant_residues) > 0) {
    c("Residues enriched >=2x relative to whole-protein frequency AND >=15% of gap region:",
      sprintf("  %s: %.2f%% of gap region (vs %.2f%% whole-protein, %.2fx enrichment)",
              dominant_residues$residue, dominant_residues$gap_region_pct,
              dominant_residues$full_protein_pct, dominant_residues$enrichment_ratio))
  } else {
    "No single residue meets the >=2x enrichment AND >=15% dominance threshold."
  },
  "",
  "Simple tandem-repeat scan (descriptive only):",
  if (length(repeat_hits_2) > 0) c("  Dipeptide repeats:", paste0("    ", repeat_hits_2)) else "  No dipeptide repeats (>=3x) found.",
  if (length(repeat_hits_3) > 0) c("  Tripeptide repeats:", paste0("    ", repeat_hits_3)) else "  No tripeptide repeats (>=3x) found.",
  "",
  "-- Full amino-acid composition table --",
  "See: tac102_gap_region_composition.csv",
  "=============================================="
)

report_out <- file.path(PATHS$reports, "tac102_gap_region_report.txt")
writeLines(report, report_out)
log_info(sprintf("Saved gap-region characterization report: %s", basename(report_out)), LOG_FILE)

# ------------------------------------------------------------------------------
# 6. Exit
# ------------------------------------------------------------------------------
log_info("Step 4B.1 TAC102 reference-region characterization completed.", LOG_FILE)
file.copy(LOG_FILE, file.path(PATHS$logs, "latest_execution.log"), overwrite = TRUE)

cat("\n[SUCCESS] Step 4B.1 completed.\n")
cat(sprintf("  Gap-region FASTA: %s\n", fasta_out))
cat(sprintf("  Composition table: %s\n", comp_out))
cat(sprintf("  Report: %s\n", report_out))