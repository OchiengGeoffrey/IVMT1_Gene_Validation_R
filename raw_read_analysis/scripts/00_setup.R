# ==============================================================================
# Step 0: Project Setup and Environment Initialization
# Project: Phase II Validation of Candidate Genes from T. equiperdum IVM-t1 Raw Reads
# File: scripts/00_setup.R
# ==============================================================================

# 1. Startup Banner ------------------------------------------------------------
cat("======================================================\n")
cat(" Phase II: Candidate Gene Recovery from Raw Reads      \n")
cat(" Organism: Trypanosoma equiperdum IVM-t1               \n")
cat(" Accession: SRR7910035                                 \n")
cat(" Pipeline Infrastructure: Setup (00_setup.R)           \n")
cat("======================================================\n\n")

# 2. Package Dependency Management & Explicit Version Recording ---------------
cran_pkgs <- c("here", "tidyverse", "seqinr", "gridExtra")
bioc_pkgs <- c("Biostrings", "ShortRead", "Rsamtools", "GenomicAlignments", 
               "GenomicRanges", "IRanges", "Rsubread", "rtracklayer")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, quiet = TRUE)
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
}

all_pkgs <- c(cran_pkgs, bioc_pkgs)
suppressPackageStartupMessages({
  lapply(all_pkgs, library, character.only = TRUE)
})

# 3. Absolute Path Definitions -------------------------------------------------
PATHS <- list(
  root      = here::here("raw_read_analysis"),
  reads     = here::here("data", "reads"),
  # Protein FASTA references are shared with the established project pipeline.
  # The annotation manifest remains in raw_read_analysis/reference (below).
  ref       = here::here("data", "reference"),
  mapped    = here::here("raw_read_analysis", "results", "mapped_reads"),
  consensus = here::here("raw_read_analysis", "results", "consensus"),
  proteins  = here::here("raw_read_analysis", "results", "proteins"),
  reports   = here::here("raw_read_analysis", "reports"),
  scripts   = here::here("raw_read_analysis", "scripts"),
  logs      = here::here("raw_read_analysis", "logs")
)

# Directory Creation
lapply(PATHS[c("ref", "mapped", "consensus", "proteins", "reports", "scripts", "logs")], function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
})

# 4. Timestamped Session Logging Setup -----------------------------------------
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE  <- file.path(PATHS$logs, sprintf("execution_%s.log", timestamp))
cat("", file = LOG_FILE)

# Source Helper Functions & Set Seed
source(file.path(PATHS$scripts, "helpers.R"))
set.seed(42)

log_info("Initialized Phase II project setup", LOG_FILE)

# 5. Project & Reference Metadata Versioning -----------------------------------
PROJECT <- list(
  name              = "IVM-t1 Raw Read Candidate Gene Validation",
  species           = "Trypanosoma equiperdum",
  strain            = "IVM-t1",
  accession         = "SRR7910035",
  date              = Sys.Date(),
  pipeline_version  = "1.0",
  reference_version = "1.0"
)

# Documented Pipeline Settings
CONFIG <- list(
  threads              = 8,      # Maximum CPU cores allocated for parallel operations (Rsubread/samtools)
  min_mapping_quality  = 20,     # Minimum MAPQ threshold to exclude ambiguous alignments
  min_base_quality     = 20,     # Minimum Phred quality score for consensus base calling
  consensus_fraction   = 0.75,   # Minimum allele frequency threshold required to call a clear consensus base
  min_depth            = 10,      # Minimum read depth required to avoid calling missing/ambiguous bases (N)
  genome_size_bp       = 35e6,    # NEW: Approx T. equiperdum nuclear genome
  qc_sample_size       = 1e5,     # Number of reads per mate sampled for QC
  expected_sequencing_chemistry = "Illumina MiSeq PE301",
  expected_read_length = 301,     # Instrument read length before adapter/quality trimming
  warn_min_read_length = 30,      # bp: shorter reads often map ambiguously
  warn_short50_fraction = 0.10,   # fraction: flags potentially over-trimmed libraries
  warn_mean_length_fraction = 0.80, # fraction of expected length retained on average
  report_read_length_fraction = 0.80 # report fraction of reads retaining this expected length
)

# Record explicit package versions
installed_df <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
pkg_versions <- installed_df[installed_df$Package %in% all_pkgs, ]
write.csv(pkg_versions, file = file.path(PATHS$reports, "package_versions.csv"), row.names = FALSE)
log_info("Saved package_versions.csv to reports directory", LOG_FILE)

# Detect input FASTQ format (.fastq or .fastq.gz)
PATHS$fastq_r1 <- detect_fastq(here::here("data", "reads", "SRR7910035_1"), LOG_FILE)
PATHS$fastq_r2 <- detect_fastq(here::here("data", "reads", "SRR7910035_2"), LOG_FILE)

# 6. Candidate Gene Metadata Table (Robust CSV Parser & Validator) ------------
annotation_csv <- here::here("raw_read_analysis", "reference", "gene_annotations.csv")

if (!file.exists(annotation_csv)) {
  log_error(sprintf("Missing annotation metadata CSV at: %s", annotation_csv), LOG_FILE)
  stop("Execution halted: Missing reference/gene_annotations.csv")
}

# Auto-detect delimiter
first_line <- readLines(annotation_csv, n = 1)

if (grepl(";", first_line)) {
  TARGETS <- read.csv(annotation_csv, stringsAsFactors = FALSE, sep = ";")
  log_info("Detected semicolon (;) delimiter in gene_annotations.csv", LOG_FILE)
} else if (grepl("\t", first_line)) {
  TARGETS <- read.delim(annotation_csv, stringsAsFactors = FALSE, sep = "\t")
  log_info("Detected tab (\\t) delimiter in gene_annotations.csv", LOG_FILE)
} else {
  TARGETS <- read.csv(annotation_csv, stringsAsFactors = FALSE)
  log_info("Detected comma (,) delimiter in gene_annotations.csv", LOG_FILE)
}

# Diagnostics & Column Assertion
cat("\n[DIAGNOSTIC] Column names imported from CSV:\n")
print(colnames(TARGETS))

if (!"reference_file" %in% colnames(TARGETS)) {
  msg <- paste("Column 'reference_file' not found. Imported columns are:", paste(colnames(TARGETS), collapse = ", "))
  log_error(msg, LOG_FILE)
  stop(msg)
}

TARGETS <- TARGETS |>
  dplyr::mutate(reference_file = file.path(PATHS$ref, reference_file))

# Save TARGETS table in RDS and CSV formats
saveRDS(TARGETS, file = file.path(PATHS$reports, "target_genes_metadata.rds"))
write.csv(TARGETS, file = file.path(PATHS$reports, "target_genes_metadata.csv"), row.names = FALSE)
log_info("Successfully loaded gene_annotations.csv and exported target_genes_metadata (RDS and CSV)", LOG_FILE)

# 7. Verification & Environment Snapshot ---------------------------------------
verify_write_permissions(PATHS, LOG_FILE)
verify_inputs_and_fastas(PATHS, TARGETS, LOG_FILE)

# Check software dependencies
check_program("samtools", required = FALSE, log_file = LOG_FILE)
check_program("diamond", required = FALSE, log_file = LOG_FILE)

# Diagnostic & Session Info Export
log_system_memory(LOG_FILE)
capture.output(sessionInfo(), file = file.path(PATHS$reports, "sessionInfo.txt"))
log_info("Session info saved to reports/sessionInfo.txt", LOG_FILE)
log_info("Step 0 Project Setup completed successfully.", LOG_FILE)

# Replicate complete session log to latest_execution.log
file.copy(LOG_FILE, file.path(PATHS$logs, "latest_execution.log"), overwrite = TRUE)

cat("\n[SUCCESS] Setup finished. Infrastructure, logging, package versions, and targets table initialized.\n")
