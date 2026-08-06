source("scripts/01_project_setup.R")

project_header("Verify ATOM40 Presence")

###############################################################
## Input files
###############################################################

reference_file <- PATHS$atom40
assembly_file  <- PATHS$assembly
assembly_db    <- file.path(PATHS$blast, "IVMT1")

###############################################################
## Output files
###############################################################

blast_csv      <- file.path(PATHS$reports, "ATOM40_tblastn_results.csv")
genomic_fasta  <- file.path(PATHS$extracted, "ATOM40_genomic.fasta")
protein_fasta  <- file.path(PATHS$extracted, "ATOM40_protein.fasta")

dir.create(PATHS$reports, recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$extracted, recursive = TRUE, showWarnings = FALSE)

###############################################################
## Step 1. Search assembly
###############################################################

cat("Searching for ATOM40 in the IVM-t1 genome assembly...\n\n")

run_tblastn(
    query    = reference_file,
    database = assembly_db,
    output   = blast_csv
)

###############################################################
## Step 2. Read BLAST results
###############################################################

blast <- read.csv(
    blast_csv,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
)

if (nrow(blast) == 0) {
    stop("No BLAST hits were recovered for ATOM40.")
}

colnames(blast) <- c(
    "qseqid","sseqid","pident","length",
    "mismatch","gapopen","qstart","qend",
    "sstart","send","evalue","bitscore"
)

best_hit <- blast[1, ]

cat("Top tblastn hit:\n")
print(best_hit[, c(
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "evalue",
    "bitscore"
)])

###############################################################
## Step 3. Extract genomic region
###############################################################

extract_tblastn_hit(
    assembly_fasta = assembly_file,
    blast_row      = best_hit,
    outfile        = genomic_fasta
)

###############################################################
## Step 4. Translate genomic fragment
###############################################################

translate_genomic_sequence(
    genomic_fasta = genomic_fasta,
    protein_fasta = protein_fasta
)

cat("\nATOM40 genomic fragment written to:\n")
cat(genomic_fasta, "\n")

cat("\nTranslated protein written to:\n")
cat(protein_fasta, "\n")