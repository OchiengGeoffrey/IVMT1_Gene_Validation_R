source("scripts/01_project_setup.R")

project_header("Extract ATPase gamma CDS")

###############################################################
## Path definitions
###############################################################

report_file   <- file.path(PATHS$reports, "ATPF1G_tblastn_results.csv")
cds_out_file  <- file.path(PATHS$extracted, "ATPF1G_CDS.fasta")
prot_out_file <- file.path(PATHS$extracted, "ATPF1G_protein.fasta")
ref_file      <- PATHS$atpf1g

###############################################################
## Create output directory
###############################################################

dir.create(PATHS$extracted, recursive = TRUE, showWarnings = FALSE)

###############################################################
## Load BLAST results
###############################################################

cat("Loading tblastn results...\n")
hits <- read.csv(report_file, stringsAsFactors = FALSE)

if (nrow(hits) == 0) {
  stop("No BLAST hits found in: ", report_file)
}

###############################################################
## Best hit
###############################################################

best_hit <- hits[which.max(hits$bitscore), ]

contig      <- best_hit$sseqid
start_coord <- best_hit$sstart
end_coord   <- best_hit$send

cat(sprintf("Best hit on contig %s: %d to %d (Bitscore: %g)\n\n", 
            contig, start_coord, end_coord, best_hit$bitscore))

###############################################################
## Load assembly
###############################################################

cat("Loading genome assembly...\n")
assembly <- read_dna_fasta(PATHS$assembly)

###############################################################
## Extract CDS
###############################################################

cat("Extracting coding sequence...\n")

cds_dna <- extract_region(
  assembly = assembly,
  contig   = contig,
  start    = start_coord,
  end      = end_coord,
  flank    = 0
)

strand_char <- if (start_coord < end_coord) "+" else "-"
names(cds_dna) <- sprintf("ATPF1G_CDS | %s:%d-%d | strand:%s", contig, start_coord, end_coord, strand_char)

cds_len <- Biostrings::width(cds_dna)[1]

if (cds_len %% 3 != 0) {
  stop("CDS length is not divisible by 3.")
}

if (cds_len != 915) {
  stop(sprintf("CDS length should be 915 bp but observed %d bp.", cds_len))
}

###############################################################
## Start codon check
###############################################################

start_codon <- substr(as.character(cds_dna[[1]]), 1, 3)

if (start_codon != "ATG") {
  warning(sprintf("Expected start codon ATG but found %s.", start_codon))
}

###############################################################
## Translation
###############################################################

cat("Translating CDS to protein...\n")

protein_seq <- translate_cds(cds_dna)

if (any(is.na(protein_seq))) {
  stop("Protein translation failed.")
}

names(protein_seq) <- sprintf("ATPF1G_protein | %s:%d-%d", contig, start_coord, end_coord)

print_sequence_info(protein_seq)

prot_len <- Biostrings::width(protein_seq)[1]

if (prot_len != 305) {
  stop(sprintf("Protein length should be 305 aa but observed %d aa.", prot_len))
}

###############################################################
## Compare against reference
###############################################################

identity_val <- NA

if (file.exists(ref_file)) {
  cat("\nComparing against reference protein...\n")
  reference <- read_protein_fasta(ref_file)
  identity_val <- round(percent_identity(protein_seq, reference), 2)
  cat("Protein identity :", identity_val, "%\n\n")
}

###############################################################
## Automatic validation
###############################################################

project_header("Automatic validation")

prot_str   <- as.character(protein_seq[[1]])
prot_chars <- unlist(strsplit(prot_str, ""))

internal_stop_codons <- sum(prot_chars == "*")
frame_valid          <- (cds_len %% 3 == 0)

cat("Expected CDS length   : 915 bp\n")
cat("Observed CDS length   :", cds_len, "bp\n")
cat("Expected protein      : 305 aa\n")
cat("Observed protein      :", prot_len, "aa\n")
cat("Start codon           :", start_codon, "\n")
cat("Internal stop codons  :", internal_stop_codons, "\n")
cat("Frame valid           :", frame_valid, "\n\n")

###############################################################
## Save outputs
###############################################################

save_fasta(cds_dna, cds_out_file)
save_fasta(protein_seq, prot_out_file)

###############################################################
## Final report
###############################################################

cat("ATPase gamma CDS extracted successfully.\n\n")
cat("Contig:\n")
cat(contig, "\n\n")
cat("Coordinates:\n")
cat(sprintf("%d-%d", start_coord, end_coord), "\n\n")
cat("CDS length:\n")
cat(cds_len, "bp\n\n")
cat("Protein length:\n")
cat(prot_len, "aa\n\n")

if (!is.na(identity_val)) {
  cat("Protein identity:\n")
  cat(identity_val, "%\n\n")
}

cat("Results written to\n")
cat(cds_out_file, "\n")
cat(prot_out_file, "\n")