source("scripts/01_project_setup.R")

project_header("Compare ATPase alpha proteins")

###############################################################
## Path definitions
###############################################################

ref_file  <- PATHS$atpf1a
qry_file  <- file.path(PATHS$extracted, "ATPF1A_protein.fasta")
blast_csv <- file.path(PATHS$reports, "ATPF1A_tblastn_results.csv")
out_csv   <- file.path(PATHS$reports, "ATPF1A_amino_acid_differences.csv")

###############################################################
## Load data & BLAST coordinates
###############################################################

cat("Loading protein sequences and BLAST coordinates...\n\n")

if (!file.exists(ref_file)) {
  stop("Reference protein file not found: ", ref_file)
}

if (!file.exists(qry_file)) {
  stop("Query protein file not found: ", qry_file)
}

if (!file.exists(blast_csv)) {
  stop("BLAST results file not found: ", blast_csv)
}

reference <- read_protein_fasta(ref_file)
ivmt1     <- read_protein_fasta(qry_file)
blast_res <- read.csv(blast_csv, stringsAsFactors = FALSE)

# Extract reference window from best tblastn hit
best_hit  <- blast_res[which.max(blast_res$bitscore), ]
ref_start <- best_hit$qstart
ref_end   <- best_hit$qend

ref_chars_full <- unlist(strsplit(as.character(reference[[1]]), ""))
qry_chars      <- unlist(strsplit(as.character(ivmt1[[1]]), ""))

ref_len <- length(ref_chars_full)
qry_len <- length(qry_chars)

###############################################################
## Coordinate-aware partial sequence comparison
###############################################################

# Slice reference sequence to match BLAST alignment coordinates
ref_chars_sub    <- ref_chars_full[ref_start:ref_end]
alignment_length <- length(ref_chars_sub)
eval_len         <- min(alignment_length, qry_len)

if (qry_len != alignment_length) {
  warning(
    sprintf(
      "Recovered protein (%d aa) does not match the BLAST alignment length (%d aa). Comparison truncated to %d residues.",
      qry_len,
      alignment_length,
      eval_len
    )
  )
}

ref_eval <- ref_chars_sub[1:eval_len]
qry_eval <- qry_chars[1:eval_len]

diff_pos_relative <- which(ref_eval != qry_eval)
diff_pos_absolute <- diff_pos_relative + ref_start - 1

mutation_table <- data.frame(
  Position  = diff_pos_absolute,
  Mutation  = paste0(ref_eval[diff_pos_relative], diff_pos_absolute, qry_eval[diff_pos_relative]),
  Reference = ref_eval[diff_pos_relative],
  IVM_t1    = qry_eval[diff_pos_relative],
  stringsAsFactors = FALSE
)

identity_val <- round(100 * (eval_len - length(diff_pos_relative)) / eval_len, 2)

###############################################################
## Save report
###############################################################

dir.create(PATHS$reports, recursive = TRUE, showWarnings = FALSE)
save_blast_results(mutation_table, out_csv)

###############################################################
## Console report
###############################################################

cat(sprintf("Reference length         : %d aa\n", ref_len))
cat(sprintf("Recovered fragment       : %d aa\n\n", qry_len))

cat(sprintf("Aligned reference region : aa %d-%d\n", ref_start, ref_end))
cat(sprintf("BLAST alignment length   : %d aa\n", alignment_length))
cat(sprintf("Protein recovered        : %.1f %%\n\n", 100 * qry_len / ref_len))

cat(sprintf("Protein identity (aligned region): %g %%\n", identity_val))
cat(sprintf("Identical residues               : %d/%d\n", eval_len - length(diff_pos_relative), eval_len))
cat(sprintf("Differences                      : %d amino acid substitutions\n\n", length(diff_pos_relative)))

if (nrow(mutation_table) > 0) {
  print(head(mutation_table, 15))
  if (nrow(mutation_table) > 15) {
    cat(sprintf("\n... [%d total differences written to CSV]\n", nrow(mutation_table)))
  }
  cat("\n")
} else {
  cat("No amino acid differences detected in evaluated region.\n\n")
}

###############################################################
## Output confirmation
###############################################################

cat("Mutation table written to:\n")
cat(out_csv, "\n")

cat("=================================================================\n")
cat("NOTE: ATPF1A contains introns in the genome assembly.\n")
cat("Direct translation of the genomic fragment incorporates intron\n")
cat("sequence and frame shifts. Residue substitutions are reported for\n")
cat("completeness only and are not biologically interpretable.\n")
cat("=================================================================\n\n")