source("scripts/01_project_setup.R")

project_header("Compare ATPase gamma proteins")

###############################################################
## Path definitions
###############################################################

ref_file <- PATHS$atpf1g
qry_file <- file.path(PATHS$extracted, "ATPF1G_protein.fasta")
out_csv  <- file.path(PATHS$reports, "ATPF1G_amino_acid_differences.csv")

###############################################################
## Load protein sequences
###############################################################

cat("Loading protein sequences...\n")

if (!file.exists(ref_file)) {
  stop("Reference protein file not found: ", ref_file)
}

if (!file.exists(qry_file)) {
  stop("Query protein file not found: ", qry_file)
}

reference <- read_protein_fasta(ref_file)
ivmt1     <- read_protein_fasta(qry_file)

ref_len <- Biostrings::width(reference)[1]
qry_len <- Biostrings::width(ivmt1)[1]

if (ref_len != qry_len) {
  stop(sprintf("Protein lengths differ! Reference: %d aa, IVM-t1: %d aa", ref_len, qry_len))
}

###############################################################
## Sequence comparison & mutation identification
###############################################################

ref_chars <- unlist(strsplit(as.character(reference[[1]]), ""))
qry_chars <- unlist(strsplit(as.character(ivmt1[[1]]), ""))

diff_pos <- which(ref_chars != qry_chars)

mutation_table <- data.frame(
  Position  = diff_pos,
  Mutation  = paste0(ref_chars[diff_pos], diff_pos, qry_chars[diff_pos]),
  Reference = ref_chars[diff_pos],
  IVM_t1    = qry_chars[diff_pos],
  stringsAsFactors = FALSE
)

identity_val <- round(percent_identity(ivmt1, reference), 2)

###############################################################
## Save mutation report
###############################################################

dir.create(PATHS$reports, recursive = TRUE, showWarnings = FALSE)
save_blast_results(mutation_table, out_csv)

###############################################################
## Console report
###############################################################

cat(sprintf("Reference length   : %d aa\n", ref_len))
cat(sprintf("IVMT1 length       : %d aa\n\n", qry_len))

cat(sprintf("Protein identity   : %g %%\n", identity_val))
cat(sprintf("Identical residues : %d/%d\n", ref_len - length(diff_pos), ref_len))
cat(sprintf("Differences        : %d amino acid substitutions\n\n", length(diff_pos)))

if (nrow(mutation_table) > 0) {
  print(mutation_table)
  cat("\n")
} else {
  cat("No amino acid differences detected (100% identity).\n\n")
}

###############################################################
## Diagnostic check for dyskinetoplastic markers
###############################################################

cat("Key ATPase gamma residues\n")
cat("-------------------------\n")

markers <- data.frame(
  Position = c(262, 273, 281, 282),
  Expected = c("L", "A", "A", "M"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(markers))) {
  pos    <- markers$Position[i]
  exp    <- markers$Expected[i]
  obs    <- qry_chars[pos]
  status <- if (obs == exp) "Wild type" else "Variant"
  
  cat(sprintf("%s%d : %s (%s)\n", exp, pos, obs, status))
}

cat(sprintf("A273P marker : %s\n\n", ifelse(qry_chars[273] == "P", "Present", "Absent")))

cat("Mutation table written to:\n")
cat(out_csv, "\n")