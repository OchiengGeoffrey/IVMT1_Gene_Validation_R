source("scripts/functions.R")

project_header("ATPase alpha validation")

###############################################################
# Input files
###############################################################

assembly_file  <- "data/assembly/Tequiperdum_IVMt1.fasta"
reference_file <- "data/reference/ATPF1A.fasta"

###############################################################
# Load data
###############################################################

assembly <- read_dna_fasta(assembly_file)

reference <- read_protein_fasta(reference_file)

cat("Assembly contigs :", length(assembly), "\n")
cat("Reference length :", width(reference), "aa\n\n")

###############################################################
# Run tblastn
###############################################################

blast_output <- "results/blast/ATPF1A_tblastn.tsv"

run_tblastn(
  query = reference_file,
  database = "results/blast/IVMT1",
  output = blast_output
)

###############################################################
# Read results
###############################################################

hits <- read_blast_table(blast_output)

cat("Number of hits :", nrow(hits), "\n\n")

###############################################################
# Top hit
###############################################################

best_hit <- hits |>
  dplyr::arrange(desc(bitscore)) |>
  dplyr::slice(1)

print(best_hit)

###############################################################
# Save table
###############################################################

write.csv(
  hits,
  "reports/ATPF1A_tblastn_results.csv",
  row.names = FALSE
)

cat("\nResults saved to:\n")
cat("reports/ATPF1A_tblastn_results.csv\n")