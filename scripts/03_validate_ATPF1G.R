source("scripts/functions.R")

project_header("ATPase gamma validation")

assembly <- readDNAStringSet(
  "data/assembly/Tequiperdum_IVMt1.fasta"
)

gamma <- read_protein_fasta(
  "data/reference/ATPF1G.fasta"
)

cat("Assembly contigs :", length(assembly), "\n")
cat("Reference length :", width(gamma), "aa\n\n")

blast_output <- "results/blast/ATPF1G_tblastn.tsv"

run_tblastn(
  query = "data/reference/ATPF1G.fasta",
  db = "results/blast/IVMT1",
  output = blast_output
)

hits <- read_blast_table(blast_output)

cat("Number of hits :", nrow(hits), "\n\n")

print(hits)

save_blast_results(
  hits,
  "reports/ATPF1G_tblastn_results.csv"
)

cat("\nResults saved to:\n")
cat("reports/ATPF1G_tblastn_results.csv\n")