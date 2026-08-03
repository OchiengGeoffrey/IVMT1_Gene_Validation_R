source("scripts/functions.R")

project_header("Creating BLAST database")

assembly <- "data/assembly/Tequiperdum_IVMt1.fasta"

dir.create("results/blast", recursive = TRUE, showWarnings = FALSE)

cmd <- paste(
  "makeblastdb",
  "-in", shQuote(assembly),
  "-dbtype nucl",
  "-out results/blast/IVMT1"
)

cat(cmd, "\n")

system(cmd)