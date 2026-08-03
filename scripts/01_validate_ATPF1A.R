source("scripts/functions.R")

project_header("ATPase alpha validation")

assembly_file <-
  "data/assembly/Tequiperdum_IVMt1.fasta"

reference_file <-
  "data/reference/ATPF1A.fasta"

blast_db <-
  "data/blastdb/IVMT1"

blast_output <-
  "results/blast/ATPF1A_tblastn.tsv"

assembly <- read_dna_fasta(assembly_file)

alpha <- read_protein_fasta(reference_file)

cat("Assembly contigs :", length(assembly), "\n")
cat("Reference length :", width(alpha), "aa\n")

create_blast_database(
  fasta = assembly_file,
  db_name = blast_db
)

run_tblastn(
  query = reference_file,
  database = blast_db,
  output = blast_output
)

hits <- read_blast_table(blast_output)

hits