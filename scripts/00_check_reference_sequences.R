source("scripts/01_project_setup.R")

project_header("Reference proteins")

alpha <- read_protein_fasta(PATHS$atpf1a)
gamma <- read_protein_fasta(PATHS$atpf1g)

print_sequence_info(alpha)
print_sequence_info(gamma)