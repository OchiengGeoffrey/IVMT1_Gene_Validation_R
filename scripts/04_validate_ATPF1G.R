source("scripts/01_project_setup.R")

validate_gene(
  protein_fasta = PATHS$atpf1g,
  gene_name = "ATPase gamma",
  database = file.path(PATHS$blast, "IVMT1"),
  output_prefix = "ATPF1G"
)