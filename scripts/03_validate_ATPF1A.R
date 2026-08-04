source("scripts/01_project_setup.R")

validate_gene(
  protein_fasta = PATHS$atpf1a,
  gene_name = "ATPase alpha",
  database = file.path(PATHS$blast, "IVMT1"),
  output_prefix = "ATPF1A"
)