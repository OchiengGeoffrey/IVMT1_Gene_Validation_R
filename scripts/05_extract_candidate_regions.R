source("scripts/01_project_setup.R")

project_header("Extract candidate genomic regions")

assembly <- readDNAStringSet(PATHS$assembly)

dir.create(
  PATHS$extracted,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################
## ATPase gamma
###############################################################

gamma_hits <- read.csv(
  file.path(PATHS$reports, "ATPF1G_tblastn_results.csv")
)

best_gamma <- gamma_hits[which.max(gamma_hits$bitscore), ]

gamma_region <- extract_region(
  assembly = assembly,
  contig = best_gamma$sseqid,
  start = best_gamma$sstart,
  end = best_gamma$send,
  flank = 1000
)

save_fasta(
  gamma_region,
  file.path(PATHS$extracted, "ATPF1G_region.fasta")
)

cat("ATPase gamma region extracted.\n")
cat(names(gamma_region), "\n")
cat("Length:", width(gamma_region), "bp\n")