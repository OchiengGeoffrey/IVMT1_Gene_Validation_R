###############################################################
## IVMT1_Gene_Validation_R
##
## Project configuration
###############################################################

library(Biostrings)

###############################################################
## Project information
###############################################################

PROJECT_NAME <- "IVMT1_Gene_Validation_R"

###############################################################
## Project paths
###############################################################

PATHS <- list(

  assembly = "data/assembly/Tequiperdum_IVMt1.fasta",

  reads1 = "data/reads/SRR7910035_1.fastq",

  reads2 = "data/reads/SRR7910035_2.fastq",

  atpf1a = "data/reference/ATPF1A.fasta",

  atpf1g = "data/reference/ATPF1G.fasta",

  blast = "results/blast",

  extracted = "results/extracted_regions",

  mapping = "results/mapping",

  consensus = "results/consensus",

  variants = "results/variants",

  reports = "reports"

)

###############################################################
## Load helper functions
###############################################################

source("scripts/functions.R")