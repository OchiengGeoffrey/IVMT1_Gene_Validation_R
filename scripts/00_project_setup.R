###############################################################
## IVMT1 Gene Validation
##
## Independent validation of ATP synthase alpha and gamma
## using R/Bioconductor.
##
## Author: Geoffrey Omondi
###############################################################

library(Biostrings)

source("scripts/functions.R")

alpha <- read_protein(alpha_reference)

gamma <- read_protein(gamma_reference)

print_sequence_info(alpha)

print_sequence_info(gamma)