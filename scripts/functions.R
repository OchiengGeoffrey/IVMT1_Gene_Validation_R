###############################################################
## Utility functions
## IVMT1 Gene Validation
###############################################################

library(Biostrings)

###############################################################
## Read a FASTA sequence
###############################################################

read_dna_fasta <- function(file){

  if(!file.exists(file))
    stop("File not found: ", file)

  readDNAStringSet(file)

}

read_protein_fasta <- function(file){

  if(!file.exists(file))
    stop("File not found: ", file)

  readAAStringSet(file)

}

###############################################################
## Sequence summary
###############################################################

sequence_summary <- function(sequence){

  tibble::tibble(

    Name = names(sequence),

    Length = width(sequence)

  )

}

###############################################################
## Print sequence information
###############################################################

print_sequence_info <- function(sequence){

  cat("\n")

  cat("---------------------------------\n")

  cat("Sequence name:\n")

  cat(names(sequence), "\n\n")

  cat("Length:\n")

  cat(width(sequence), "\n\n")

  cat("First 60 residues:\n")

  cat(substr(as.character(sequence),1,60),"...\n")

  cat("---------------------------------\n")

}

###############################################################
## Amino-acid composition
###############################################################

aa_composition <- function(sequence){

  x <- strsplit(as.character(sequence), "")[[1]]

  sort(table(x), decreasing = TRUE)

}