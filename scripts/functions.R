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

  dna <- readDNAStringSet(file)

  return(dna)

}

read_protein_fasta <- function(file){

  if(!file.exists(file))
    stop("File not found: ", file)

  aa <- readAAStringSet(file)

  return(aa)

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

  cat(substr(as.character(sequence)[1], 1, 60),"...\n")

  cat("---------------------------------\n")

}

###############################################################
## Amino-acid composition
###############################################################

aa_composition <- function(sequence){

  x <- strsplit(as.character(sequence), "")[[1]]

  sort(table(x), decreasing = TRUE)

}

###############################################################
## Print a nice header
###############################################################
project_header <- function(title){

  cat("\n")
  cat(rep("=",70), sep="")
  cat("\n")
  cat(title)
  cat("\n")
  cat(rep("=",70), sep="")
  cat("\n\n")

}

###############################################################
## save FASTA
###############################################################
save_fasta <- function(sequence, filename){

  Biostrings::writeXStringSet(
    sequence,
    filepath = filename,
    format = "fasta"
  )

}

###############################################################
## Calculate Identity
###############################################################
percent_identity <- function(a, b){

  a <- strsplit(as.character(a)[1], "")[[1]]
  b <- strsplit(as.character(b)[1], "")[[1]]

  if(length(a) != length(b))
    stop("Sequences must have equal length.")

  identical_positions <- sum(a == b)

  return(100 * identical_positions / length(a))

}

###############################################################
## Reverse Complement
###############################################################
reverse_complement <- function(seq){

  Biostrings::reverseComplement(seq)

}

###############################################################
## Translate CDS
###############################################################
translate_cds <- function(cds){

  aa <- Biostrings::translate(cds)

  return(aa)

}

###############################################################
## Helper
###############################################################
save_single_fasta <- function(sequence,
                              filename,
                              header){

  writeLines(
    c(
      paste0(">", header),
      as.character(sequence)
    ),
    filename
  )

}

assembly_statistics <- function(assembly){

  tibble::tibble(

    Contigs = length(assembly),

    Total_bp = sum(width(assembly)),

    Longest_contig = max(width(assembly))

  )

}

###############################################################
## Create BLAST database
###############################################################

create_blast_database <- function(
    fasta,
    db_name,
    dbtype = "nucl"
){

  if(!file.exists(fasta))
    stop("Assembly not found.")

  marker <- paste0(db_name, ".nin")

  if(file.exists(marker)){

    message("BLAST database already exists.")

    return(invisible(TRUE))

  }

  dir.create(dirname(db_name),
             recursive = TRUE,
             showWarnings = FALSE)

  system2(
    "makeblastdb",
    args = c(
      "-in", fasta,
      "-dbtype", dbtype,
      "-out", db_name
    )
  )

}

###############################################################
## Run tblastn
###############################################################

run_tblastn <- function(
    query,
    database,
    output
){

  system2(
    "tblastn",
    args = c(
      "-query", query,
      "-db", database,
      "-outfmt", "6",
      "-evalue", "1e-10",
      "-max_target_seqs", "5",
      "-out", output
    )
  )

}

###############################################################
## Read BLAST table
###############################################################

read_blast_table <- function(file){

  cols <- c(

    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore"

  )

  read.delim(
    file,
    header = FALSE,
    sep = "\t",
    col.names = cols
  )

}