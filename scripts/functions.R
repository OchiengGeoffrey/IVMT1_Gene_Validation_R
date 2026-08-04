###############################################################
## Utility Functions
##
## IVMT1_Gene_Validation_R
##
## This file contains reusable helper functions for:
##   - FASTA input/output
##   - Sequence summaries
##   - BLAST automation
##   - Genomic region extraction
##
## Author: Geoffrey Omondi
###############################################################

library(Biostrings)

save_blast_results <- function(tbl, filename){

  directory <- dirname(filename)

  if(!dir.exists(directory)){
    dir.create(directory,
               recursive = TRUE)
  }

  write.csv(
    tbl,
    filename,
    row.names = FALSE
  )

}

###############################################################
## Read a FASTA sequence
###############################################################

read_dna_fasta <- function(file){

  if(!file.exists(file)){
    stop("File not found: ", file)
  }

  dna <- readDNAStringSet(file)

  if(length(dna) == 0){
    stop("No DNA sequences found in ", file)
  }

  dna

}

read_protein_fasta <- function(file){

  if(!file.exists(file)){
    stop("File not found: ", file)
  }

  aa <- readAAStringSet(file)

  if(length(aa) == 0){
    stop("No protein sequences found in ", file)
  }

  aa

}

save_fasta <- function(sequence, filename){

  directory <- dirname(filename)

  if(!dir.exists(directory)){
    dir.create(directory, recursive = TRUE)
  }

  Biostrings::writeXStringSet(
    sequence,
    filepath = filename,
    format = "fasta"
  )

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

  if(length(sequence) == 0){
    stop("Sequence object is empty.")
  }

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

  x <- strsplit(as.character(sequence)[1], "")[[1]]

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
## Calculate Identity
###############################################################
percent_identity <- function(a, b){

  a <- strsplit(as.character(a)[1], "")[[1]]
  b <- strsplit(as.character(b)[1], "")[[1]]

  if(length(a) != length(b))
    stop("Sequences must have equal length.")

  identical_positions <- sum(a == b)

  100 * identical_positions / length(a)

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

  aa

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
## Extract genomic region
###############################################################

extract_region <- function(assembly, contig, start, end, flank = 0){

  matches <- startsWith(names(assembly), contig)

  if(sum(matches) == 0)
    stop("Contig not found: ", contig)

  if(sum(matches) > 1)
    stop("Multiple contigs matched: ", contig)

  chr <- assembly[matches][[1]]

  chr_length <- nchar(as.character(chr))

  left  <- max(1, min(start, end) - flank)
  right <- min(chr_length, max(start, end) + flank)

  region <- Biostrings::subseq(
    chr,
    start = left,
    end = right
  )

  if(start > end){
    region <- Biostrings::reverseComplement(region)
  }

  region_set <- Biostrings::DNAStringSet(region)

  names(region_set) <- paste0(
    contig,
    ":",
    left,
    "-",
    right
  )

  return(region_set)

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

  status <- system2(
  "makeblastdb",
  args = c(
    "-in", fasta,
    "-dbtype", dbtype,
    "-out", db_name
  )
)

if(status != 0){
  stop("makeblastdb failed.")
}

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

  if(!file.exists(file)){
    stop("BLAST output not found: ", file)
  }

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
    col.names = cols,
    stringsAsFactors = FALSE
  )

}

###############################################################
## Run tblastn
###############################################################

run_tblastn <- function(query, db, output){

  if(!file.exists(query)){
    stop("Query FASTA not found: ", query)
  }

  cmd <- paste(
    "tblastn",
    "-query", shQuote(query),
    "-db", shQuote(db),
    "-out", shQuote(output),
    "-outfmt 6",
    "-max_target_seqs 5",
    "-evalue 1e-5"
  )

  cat(cmd, "\n\n")

  status <- system(cmd)

  if(status != 0){
    stop("tblastn failed.")
  }

  if(!file.exists(output)){
    stop("tblastn completed but no output file was produced.")
  }

}

###############################################################
## Save BLAST results
###############################################################

save_blast_results <- function(tbl, filename){

  directory <- dirname(filename)

  if(!dir.exists(directory)){
    dir.create(directory,
               recursive = TRUE)
  }

  write.csv(
    tbl,
    filename,
    row.names = FALSE
  )

}

###############################################################
## Gene Validation
###############################################################
validate_gene <- function(
    protein_fasta,
    gene_name,
    database,
    output_prefix
){

  ###############################################################
  ## Header
  ###############################################################

  project_header(paste(gene_name, "validation"))

  ###############################################################
  ## Input checks
  ###############################################################

  if (!file.exists(protein_fasta)) {
    stop("Protein FASTA not found: ", protein_fasta)
  }

  ###############################################################
  ## Load protein
  ###############################################################

  protein <- read_protein_fasta(protein_fasta)

  ###############################################################
  ## Assembly summary
  ###############################################################

  assembly <- read_dna_fasta(PATHS$assembly)

  cat("Assembly contigs :", length(assembly), "\n")
  cat("Reference length :", width(protein), "aa\n\n")

  ###############################################################
  ## Run tblastn
  ###############################################################

  blast_output <- file.path(
    PATHS$blast,
    paste0(output_prefix, "_tblastn.tsv")
  )

  run_tblastn(
    query = protein_fasta,
    db = database,
    output = blast_output
  )

  ###############################################################
  ## Read BLAST results
  ###############################################################

  hits <- read_blast_table(blast_output)

  cat("Number of hits :", nrow(hits), "\n\n")

  ###############################################################
  ## Best hit
  ###############################################################

  best_hit <- hits |>
    dplyr::arrange(dplyr::desc(bitscore)) |>
    dplyr::slice(1)

  cat("Best hit\n")
  cat("-------------------------\n")
  cat("Contig      :", best_hit$sseqid, "\n")
  cat("Identity    :", round(best_hit$pident, 2), "%\n")
  cat("Alignment   :", best_hit$length, "aa\n")
  cat("E-value     :", best_hit$evalue, "\n")
  cat("Bit score   :", best_hit$bitscore, "\n")
  cat(
    "Coordinates :",
    best_hit$sstart,
    "-",
    best_hit$send,
    "\n\n"
  )

  ###############################################################
  ## Save results
  ###############################################################

  report_file <- file.path(
    PATHS$reports,
    paste0(output_prefix, "_tblastn_results.csv")
  )

  save_blast_results(
    hits,
    report_file
  )

  cat("Results saved to:\n")
  cat(report_file, "\n")

  invisible(best_hit)
}