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
    output,
    max_target_seqs = 5,
    evalue = "1e-10"
){

    system2(
        "tblastn",
        args = c(
            "-query", query,
            "-db", database,
            "-outfmt", "6",
            "-evalue", evalue,
            "-max_target_seqs", max_target_seqs,
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
## Extract genomic region from top tblastn hit
###############################################################

extract_tblastn_hit <- function(
    assembly_fasta,
    blast_row,
    outfile,
    flank = 0
){

    genome <- Biostrings::readDNAStringSet(assembly_fasta)

    seqname <- blast_row$sseqid

    idx <- grep(
        paste0("^", seqname),
        names(genome)
    )

    if(length(idx) == 0){
        stop("Contig ", seqname, " not found in assembly.")
    }

    idx <- idx[1]

    start <- min(blast_row$sstart, blast_row$send)
    end   <- max(blast_row$sstart, blast_row$send)

    seq_length <- length(genome[[idx]])

    start <- max(1, start - flank)
    end   <- min(seq_length, end + flank)

    region <- Biostrings::subseq(
        genome[[idx]],
        start = start,
        end   = end
    )

    region <- Biostrings::DNAStringSet(region)

    names(region) <- paste0(
        blast_row$qseqid,
        "_",
        seqname,
        "_",
        start,
        "_",
        end
    )

    Biostrings::writeXStringSet(region, outfile)

    invisible(region)
}

###############################################################
## Translate a single reading frame
###############################################################

translate_reading_frame <- function(dna, frame){

  dna <- Biostrings::DNAString(dna)

  seq_length <- length(dna)

  remainder <- (seq_length - frame + 1) %% 3

  if(remainder > 0){
    dna <- Biostrings::subseq(
      dna,
      start = 1,
      end   = seq_length - remainder
    )
  }

  dna <- Biostrings::subseq(dna, start = frame)

  Biostrings::translate(
    dna,
    if.fuzzy.codon = "X"
  )

}

###############################################################
## Extract candidate ORFs from one translated frame
###############################################################

extract_candidate_orfs <- function(
    translated_aa,
    min_length = 30
){

  residues <- as.character(translated_aa)

  segments <- strsplit(residues, "\\*", fixed = FALSE)[[1]]

  segments <- segments[nchar(segments) >= min_length]

  if(length(segments) == 0){
    return(character(0))
  }

  segments

}

###############################################################
## Local alignment metrics (base R)
###############################################################

local_alignment_metrics <- function(
    pattern,
    subject,
    match_score = 1,
    mismatch_score = -1,
    gap_score = -1
){

  pattern <- strsplit(as.character(pattern), "")[[1]]
  subject <- strsplit(as.character(subject), "")[[1]]

  n <- length(pattern)
  m <- length(subject)

  if(n == 0 || m == 0){
    return(list(
      identity   = 0,
      coverage   = 0,
      aligned_ref = 0
    ))
  }

  score_matrix <- matrix(
    0,
    nrow = n + 1,
    ncol = m + 1
  )

  best_score <- 0
  best_i <- 0
  best_j <- 0

  for(i in seq_len(n)){
    for(j in seq_len(m)){

      diagonal <- score_matrix[i, j] +
        if(pattern[i] == subject[j]) match_score else mismatch_score

      current <- max(
        0,
        diagonal,
        score_matrix[i, j + 1] + gap_score,
        score_matrix[i + 1, j] + gap_score
      )

      score_matrix[i + 1, j + 1] <- current

      if(current > best_score){
        best_score <- current
        best_i <- i
        best_j <- j
      }

    }
  }

  if(best_score <= 0){
    return(list(
      identity    = 0,
      coverage    = 0,
      aligned_ref = 0
    ))
  }

  aligned_ref <- 0
  matches <- 0
  i <- best_i
  j <- best_j

  while(i > 0 && j > 0 && score_matrix[i + 1, j + 1] > 0){

    aligned_ref <- aligned_ref + 1

    if(pattern[i] == subject[j]){
      matches <- matches + 1
    }

    diagonal <- score_matrix[i, j]
    up <- score_matrix[i, j + 1]
    left <- score_matrix[i + 1, j]
    current <- score_matrix[i + 1, j + 1]

    if(current == diagonal + if(pattern[i] == subject[j])
      match_score else mismatch_score){
      i <- i - 1
      j <- j - 1
    } else if(current == up + gap_score){
      j <- j - 1
    } else {
      i <- i - 1
    }

  }

  list(
    identity    = if(aligned_ref > 0) 100 * matches / aligned_ref else 0,
    coverage    = 100 * aligned_ref / n,
    aligned_ref = aligned_ref
  )

}

###############################################################
## Score one ORF against a reference protein
##
## Composite ORF score (higher is better):
##
##   score = (identity / 100) * (coverage / 100) * length_factor
##           - (0.10 * internal_stops)
##
## where:
##   identity      = local pairwise alignment identity to reference (%)
##   coverage      = aligned reference residues / reference length (%)
##   length_factor = min(orf_length / reference_length, 1)
##   internal_stops = number of in-frame stop codons within the ORF
##
## Tie-breaking (applied by find_best_orf):
##   higher identity, then higher coverage, then longer ORF.
###############################################################

score_orf_against_reference <- function(
    orf,
    reference
){

  orf <- Biostrings::AAString(orf)

  reference <- Biostrings::AAString(reference)

  orf_length <- length(orf)

  ref_length <- length(reference)

  if(orf_length == 0 || ref_length == 0){
    return(NULL)
  }

  orf_chars <- strsplit(as.character(orf), "")[[1]]

  internal_stops <- sum(orf_chars == "*")

  alignment <- local_alignment_metrics(
    pattern = reference,
    subject = orf
  )

  identity <- alignment$identity
  coverage <- alignment$coverage

  length_factor <- min(orf_length / ref_length, 1)

  score <- (identity / 100) *
    (coverage / 100) *
    length_factor -
    (0.10 * internal_stops)

  list(
    protein        = orf,
    identity       = identity,
    coverage       = coverage,
    orf_length     = as.integer(orf_length),
    internal_stops = as.integer(internal_stops),
    score          = score
  )

}

###############################################################
## Evaluate translated ORFs and return the best match
###############################################################

find_best_orf <- function(
    genomic_sequence,
    reference_protein,
    min_orf_length = 30
){

  genomic_sequence <- Biostrings::DNAString(genomic_sequence)

  reference_protein <- Biostrings::AAString(reference_protein)

  ref_length <- length(reference_protein)

  min_length <- max(
    min_orf_length,
    floor(ref_length * 0.10)
  )

  candidates <- list()

  orientations <- c("plus", "minus")

  frames <- 1:3

  for(orientation in orientations){

    if(orientation == "plus"){
      template <- genomic_sequence
    } else {
      template <- Biostrings::reverseComplement(genomic_sequence)
    }

    for(frame in frames){

      translated <- translate_reading_frame(
        template,
        frame
      )

      orfs <- extract_candidate_orfs(
        translated,
        min_length = min_length
      )

      if(length(orfs) == 0){
        next
      }

      for(orf in orfs){

        scored <- score_orf_against_reference(
          orf,
          reference_protein
        )

        if(is.null(scored)){
          next
        }

        scored$frame <- frame
        scored$orientation <- orientation

        candidates[[length(candidates) + 1]] <- scored

      }

    }

  }

  if(length(candidates) == 0){
    stop("No candidate ORFs met the minimum length threshold.")
  }

  candidate_table <- do.call(
    rbind,
    lapply(candidates, function(x){
      data.frame(
        frame          = x$frame,
        orientation    = x$orientation,
        identity       = x$identity,
        coverage       = x$coverage,
        orf_length     = x$orf_length,
        internal_stops = x$internal_stops,
        score          = x$score,
        stringsAsFactors = FALSE
      )
    })
  )

  best_idx <- order(
    -candidate_table$score,
    -candidate_table$identity,
    -candidate_table$coverage,
    -candidate_table$orf_length
  )[1]

  best <- candidates[[best_idx]]

  list(
    protein        = best$protein,
    frame          = best$frame,
    orientation    = best$orientation,
    identity       = best$identity,
    coverage       = best$coverage,
    orf_length     = best$orf_length,
    internal_stops = best$internal_stops,
    score          = best$score
  )

}

###############################################################
## Recover a gene from the best tblastn hit using six-frame ORFs
###############################################################

recover_gene_from_tblastn <- function(
    reference_protein,
    assembly,
    blast_database,
    gene_name,
    output_directory,
    flank = 0,
    max_target_seqs = 5,
    evalue = "1e-10"
){

  if(!file.exists(reference_protein)){
    stop("Reference protein not found: ", reference_protein)
  }

  if(!file.exists(assembly)){
    stop("Assembly not found: ", assembly)
  }

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  tblastn_file <- file.path(
    output_directory,
    paste0(gene_name, "_tblastn.tsv")
  )

  genomic_file <- file.path(
    output_directory,
    paste0(gene_name, "_genomic.fasta")
  )

  protein_file <- file.path(
    output_directory,
    paste0(gene_name, "_protein.fasta")
  )

  ref_protein <- read_protein_fasta(reference_protein)

  run_tblastn(
    query           = reference_protein,
    database        = blast_database,
    output          = tblastn_file,
    max_target_seqs = max_target_seqs,
    evalue          = evalue
  )

  hits <- read_blast_table(tblastn_file)

  if(nrow(hits) == 0){
    stop("No tblastn hits recovered for ", gene_name, ".")
  }

  best_hit <- hits[
    order(-hits$bitscore),
  ][1, , drop = FALSE]

  genomic <- extract_tblastn_hit(
    assembly_fasta = assembly,
    blast_row      = best_hit,
    outfile        = genomic_file,
    flank          = flank
  )

  minus_strand <- best_hit$sstart > best_hit$send

  if(minus_strand){
    genomic <- Biostrings::reverseComplement(genomic)
    names(genomic) <- paste0(names(genomic), "_oriented")
  }

  best_orf <- find_best_orf(
    genomic_sequence  = genomic[[1]],
    reference_protein = ref_protein[[1]]
  )

  save_fasta(genomic, genomic_file)

  recovered_protein <- Biostrings::AAStringSet(best_orf$protein)

  names(recovered_protein) <- paste0(gene_name, "_recovered")

  save_fasta(recovered_protein, protein_file)

  statistics <- list(
    frame          = best_orf$frame,
    orientation    = best_orf$orientation,
    identity       = best_orf$identity,
    coverage       = best_orf$coverage,
    orf_length     = best_orf$orf_length,
    internal_stops = best_orf$internal_stops,
    score          = best_orf$score,
    blast_contig   = best_hit$sseqid,
    blast_strand   = ifelse(minus_strand, "minus", "plus"),
    blast_pident   = best_hit$pident,
    blast_evalue   = best_hit$evalue,
    blast_bitscore = best_hit$bitscore
  )

  list(
    protein    = recovered_protein,
    genomic    = genomic,
    tblastn    = hits,
    statistics = statistics
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

###############################################################
## Summarize recovered gene
###############################################################

summarize_gene <- function(
    gene,
    ref_protein,
    qry_protein,
    qry_cds = NULL,
    mutations = NULL,
    sequence_source = NA_character_,
    notes = "",
    status = "Complete CDS recovered",
    evaluate_start_codon = TRUE,
    coverage_pct = NULL
){

    ref_len <- Biostrings::width(ref_protein)[1]
    qry_len <- Biostrings::width(qry_protein)[1]

    ref_cds_bp <- ref_len * 3

    if(is.null(qry_cds)){
    cds_len <- NA_integer_
    start_codon <- NA_character_

    } else {

        cds_len <- Biostrings::width(qry_cds)[1]

        if(evaluate_start_codon){

            start_codon <- as.character(
                Biostrings::subseq(qry_cds[[1]], start = 1, width = 3)
            )

        } else {

            start_codon <- NA_character_

        }

    }

    ## Coverage_pct definition
    ##
    ## When `coverage_pct` is supplied (e.g. by recover_gene_from_tblastn(),
    ## which computes true local-alignment coverage in
    ## local_alignment_metrics()/score_orf_against_reference()/find_best_orf()
    ## as 100 * aligned reference residues / reference length), that value is
    ## used directly so Coverage_pct always reflects genuine sequence
    ## homology rather than raw sequence length.
    ##
    ## When no alignment-based coverage is available (curated CDS extraction
    ## scripts that call summarize_gene() directly without an ORF search),
    ## coverage falls back to the length ratio recovered/reference protein
    ## length. This fallback is only a valid proxy for true coverage when the
    ## recovered sequence is a clean, in-frame CDS translation; it will
    ## overstate coverage if the extracted region contains undetected
    ## frameshifts (e.g. uncorrected introns).
    if(is.null(coverage_pct)){

        coverage <- round(100 * qry_len / ref_len, 1)

    } else {

        coverage <- round(coverage_pct, 1)

    }

    prot_chars <- strsplit(as.character(qry_protein[[1]]), "")[[1]]
    internal_stops <- sum(head(prot_chars, -1) == "*")

    if(is.null(mutations)){

        n_mut <- NA_integer_
        ident <- NA_real_
        mut_string <- NA_character_

    } else {

        n_mut <- nrow(mutations)

        ident <- round(
            100 * (qry_len - n_mut) / qry_len,
            2
        )

        if(n_mut == 0){

            mut_string <- "None"

        } else {

            mut_string <- paste(
                mutations$Mutation,
                collapse = "; "
            )

        }

    }

    data.frame(

        Gene                     = gene,
        ReferenceCDS_bp          = ref_cds_bp,
        RecoveredCDS_bp          = cds_len,
        FullLength_aa            = ref_len,
        RecoveredProtein_aa      = qry_len,
        Coverage_pct             = coverage,
        ProteinIdentity_pct      = ident,
        Substitutions            = n_mut,
        ObservedSubstitutions    = mut_string,
        StartCodon               = start_codon,
        InternalStops_evalRegion = internal_stops,
        SequenceSource           = sequence_source,
        Status                   = status,
        Notes                    = notes,

        stringsAsFactors = FALSE
    )

}
