source("scripts/01_project_setup.R")

project_header("Extract ATPase alpha CDS")

###############################################################
## Path definitions
###############################################################

report_file   <- file.path(PATHS$reports, "ATPF1A_tblastn_results.csv")
cds_out_file  <- file.path(PATHS$extracted, "ATPF1A_CDS.fasta")
prot_out_file <- file.path(PATHS$extracted, "ATPF1A_protein.fasta")
ref_file      <- PATHS$atpf1a

###############################################################
## Ensure output directory exists
###############################################################

dir.create(PATHS$extracted,
           recursive = TRUE,
           showWarnings = FALSE)

###############################################################
## Load BLAST results
###############################################################

cat("Loading tblastn results...\n")

hits <- read.csv(report_file)

if(nrow(hits)==0)
    stop("No BLAST hits found.")

best_hit <- hits[which.max(hits$bitscore),]

contig      <- best_hit$sseqid
start_coord <- best_hit$sstart
end_coord   <- best_hit$send

cat(sprintf(
    "Best hit on contig %s: %d to %d (Bitscore: %.1f)\n\n",
    contig,
    start_coord,
    end_coord,
    best_hit$bitscore
))

###############################################################
## Load assembly
###############################################################

assembly <- read_dna_fasta(PATHS$assembly)

###############################################################
## Extract region
###############################################################

cat("Extracting coding sequence...\n")

cds_dna <- extract_region(
    assembly = assembly,
    contig = contig,
    start = start_coord,
    end = end_coord,
    flank = 0
)

strand_char <- ifelse(start_coord < end_coord, "+", "-")

names(cds_dna) <-
    sprintf(
        "ATPF1A_CDS | %s:%d-%d | strand:%s",
        contig,
        start_coord,
        end_coord,
        strand_char
    )

cds_len <- width(cds_dna)

if(cds_len %% 3 != 0)
    stop("CDS length is not divisible by 3.")

###############################################################
## Start codon check
###############################################################

start_codon <- substr(
    as.character(cds_dna[[1]]),
    1,
    3
)

if(cds_len == 1752){

    if(start_codon != "ATG")
        warning("Complete CDS does not begin with ATG.")

} else {

    message("Partial CDS recovered; start codon cannot be evaluated.")

}

###############################################################
## Translation
###############################################################

cat("Translating CDS...\n")

protein_seq <- translate_cds(cds_dna)

names(protein_seq) <-
    sprintf(
        "ATPF1A_protein | %s:%d-%d",
        contig,
        start_coord,
        end_coord
    )

print_sequence_info(protein_seq)

prot_len <- width(protein_seq)

###############################################################
## Compare to reference
###############################################################

identity_val <- NA

if(file.exists(ref_file)){

    cat("\nComparing against reference protein...\n")

    reference <- read_protein_fasta(ref_file)

    if(width(reference)==prot_len){

        identity_val <-
            round(
                percent_identity(
                    protein_seq,
                    reference
                ),
                2
            )

        cat("Protein identity :", identity_val,"%\n\n")

    }else{

        cat("Reference comparison skipped (different lengths).\n\n")

    }

}

###############################################################
## Automatic validation
###############################################################

project_header("Automatic validation")

prot_str <- as.character(protein_seq[[1]])

prot_chars <- strsplit(prot_str,"")[[1]]

internal_stop_codons <- sum(prot_chars=="*")

cat("Expected CDS length   : 1752 bp\n")
cat("Observed CDS length   :", cds_len,"bp\n")

cat("Expected protein      : 584 aa\n")
cat("Observed protein      :", prot_len,"aa\n")
cat(
    "Protein recovered    :",
    round(100 * prot_len / 584, 1),
    "%\n"
)

if(cds_len == 1752){

    cat("Start codon           :", start_codon, "\n")

}else{

    cat("Start codon           : Not evaluated (partial CDS)\n")

}

cat("Internal stop codons  :", internal_stop_codons,"\n")

cat("Frame valid           :", cds_len %% 3 == 0,"\n\n")

###############################################################
## Expected completeness
###############################################################

if(cds_len < 1752){

    warning(
        "Only a partial ATPase alpha CDS was recovered.\n",
        "This is expected from the current assembly validation ",
        "because tblastn aligned only part of the protein."
    )

}

###############################################################
## Save outputs
###############################################################

save_fasta(cds_dna, cds_out_file)
save_fasta(protein_seq, prot_out_file)

###############################################################
## Final report
###############################################################

cat("ATPase alpha CDS extracted.\n\n")

cat("Results written to\n")

cat(cds_out_file,"\n")
cat(prot_out_file,"\n")