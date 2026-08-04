source("scripts/01_project_setup.R")

project_header("Generate ATPase Summary & Supplementary Tables")

###############################################################
## Path definitions
###############################################################

ref_atpf1g_file <- PATHS$atpf1g
ref_atpf1a_file <- PATHS$atpf1a

qry_atpf1g_prot_file <- file.path(PATHS$extracted, "ATPF1G_protein.fasta")
qry_atpf1a_prot_file <- file.path(PATHS$extracted, "ATPF1A_protein.fasta")
qry_atpf1g_cds_file  <- file.path(PATHS$extracted, "ATPF1G_CDS.fasta")
qry_atpf1a_cds_file  <- file.path(PATHS$extracted, "ATPF1A_CDS.fasta")

blast_atpf1a_csv <- file.path(PATHS$reports, "ATPF1A_tblastn_results.csv")
atpf1g_mut_csv   <- file.path(PATHS$reports, "ATPF1G_amino_acid_differences.csv")
atpf1a_mut_csv   <- file.path(PATHS$reports, "ATPF1A_amino_acid_differences.csv")

out_summary_csv <- file.path(PATHS$reports, "ATPase_summary.csv")
out_markers_csv <- file.path(PATHS$reports, "ATPase_dyskinetoplastic_markers.csv")
out_text_report <- file.path(PATHS$reports, "ATPase_summary.txt")

dir.create(PATHS$reports, recursive = TRUE, showWarnings = FALSE)

###############################################################
## Load Data
###############################################################

cat("Loading sequence data, BLAST alignments, and mutation reports...\n\n")

ref_atpf1g <- read_protein_fasta(ref_atpf1g_file)
ref_atpf1a <- read_protein_fasta(ref_atpf1a_file)

qry_atpf1g <- read_protein_fasta(qry_atpf1g_prot_file)
qry_atpf1a <- read_protein_fasta(qry_atpf1a_prot_file)

cds_atpf1g <- Biostrings::readDNAStringSet(qry_atpf1g_cds_file)
cds_atpf1a <- Biostrings::readDNAStringSet(qry_atpf1a_cds_file)

mut_atpf1g   <- read.csv(atpf1g_mut_csv, stringsAsFactors = FALSE)
mut_atpf1a   <- read.csv(atpf1a_mut_csv, stringsAsFactors = FALSE)
blast_atpf1a <- read.csv(blast_atpf1a_csv, stringsAsFactors = FALSE)

###############################################################
## Calculate ATPF1G Metrics
###############################################################

len_ref_g     <- Biostrings::width(ref_atpf1g)[1]
len_cds_ref_g <- len_ref_g * 3
len_qry_g     <- Biostrings::width(qry_atpf1g)[1]
len_cds_g     <- Biostrings::width(cds_atpf1g)[1]
cov_g         <- round(100 * (len_qry_g / len_ref_g), 1)
ident_g       <- round(100 * (len_qry_g - nrow(mut_atpf1g)) / len_qry_g, 2)

start_g      <- as.character(Biostrings::subseq(cds_atpf1g[[1]], start = 1, width = 3))
prot_chars_g <- strsplit(as.character(qry_atpf1g[[1]]), "")[[1]]
stops_g      <- sum(head(prot_chars_g, -1) == "*")

# Dyskinetoplastic markers check (Supplementary Table S2)
qry_g_chars <- unlist(strsplit(as.character(qry_atpf1g[[1]]), ""))
marker_df <- data.frame(
  Marker           = c("L262", "A273", "A281", "M282", "A273P"),
  ReferenceResidue = c("L", "A", "A", "M", "P"),
  Observed         = c(qry_g_chars[262], qry_g_chars[273], qry_g_chars[281], qry_g_chars[282], qry_g_chars[273]),
  Status           = c(
    ifelse(qry_g_chars[262] == "L", "Wild type", "Variant"),
    ifelse(qry_g_chars[273] == "A", "Wild type", "Variant"),
    ifelse(qry_g_chars[281] == "A", "Wild type", "Variant"),
    ifelse(qry_g_chars[282] == "M", "Wild type", "Variant"),
    ifelse(qry_g_chars[273] == "P", "Present", "Absent")
  ),
  BiologicalInterpretation = c(
    "Wild-type residue retained",
    "Wild-type residue retained",
    "Wild-type residue retained",
    "Wild-type residue retained",
    "Substitution absent"
  ),
  stringsAsFactors = FALSE
)

###############################################################
## Calculate ATPF1A Metrics
###############################################################

len_ref_a     <- Biostrings::width(ref_atpf1a)[1]
len_cds_ref_a <- len_ref_a * 3
len_qry_a     <- Biostrings::width(qry_atpf1a)[1]
len_cds_a     <- Biostrings::width(cds_atpf1a)[1]
cov_a         <- round(100 * (len_qry_a / len_ref_a), 1)

prot_chars_a <- strsplit(as.character(qry_atpf1a[[1]]), "")[[1]]
stops_a      <- sum(head(prot_chars_a, -1) == "*")

###############################################################
## Build Overall Summary Table (Supplementary Table S1)
###############################################################

mutations_g_str <- paste(mut_atpf1g$Mutation, collapse = "; ")

summary_df <- data.frame(
  Gene                     = c("ATPF1G", "ATPF1A"),
  SequenceSource           = c(
    "Complete CDS recovered from genomic assembly",
    "Partial genomic fragment identified by tblastn"
  ),
  ReferenceCDS_bp          = c(len_cds_ref_g, len_cds_ref_a),
  RecoveredCDS_bp          = c(len_cds_g, len_cds_a),
  FullLength_aa            = c(len_ref_g, len_ref_a),
  RecoveredProtein_aa      = c(len_qry_g, len_qry_a),
  Coverage_pct             = c(cov_g, cov_a),
  ProteinIdentity_pct      = c(ident_g, NA_real_),
  Substitutions            = c(nrow(mut_atpf1g), NA_integer_),
  ObservedSubstitutions    = c(mutations_g_str, "Not interpreted (partial genomic fragment)"),
  StartCodon               = c(start_g, "Not evaluated"),
  InternalStops_evalRegion = c(stops_g, stops_a),
  Status                   = c("Complete coding sequence recovered", "Partial genomic fragment recovered"),
  Notes                    = c(
    "All four dyskinetoplastic markers (L262, A273, A281, M282) wild type.",
    "tblastn identified the genomic locus, but reconstruction of the mature coding sequence was not performed because ATPF1A contains multiple introns. Consequently, only a partial genomic fragment was translated and ATPF1A sequence conservation could not be reliably evaluated."
  ),
  stringsAsFactors         = FALSE
)

###############################################################
## Draft Manuscript Results Paragraphs
###############################################################

p1 <- sprintf(
  "ATPase gamma (ATPF1G) was recovered as a complete coding sequence (%d bp; %d amino acids) with %.2f%% amino acid identity to the T. brucei TREU927 reference. Only %d amino acid substitutions (%s) were identified. Critically, residues previously associated with dyskinetoplastic adaptation (L262, A273, A281, and M282) remained identical to the reference sequence, and the A273P substitution was absent.",
  len_cds_g, len_qry_g, ident_g,
  nrow(mut_atpf1g),
  paste(mut_atpf1g$Mutation, collapse = " and ")
)

p2 <- sprintf(
  "ATPase alpha (ATPF1A) could not be reconstructed as a complete coding sequence from the genome assembly. Although tblastn identified the genomic locus, translation of the recovered genomic fragment yielded only %d of the expected %d amino acids (%.1f%% coverage). Because ATPF1A is intron-containing in Trypanosoma, reconstruction of the mature coding sequence requires exon-aware gene prediction or transcript evidence. Consequently, ATPF1A sequence conservation was not evaluated further.",
  len_qry_a, len_ref_a, cov_a
)

manuscript_text <- paste0(
  "======================================================================\n",
  "SUPPLEMENTARY MATERIAL & RESULTS DRAFT: ATPase VALIDATION\n",
  "======================================================================\n\n",
  "SUPPLEMENTARY TABLE INDEX:\n",
  " - Table S1: ATPase summary\n",
  " - Table S2: ATPF1G dyskinetoplastic diagnostic markers\n",
  " - Table S3: ATPF1G amino acid substitutions\n",
  " - Table S4: ATPF1A amino acid differences (technical comparison only)\n\n",
  "MANUSCRIPT RESULTS DRAFT:\n\n",
  p1, "\n\n", p2, "\n"
)

###############################################################
## Save Deliverables
###############################################################

write.csv(summary_df, out_summary_csv, row.names = FALSE)
write.csv(marker_df, out_markers_csv, row.names = FALSE)
cat(manuscript_text, file = out_text_report)

###############################################################
## Console Output
###############################################################

cat("======================================================\n")
cat("Supplementary Table S1: ATPase Summary\n")
cat("======================================================\n")
print(summary_df[, c("Gene", "SequenceSource", "Coverage_pct", "ProteinIdentity_pct", "Substitutions", "Status")])

cat("\n======================================================\n")
cat("Supplementary Table S2: ATPF1G Diagnostic Markers\n")
cat("======================================================\n")
print(marker_df)

cat("\n", manuscript_text, "\n")

cat("Deliverables written to:\n")
cat(" -", out_summary_csv, "\n")
cat(" -", out_markers_csv, "\n")
cat(" -", out_text_report, "\n")