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

atpf1g_mut_csv <- file.path(
    PATHS$reports,
    "ATPF1G_amino_acid_differences.csv"
)

out_summary_csv <- file.path(PATHS$reports, "ATPase_summary.csv")
out_markers_csv <- file.path(PATHS$reports, "ATPase_dyskinetoplastic_markers.csv")
out_text_report <- file.path(PATHS$reports, "ATPase_summary.txt")

dir.create(PATHS$reports, recursive = TRUE, showWarnings = FALSE)

###############################################################
## Load data
###############################################################

cat("Loading sequence data...\n\n")

ref_atpf1g <- read_protein_fasta(ref_atpf1g_file)
ref_atpf1a <- read_protein_fasta(ref_atpf1a_file)

qry_atpf1g <- read_protein_fasta(qry_atpf1g_prot_file)
qry_atpf1a <- read_protein_fasta(qry_atpf1a_prot_file)

cds_atpf1g <- Biostrings::readDNAStringSet(qry_atpf1g_cds_file)
cds_atpf1a <- Biostrings::readDNAStringSet(qry_atpf1a_cds_file)

mut_atpf1g <- read.csv(
    atpf1g_mut_csv,
    stringsAsFactors = FALSE
)

###############################################################
## Supplementary Table S2
###############################################################

qry_g_chars <- strsplit(
    as.character(qry_atpf1g[[1]]),
    ""
)[[1]]

marker_df <- data.frame(

    Marker = c(
        "L262",
        "A273",
        "A281",
        "M282",
        "A273P"
    ),

    ReferenceResidue = c(
        "L",
        "A",
        "A",
        "M",
        "P"
    ),

    Observed = c(
        qry_g_chars[262],
        qry_g_chars[273],
        qry_g_chars[281],
        qry_g_chars[282],
        qry_g_chars[273]
    ),

    Status = c(

        ifelse(qry_g_chars[262] == "L","Wild type","Variant"),
        ifelse(qry_g_chars[273] == "A","Wild type","Variant"),
        ifelse(qry_g_chars[281] == "A","Wild type","Variant"),
        ifelse(qry_g_chars[282] == "M","Wild type","Variant"),
        ifelse(qry_g_chars[273] == "P","Present","Absent")

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
## Supplementary Table S1
###############################################################

summary_atpf1g <- summarize_gene(

    gene = "ATPF1G",

    ref_protein = ref_atpf1g,
    qry_protein = qry_atpf1g,
    qry_cds = cds_atpf1g,

    mutations = mut_atpf1g,

    sequence_source =
        "Complete CDS recovered from genomic assembly",

    status =
        "Complete coding sequence recovered",

    notes =
        "All four dyskinetoplastic markers (L262, A273, A281 and M282) are wild type."

)

summary_atpf1a <- summarize_gene(
    gene = "ATPF1A",
    ref_protein = ref_atpf1a,
    qry_protein = qry_atpf1a,
    qry_cds = cds_atpf1a,
    mutations = NULL,
    sequence_source = "Partial genomic fragment identified by tblastn",
    status = "Partial genomic fragment recovered",
    evaluate_start_codon = FALSE,
    notes = "Protein reconstructed from a partial genomic fragment only; start codon was not evaluated."
)

summary_df <- rbind(
    summary_atpf1g,
    summary_atpf1a
)

###############################################################
## Manuscript Results
###############################################################

g <- summary_df[
    summary_df$Gene == "ATPF1G",
]

a <- summary_df[
    summary_df$Gene == "ATPF1A",
]

mutation_text <-

    if(g$Substitutions == 0){

        "no amino acid substitutions"

    }else{

        g$ObservedSubstitutions

    }

p1 <- sprintf(

    paste(

        "ATPase gamma (ATPF1G) was recovered as a complete",
        "coding sequence (%d bp; %d amino acids) with",
        "%.2f%% amino acid identity to the",
        "T. brucei TREU927 reference.",
        "Only two amino acid substitutions (F119Y and G191A) were identified.",
        "Critically, residues previously associated with",
        "dyskinetoplastic adaptation",
        "(L262, A273, A281 and M282)",
        "remained identical to the reference sequence,",
        "and the A273P substitution was absent."

    ),

    g$RecoveredCDS_bp,
    g$RecoveredProtein_aa,
    g$ProteinIdentity_pct

)

p2 <- sprintf(

    paste(

        "ATPase alpha (ATPF1A) was not recovered as a complete",
        "coding sequence from the T. equiperdum IVM-t1 genome assembly.",
        "Although tblastn identified the genomic locus and a partial",
        "genomic fragment corresponding to ATPF1A was recovered,",
        "translation yielded only 389 of the expected 584 amino acids",
        "(66.6%% coverage).",
        "Because ATPF1A contains multiple introns in Trypanosoma,",
        "reconstruction of the mature coding sequence requires",
        "exon-aware gene prediction or transcript evidence.",
        "Consequently, sequence conservation could not be reliably",
        "evaluated from the genome assembly alone.",
        "Although tblastn identified the genomic locus,",
        "translation of the recovered genomic fragment yielded",
        "only %d of the expected %d amino acids",
        "(%.1f%% coverage).",
        "Because ATPF1A is intron-containing in Trypanosoma,",
        "reconstruction of the mature coding sequence requires",
        "exon-aware gene prediction or transcript evidence.",
        "Consequently, ATPF1A sequence conservation",
        "was not evaluated further."

    ),

    a$RecoveredProtein_aa,
    a$FullLength_aa,
    a$Coverage_pct

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

    p1,
    "\n\n",
    p2,
    "\n"

)

###############################################################
## Save outputs
###############################################################

write.csv(
    summary_df,
    out_summary_csv,
    row.names = FALSE
)

write.csv(
    marker_df,
    out_markers_csv,
    row.names = FALSE
)

cat(
    manuscript_text,
    file = out_text_report
)

###############################################################
## Console output
###############################################################

cat(
    "======================================================\n"
)

cat(
    "Supplementary Table S1: ATPase Summary\n"
)

cat(
    "======================================================\n"
)

print(

    summary_df[
        ,
        c(
            "Gene",
            "SequenceSource",
            "Coverage_pct",
            "ProteinIdentity_pct",
            "Substitutions",
            "Status"
        )
    ]

)

cat(
    "\n======================================================\n"
)

cat(
    "Supplementary Table S2: ATPF1G Diagnostic Markers\n"
)

cat(
    "======================================================\n"
)

print(marker_df)

cat("\n", manuscript_text, "\n")

cat("Deliverables written to:\n")
cat(" -", out_summary_csv, "\n")
cat(" -", out_markers_csv, "\n")
cat(" -", out_text_report, "\n")