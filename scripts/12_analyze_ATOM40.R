source("scripts/01_project_setup.R")

project_header("Verify ATOM40 Presence")

dir.create(PATHS$reports,   recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$extracted, recursive = TRUE, showWarnings = FALSE)

result <- recover_gene_from_tblastn(
  reference_protein = PATHS$atom40,
  assembly          = PATHS$assembly,
  blast_database    = file.path(PATHS$blast, "IVMT1"),
  gene_name         = "ATOM40",
  output_directory  = PATHS$extracted
)

cat("\nRecovered ORF statistics:\n")
cat("-------------------------\n")
cat("Frame          :", result$statistics$frame, "\n")
cat("Orientation    :", result$statistics$orientation, "\n")
cat("Identity       :", round(result$statistics$identity, 2), "%\n")
cat("Coverage       :", round(result$statistics$coverage, 2), "%\n")
cat("ORF length     :", result$statistics$orf_length, "aa\n")
cat("Internal stops :", result$statistics$internal_stops, "\n")
cat("Score          :", round(result$statistics$score, 4), "\n\n")

reference_protein <- read_protein_fasta(PATHS$atom40)

summary_row <- summarize_gene(
  gene                 = "ATOM40",
  ref_protein          = reference_protein,
  qry_protein          = result$protein,
  sequence_source      = paste0(
    result$statistics$blast_contig,
    " (",
    result$statistics$blast_strand,
    " strand; ORF frame ",
    result$statistics$frame,
    ")"
  ),
  status               = "ORF recovered via six-frame translation",
  evaluate_start_codon = FALSE,
  notes                = paste0(
    "ORF score=",
    round(result$statistics$score, 4),
    "; alignment identity=",
    round(result$statistics$identity, 2),
    "%; reference coverage=",
    round(result$statistics$coverage, 2),
    "%"
  )
)

save_blast_results(
  result$tblastn,
  file.path(PATHS$reports, "ATOM40_tblastn_results.csv")
)

cat("Summary:\n")
print(summary_row, row.names = FALSE)

cat("\nOutputs written:\n")
cat("  Genomic :", file.path(PATHS$extracted, "ATOM40_genomic.fasta"), "\n")
cat("  Protein :", file.path(PATHS$extracted, "ATOM40_protein.fasta"), "\n")
cat("  BLAST   :", file.path(PATHS$extracted, "ATOM40_tblastn.tsv"), "\n")
cat("  Report  :", file.path(PATHS$reports, "ATOM40_tblastn_results.csv"), "\n")
