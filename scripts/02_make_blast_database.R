source("scripts/01_project_setup.R")

project_header("Creating BLAST database")

assembly <- PATHS$assembly

dir.create(PATHS$blast,
           recursive = TRUE,
           showWarnings = FALSE)

blast_db <- file.path(PATHS$blast, "IVMT1")

cmd <- paste(
  "makeblastdb",
  "-in", shQuote(assembly),
  "-dbtype nucl",
  "-out", shQuote(blast_db)
)

cat(cmd, "\n")

system(cmd)