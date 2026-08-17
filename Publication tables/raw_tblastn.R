# ============================================================
# Complete Raw tBLASTn Evidence Table
# Trypanosoma equiperdum IVM-t1
# ============================================================

library(tidyverse)
library(gt)
library(pagedown)

# ------------------------------------------------------------
# 1. Define input directory and genes
# ------------------------------------------------------------

results_dir <- "reports"

genes <- c(
  "ATPF1A",
  "ATPF1B",
  "ATPF1G",
  "ATOM40",
  "POLIB",
  "TAC102"
)

# ------------------------------------------------------------
# 2. Read all tBLASTn result files
# ------------------------------------------------------------

read_tblastn <- function(gene) {
  file <- file.path(
    results_dir,
    paste0(gene, "_tblastn_results.csv")
  )
  
  if (!file.exists(file)) {
    stop("Missing BLAST result file: ", file)
  }
  
  read_csv(file, show_col_types = FALSE) |>
    mutate(Gene = gene, .before = 1)
}

tblastn_raw <- map_dfr(genes, read_tblastn)

# ------------------------------------------------------------
# 3. Calculate derived statistics
# ------------------------------------------------------------

query_lengths <- tibble(
  Gene = genes,
  query_length_aa = c(519, 519, 305, 351, 1404, 951)
)

tblastn_evidence <- tblastn_raw |>
  left_join(query_lengths, by = "Gene") |>
  mutate(
    query_coverage = 100 * length / query_length_aa,
    query_coordinates = paste0(qstart, "–", qend),
    subject_coordinates = paste0(sstart, "–", send)
  )

# ------------------------------------------------------------
# 4. Order HSPs
# ------------------------------------------------------------

tblastn_evidence <- tblastn_evidence |>
  group_by(Gene) |>
  mutate(HSP = row_number()) |>
  ungroup()

# ------------------------------------------------------------
# 5. Select and rename columns for presentation
# ------------------------------------------------------------

tblastn_table_data <- tblastn_evidence |>
  select(
    Gene,
    HSP,
    qseqid,
    sseqid,
    length,
    query_coverage,
    pident,
    mismatch,
    gapopen,
    query_coordinates,
    subject_coordinates,
    evalue,
    bitscore
  ) |>
  rename(
    `Gene` = Gene,
    `HSP` = HSP,
    `Query protein` = qseqid,
    `IVM-t1 contig / accession` = sseqid,
    `Alignment length (aa)` = length,
    `Query coverage (%)` = query_coverage,
    `Identity (%)` = pident,
    `Mismatches` = mismatch,
    `Gap opens` = gapopen,
    `Query coordinates` = query_coordinates,
    `IVM-t1 coordinates` = subject_coordinates,
    `E-value` = evalue,
    `Bit score` = bitscore
  )
# ------------------------------------------------------------
# 6. Create publication-style gt table
# ------------------------------------------------------------

tblastn_hsp_gt <- tblastn_table_data |>
  gt(
    rowname_col = "HSP",
    groupname_col = "Gene"
  ) |>
  tab_stubhead(label = "HSP") |>
  tab_header(
    title = md("**Raw tBLASTn evidence underlying assembly-level gene recovery**"),
    subtitle = md("*Trypanosoma equiperdum* IVM-t1 mitochondrial-maintenance gene panel")
  ) |>
  cols_label(
    `Query protein`             = "Query protein",
    `IVM-t1 contig / accession` = md("IVM-t1<br>contig / accession"),
    `Alignment length (aa)`     = md("Alignment<br>length (aa)"),
    `Query coverage (%)`        = md("Query<br>coverage (%)"),
    `Identity (%)`              = md("Identity<br>(%)"),
    `Mismatches`                = "Mismatches",
    `Gap opens`                 = md("Gap<br>opens"),
    `Query coordinates`         = md("Query<br>coordinates"),
    `IVM-t1 coordinates`        = md("IVM-t1<br>coordinates"),
    `E-value`                   = "E-value",
    `Bit score`                 = md("Bit<br>score")
  ) |>
  fmt_number(
    columns = `Alignment length (aa)`,
    decimals = 0,
    use_seps = TRUE
  ) |>
  fmt_number(
    columns = `Query coverage (%)`,
    decimals = 1
  ) |>
  fmt_number(
    columns = `Identity (%)`,
    decimals = 2
  ) |>
  fmt_number(
    columns = c(Mismatches, `Gap opens`),
    decimals = 0,
    use_seps = TRUE
  ) |>
  fmt_scientific(
    columns = `E-value`,
    decimals = 2
  ) |>
  fmt_number(
    columns = `Bit score`,
    decimals = 1,
    use_seps = TRUE
  ) |>
  # Bold ATPF1A row group and identities
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups(groups = "ATPF1A")
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `Identity (%)`,
      rows = Gene == "ATPF1A"
    )
  ) |>
  # Highlight high-identity primary HSPs (Fixed column name)
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = c(`Identity (%)`, `Query coverage (%)`),
      rows = `Identity (%)` >= 90
    )
  ) |>
  # Styling options
  tab_options(
    table.font.names = "Arial",
    table.font.size = px(10.5),
    heading.title.font.size = px(15),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight = "bold",
    column_labels.font.size = px(10),
    data_row.padding = px(6),
    row_group.font.weight = "bold",
    row_group.font.size = px(11),
    table.border.top.width = px(1.5),
    table.border.bottom.width = px(1.5),
    column_labels.border.top.width = px(1),
    column_labels.border.bottom.width = px(1),
    row.striping.include_table_body = TRUE,
    row.striping.background_color = "#F7F7F7"
  ) |>
  cols_width(
    `Query protein`             ~ px(130),
    `IVM-t1 contig / accession` ~ px(120),
    `Alignment length (aa)`     ~ px(80),
    `Query coverage (%)`        ~ px(80),
    `Identity (%)`              ~ px(65),
    `Mismatches`                ~ px(60),
    `Gap opens`                 ~ px(55),
    `Query coordinates`         ~ px(85),
    `IVM-t1 coordinates`        ~ px(100),
    `E-value`                   ~ px(80),
    `Bit score`                 ~ px(65)
  ) |>
  tab_source_note(
    source_note = md(
      "**HSP**, high-scoring segment pair. Query coverage is calculated as `alignment length / full-length query protein × 100` for each HSP."
    )
  ) |>
  tab_source_note(
    source_note = md(
      "Coordinates are reported exactly as returned by tBLASTn. Because tBLASTn reports translated alignments, coordinate orientation may be decreasing for reverse-strand matches."
    )
  ) |>
  tab_source_note(
    source_note = md(
      "ATPF1A produced three low-identity partial HSPs (24.55–27.45%), whereas ATPF1B, ATPF1G, ATOM40 and POLIB each produced a near-complete/full-length high-identity primary match. TAC102 was represented by two high-identity HSPs covering residues 1–648 and 733–951 of the 951-aa query."
    )
  ) |>
  tab_source_note(
    source_note = md(
      "Source: tBLASTn searches of *T. equiperdum* IVM-t1 assembled genome sequences using *T. brucei* reference protein queries."
    )
  )

# ------------------------------------------------------------
# 7. Display & Export (fixed)
# ------------------------------------------------------------

# 1. Define output paths
html_path <- file.path(results_dir, "IVMt1_tBLASTn_HSP_evidence.html")
pdf_path  <- file.path(results_dir, "IVMt1_tBLASTn_HSP_evidence.pdf")

# 2. Save table as HTML
gtsave(tblastn_hsp_gt, html_path)

# 3. Convert HTML to Landscape PDF with explicit page size
chrome_print(
  input = html_path,
  output = pdf_path,
  options = list(
    landscape = TRUE,
    printBackground = TRUE,
    paperWidth  = 8.5,   # standard portrait values — Chrome will swap them
    paperHeight = 11,
    preferCSSPageSize = FALSE
  )
)

# 4. (Optional) PNG export
gtsave(
  tblastn_hsp_gt,
  file.path(results_dir, "IVMt1_tBLASTn_HSP_evidence.png"),
  vwidth = 1000,
  vheight = 1000,
  expand = 20
)