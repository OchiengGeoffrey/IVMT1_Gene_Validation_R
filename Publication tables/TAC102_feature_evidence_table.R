# ============================================================
# TAC102 feature evidence table – publication‑ready version
# ============================================================

library(tidyverse)
library(gt)
library(pagedown)

# ------------------------------------------------------------
# 1. Input data
# ------------------------------------------------------------

tac102_data <- tibble(
  `TAC102 feature` = c(
    "Residue 653",
    "Residue 698",
    "653–698 interval",
    "Paired-read evidence at 653",
    "Paired-read evidence at 698"
  ),
  Reference = c("K", "L", "—", "K", "L"),
  `Raw-read evidence` = c(
    "268 positional hits / 237 unique reads",
    "240 positional hits / 218 unique reads",
    "Reads recovered across both diagnostic positions",
    "Candidate pairs with positional evidence",
    "Candidate pairs with positional evidence"
  ),
  `Observed allele` = c(
    "M: 253; K: 3; other: 12",
    "L: 234; P: 3; other: 3",
    "M653–L698 observed",
    "Concordant M/M and discordant observations",
    "Predominantly L/L"
  ),
  `Read depth / observations` = c(
    "268 observations",
    "240 observations",
    "1+ reads spanning both positions",
    "Full paired-read dataset",
    "Full paired-read dataset"
  ),
  `Strand evidence` = c(
    "Fwd + Rev",
    "Fwd + Rev",
    "Fwd/Rev evidence",
    "Both mates",
    "Both mates"
  ),
  Quality = c(
    "High-quality codon support",
    "High-quality codon support",
    "Q27–34 in pilot; full scan available",
    "To be quantified exactly",
    "To be quantified exactly"
  ),
  Interpretation = c(
    "M653 is the dominant IVM-t1 sequence state",
    "L698 is predominantly reference-matching",
    "Supports sequence continuity across the assembly-unresolved region",
    "Independent mate-level support for dominant M653 state",
    "Supports reference-like L698 state"
  )
)

# ------------------------------------------------------------
# 2. Build the gt table with A4 Landscape responsive layout
# ------------------------------------------------------------

tac102_gt <- tac102_data |>
  gt() |>
  tab_header(
    title = md("**TAC102 diagnostic features: raw‑read evidence summary**"),
    subtitle = md("*Trypanosoma equiperdum* IVM‑t1")
  ) |>
  cols_label(
    `TAC102 feature`            = md("TAC102<br>feature"),
    Reference                   = "Ref",
    `Raw-read evidence`         = md("Raw‑read<br>evidence"),
    `Observed allele`           = md("Observed<br>allele"),
    `Read depth / observations` = md("Read depth /<br>observations"),
    `Strand evidence`           = md("Strand<br>evidence"),
    Quality                     = "Quality",
    Interpretation              = "Interpretation"
  ) |>
  cols_align(
    align = "left",
    columns = everything()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) |>
  tab_style(
    style = cell_fill(color = "#F0F0F0"),
    locations = cells_body(rows = seq(1, nrow(tac102_data), 2))
  ) |>
  fmt_markdown(columns = everything()) |>
  tab_options(
    table.width = pct(100),
    table.font.names = "Arial",
    table.font.size = px(9.5),
    heading.title.font.size = px(13),
    heading.subtitle.font.size = px(10.5),
    column_labels.font.weight = "bold",
    column_labels.font.size = px(9.5),
    data_row.padding = px(4),
    row.striping.include_table_body = FALSE,
    table.border.top.width = px(1.5),
    table.border.bottom.width = px(1.5),
    column_labels.border.top.width = px(1),
    column_labels.border.bottom.width = px(1)
  ) |>
  cols_width(
    `TAC102 feature`            ~ pct(12),
    Reference                   ~ pct(4),
    `Raw-read evidence`         ~ pct(17),
    `Observed allele`           ~ pct(15),
    `Read depth / observations` ~ pct(13),
    `Strand evidence`           ~ pct(8),
    Quality                     ~ pct(13),
    Interpretation              ~ pct(18)
  ) |>
  tab_source_note(
    source_note = md(
      "**Raw‑read evidence** summarises variant calls, read counts, and paired‑end support at two key residues (653 and 698) as well as the intervening interval."
    )
  ) |>
  tab_source_note(
    source_note = md(
      "**Reference** indicates the amino acid in the *T. brucei* reference sequence. **M653** and **L698** represent the dominant IVM‑t1 alleles."
    )
  ) |>
  tab_source_note(
    source_note = md(
      "Data derived from pilot sequencing; full scan of the TAC102 locus is available."
    )
  ) |>
  # Inject A4 Landscape print dimensions
  opt_css(
    css = "
    @page {
      size: A4 landscape;
      margin: 8mm;
    }
    body {
      margin: 0;
      padding: 0;
    }
    table {
      width: 100% !important;
    }
    "
  )

# ------------------------------------------------------------
# 3. Export PDF & PNG
# ------------------------------------------------------------

html_path <- "TAC102_feature_table.html"
pdf_path  <- "TAC102_feature_table.pdf"
png_path  <- "TAC102_feature_table.png"

# Save HTML intermediate
gtsave(tac102_gt, html_path)

# Render A4 Landscape PDF
chrome_print(
  input = html_path,
  output = pdf_path,
  options = list(
    landscape = TRUE,
    printBackground = TRUE,
    preferCSSPageSize = TRUE
  )
)

# Render high-res PNG fallback
gtsave(
  tac102_gt,
  png_path,
  vwidth = 1600,
  vheight = 900,
  expand = 10
)