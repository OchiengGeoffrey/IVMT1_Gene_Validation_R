# ============================================================
# Paired-read support for TAC102 positional observations table
# ============================================================

library(tidyverse)
library(gt)
library(pagedown)

# ------------------------------------------------------------
# 1. Input data
# ------------------------------------------------------------

paired_read_data <- tibble::tribble(
  ~Position, ~`Both mates concordant`, ~`Both mates discordant`, ~`Dual-evidence pairs`, ~Concordance, ~`R1-only`, ~`R2-only`,
  653,       27,                       4,                        31,                     0.871,        105,        101,
  698,       20,                       2,                        22,                     0.909,        100,        96
)

# ------------------------------------------------------------
# 2. Build the gt table
# ------------------------------------------------------------

paired_read_gt <- paired_read_data %>%
  gt() %>%

  tab_header(
    title = md("**Paired-read support for TAC102 positional observations**")
  ) %>%

  cols_label(
    Position                = "Position",
    `Both mates concordant` = md("Both mates<br>concordant"),
    `Both mates discordant` = md("Both mates<br>discordant"),
    `Dual-evidence pairs`   = md("Dual-evidence<br>pairs"),
    Concordance             = "Concordance",
    `R1-only`               = "R1-only",
    `R2-only`               = "R2-only"
  ) %>%

  cols_align(
    align = "center",
    columns = everything()
  ) %>%

  fmt_percent(
    columns = Concordance,
    decimals = 1
  ) %>%

  fmt_number(
    columns = c(Position, `Both mates concordant`, `Both mates discordant`, `Dual-evidence pairs`, `R1-only`, `R2-only`),
    decimals = 0,
    use_seps = FALSE
  ) %>%

  cols_width(
    Position                ~ pct(10),
    `Both mates concordant` ~ pct(18),
    `Both mates discordant` ~ pct(18),
    `Dual-evidence pairs`   ~ pct(18),
    Concordance             ~ pct(12),
    `R1-only`               ~ pct(12),
    `R2-only`               ~ pct(12)
  ) %>%

  tab_options(
    table.width = pct(100),
    table.font.names = "Arial",
    table.font.size = px(11.5),
    heading.title.font.size = px(14),
    heading.align = "center",             # Centered title
    column_labels.font.weight = "bold",
    column_labels.font.size = px(11),
    data_row.padding = px(8),
    table.border.top.width = px(0),
    table.border.bottom.width = px(1.5),
    column_labels.border.top.width = px(1),
    column_labels.border.bottom.width = px(1)
  ) %>%

  tab_source_note(
    source_note = md(
      "Paired-read concordance was assessed among read pairs for which both mates independently yielded positional evidence at the indicated TAC102 residue. Concordant pairs carried the same amino-acid call in both mates, whereas discordant pairs carried different calls. Dual-evidence pairs denote pairs in which both mates produced a validated positional observation. R1-only and R2-only counts represent read pairs for which positional evidence was recovered from only the indicated mate. Concordance percentages were calculated as the number of concordant pairs divided by the total number of dual-evidence pairs (concordant + discordant). Paired-read support provides independent mate-level evidence but does not by itself establish genomic fixation, zygosity, or population allele frequency."
    )
  ) %>%

  opt_css(
    css = "
    @page {
      size: A4 portrait;
      margin: 15mm;
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
# 3. Export
# ------------------------------------------------------------

html_path <- "TAC102_paired_read_support.html"
pdf_path  <- "TAC102_paired_read_support.pdf"
png_path  <- "TAC102_paired_read_support.png"

gtsave(paired_read_gt, html_path)

chrome_print(
  input = html_path,
  output = pdf_path,
  options = list(
    landscape = FALSE,
    printBackground = TRUE,
    preferCSSPageSize = TRUE
  )
)

gtsave(
  paired_read_gt,
  png_path,
  vwidth = 850,
  vheight = 400,
  expand = 10
)