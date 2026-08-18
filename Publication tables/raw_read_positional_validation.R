# ============================================================
# TAC102 positional raw-read validation table (Sized for Portrait)
# ============================================================

library(tidyverse)
library(gt)
library(pagedown)

# ------------------------------------------------------------
# 1. Input data and percentage calculations
# ------------------------------------------------------------

tac102_depth <- tibble::tribble(
  ~Position, ~Reference, ~Validated, ~K, ~M, ~L, ~P, ~Other, ~Major_state,
  653,       "K",        268,       3, 253, NA, NA, 12,     "M",
  698,       "L",        240,      NA,  NA, 234,  3,  3,     "L"
) %>%
  mutate(
    across(
      c(K, M, L, P),
      ~ if_else(is.na(.x), "—", sprintf("%d (%.1f%%)", .x, 100 * .x / Validated))
    ),
    Other = sprintf("%d (%.1f%%)", Other, 100 * Other / Validated)
  )

# ------------------------------------------------------------
# 2. Build the gt table
# ------------------------------------------------------------

tac102_gt <- tac102_depth %>%
  gt() %>%

  tab_header(
    title = md(
      "**TAC102 positional raw-read validation in *T. equiperdum* IVM-t1**"
    ),
    subtitle = md(
      "Complete paired-end FASTQ scan; validated positional observations"
    )
  ) %>%

  cols_label(
    Position    = md("TAC102<br>position"),
    Reference   = "Reference",
    Validated   = md("Validated<br>observations"),
    K           = "K",
    M           = "M",
    L           = "L",
    P           = "P",
    Other       = "Other",
    Major_state = md("Dominant<br>state")
  ) %>%

  cols_align(
    align = "center",
    columns = everything()
  ) %>%

  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Major_state
    )
  ) %>%

  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = Position == 653,
      columns = M
    )
  ) %>%

  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = Position == 698,
      columns = L
    )
  ) %>%

  # Column widths tuned to eliminate excess horizontal stretching
  cols_width(
    Position    ~ px(85),
    Reference   ~ px(75),
    Validated   ~ px(100),
    K           ~ px(85),
    M           ~ px(85),
    L           ~ px(85),
    P           ~ px(85),
    Other       ~ px(90),
    Major_state ~ px(85)
  ) %>%

  tab_options(
    table.width = pct(100),
    table.font.names = "Arial",
    table.font.size = px(13),             # Increased from 12px for better readability
    heading.title.font.size = px(16),     # Increased title scale
    heading.subtitle.font.size = px(12),
    column_labels.font.weight = "bold",
    column_labels.font.size = px(12),
    data_row.padding = px(10),            # Increased padding to give rows vertical prominence
    table.border.top.width = px(1.5),
    table.border.bottom.width = px(1.5),
    column_labels.border.top.width = px(1),
    column_labels.border.bottom.width = px(1)
  ) %>%

  tab_source_note(
    source_note = md(
      "**Note:** Percentages are calculated among validated positional observations, not conventional nucleotide sequencing depth or formal variant allele frequencies. The reference residues are K653 and L698."
    )
  ) %>%

  # CSS page rules for A4 Portrait print rendering
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
      margin: 0 auto;
    }
    "
  )

# ------------------------------------------------------------
# 3. Export
# ------------------------------------------------------------

html_path <- "TAC102_positional_validation.html"
pdf_path  <- "TAC102_positional_validation.pdf"
png_path  <- "TAC102_positional_validation.png"

# Save intermediate HTML
gtsave(tac102_gt, html_path)

# PDF Render (A4 Portrait)
chrome_print(
  input = html_path,
  output = pdf_path,
  options = list(
    landscape = FALSE,
    printBackground = TRUE,
    preferCSSPageSize = TRUE
  )
)

# PNG Render (Narrower 850px viewport prevents overly stretched columns)
gtsave(
  tac102_gt,
  png_path,
  vwidth = 850,
  vheight = 450,
  expand = 10
)