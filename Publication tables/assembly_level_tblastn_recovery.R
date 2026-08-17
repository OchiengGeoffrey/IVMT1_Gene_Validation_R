# ============================================================
# Publication Table: Assembly-level tBLASTn recovery
# Trypanosoma equiperdum IVM-t1 mitochondrial-maintenance panel
# ============================================================

library(tidyverse)
library(gt)

# ------------------------------------------------------------
# 1. Enter the primary / strongest tBLASTn evidence
# ------------------------------------------------------------

tblastn_panel <- tibble(
  Gene = c(
    "ATPF1A",
    "ATPF1B",
    "ATPF1G",
    "ATOM40",
    "POLIB",
    "TAC102"
  ),

  `Functional role` = c(
    "F₁-ATPase α subunit",
    "F₁-ATPase β subunit",
    "F₁-ATPase γ subunit",
    "Mitochondrial outer-membrane protein",
    "Mitochondrial DNA polymerase",
    "Tripartite attachment complex protein"
  ),

  `Best hit` = c(
    "CM010682.1",
    "CM010674.1",
    "CM010681.1",
    "CM010680.1",
    "QSBY01000015.1",
    "CM010678.1"
  ),

  `Alignment (aa)` = c(
    "404",
    "519",
    "305",
    "351",
    "1,404",
    "648 + 219"
  ),

  `Query coverage (%)` = c(
    77.8,
    100.0,
    100.0,
    100.0,
    100.0,
    91.2
  ),

  `Identity (%)` = c(
    25.99,
    100.00,
    99.344,
    100.00,
    99.929,
    99.383
  ),

  `E-value` = c(
    1.67e-20,
    0,
    0,
    0,
    0,
    0
  ),

  `Bit score` = c(
    96.3,
    969,
    625,
    700,
    2882,
    1300
  ),

  `Assembly interpretation` = c(
    "Anomalous / unresolved",
    "Strong full-length recovery",
    "Strong full-length recovery",
    "Strong full-length recovery",
    "Strong full-length recovery",
    "Strong high-identity recovery; incomplete reference coverage"
  )
)


# ------------------------------------------------------------
# 2. Format the table
# ------------------------------------------------------------

tblastn_gt <- tblastn_panel |>

  gt() |>

  tab_header(
    title = md(
      "**Assembly-level recovery of the IVM-t1 mitochondrial-maintenance gene panel**"
    ),
    subtitle = md(
      "*tBLASTn comparison of T. brucei reference proteins against the assembled T. equiperdum IVM-t1 genome*"
    )
  ) |>

  # ----------------------------------------------------------
  # Column labels
  # Use html() so <br> is interpreted as a line break
  # ----------------------------------------------------------

  cols_label(
    Gene = "Gene",
    `Functional role` = "Functional role",
    `Best hit` = "Best IVM-t1 hit",
    `Alignment (aa)` = html("Alignment<br>(aa)"),
    `Query coverage (%)` = html("Query<br>coverage (%)"),
    `Identity (%)` = html("Identity<br>(%)"),
    `E-value` = "E-value",
    `Bit score` = html("Bit<br>score"),
    `Assembly interpretation` = "Assembly-level interpretation"
  ) |>

  # ----------------------------------------------------------
  # Number formatting
  # vars() helps with non-standard column names
  # ----------------------------------------------------------

  fmt_number(
    columns = vars(`Query coverage (%)`),
    decimals = 1
  ) |>

  fmt_number(
    columns = vars(`Identity (%)`),
    decimals = 2
  ) |>

  fmt_scientific(
    columns = vars(`E-value`),
    decimals = 2
  ) |>

  fmt_number(
    columns = vars(`Bit score`),
    decimals = 1
  ) |>

  # ----------------------------------------------------------
  # Emphasize ATPF1A anomaly
  # ----------------------------------------------------------

  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = vars(Gene),
      rows = Gene == "ATPF1A"
    )
  ) |>

  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = vars(`Identity (%)`),
      rows = Gene == "ATPF1A"
    )
  ) |>

  # ----------------------------------------------------------
  # Emphasize strong control loci
  # ----------------------------------------------------------

  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = vars(`Identity (%)`),
      rows = Gene != "ATPF1A"
    )
  ) |>

  # ----------------------------------------------------------
  # Table-wide styling
  # Corrected: table.font.name, not table.font.names
  # ----------------------------------------------------------

  tab_options(
    table.font.name = "Arial",
    table.font.size = px(9),
    table.width = pct(100),          # fill the page width, never overflow

    heading.title.font.size = px(12),
    heading.subtitle.font.size = px(9),

    column_labels.font.weight = "bold",
    column_labels.font.size = px(8.5),

    data_row.padding = px(4),

    table.border.top.width = px(1.2),
    table.border.bottom.width = px(1.2),
    column_labels.border.top.width = px(1),
    column_labels.border.bottom.width = px(1),

    row.striping.include_table_body = TRUE,
    row.striping.background_color = "#F7F7F7",

    source_notes.font.size = px(8)
  ) |>
  # ----------------------------------------------------------
  # Column widths
  # ----------------------------------------------------------

  cols_width(
    vars(Gene) ~ pct(7),
    vars(`Functional role`) ~ pct(16),
    vars(`Best hit`) ~ pct(12),
    vars(`Alignment (aa)`) ~ pct(8),
    vars(`Query coverage (%)`) ~ pct(9),
    vars(`Identity (%)`) ~ pct(7),
    vars(`E-value`) ~ pct(9),
    vars(`Bit score`) ~ pct(7),
    vars(`Assembly interpretation`) ~ pct(25)
  ) |>

  # ----------------------------------------------------------
  # Caption / source note
  # ----------------------------------------------------------

  tab_source_note(
    source_note = md(
      "**Note:** Query coverage is calculated relative to the full-length reference protein."
    )
  ) |>

  tab_source_note(
    source_note = md(
      "ATPF1A lacked a convincing full-length, high-identity genomic match; its strongest alignment was only 25.99% identical over 404 aa. In contrast, ATPF1B, ATPF1G, ATOM40 and POLIB showed near-complete or complete high-identity recovery. TAC102 was represented by two high-identity alignments covering 867 of 951 aa (91.2% combined query coverage)."
    )
  ) |>

  tab_source_note(
    source_note = md(
      "tBLASTn, translated nucleotide BLAST; IVM-t1, *Trypanosoma equiperdum* IVM-t1."
    )
  )


# ------------------------------------------------------------
# 3. Display
# ------------------------------------------------------------

tblastn_gt


# ------------------------------------------------------------
# 4. Export
# ------------------------------------------------------------

# Create output directory if it does not already exist
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

gtsave(
  tblastn_gt,
  "reports/IVMt1_tBLASTn_gene_panel_A4.pdf",
  vwidth = 794, vheight = 1123
)

# PNG export may require the webshot2 package and a Chromium/Chrome browser.
# If PNG export fails, try:
# install.packages("webshot2")
# webshot2::install_chrome()
gtsave(
  tblastn_gt,
  "reports/IVMt1_tBLASTn_gene_panel_A4.png",
  vwidth = 794, zoom = 3
)