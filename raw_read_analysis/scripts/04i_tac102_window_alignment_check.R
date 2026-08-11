library(dplyr)
library(readr)
library(stringr)

cat("==============================================================\n")
cat(" Step 4B.3i: Targeted Read Window & Phase Alignment Check\n")
cat("==============================================================\n\n")

# ==============================================================
# 1. PATHS & DATA LOADING
# ==============================================================

csv_tails      <- "raw_read_analysis/reports/tac102_41pairs_unaligned_tails.csv"
csv_04h_55     <- "raw_read_analysis/reports/tac102_4B3h_55_anomaly_audit.csv"

if (!file.exists(csv_tails)) {
  stop("Missing tail sequence report: ", csv_tails)
}

tails_df <- read_csv(csv_tails, show_col_types = FALSE)

# ==============================================================
# 2. HELPER: WINDOW & PHASE CALCULATOR
# ==============================================================

audit_row_window <- function(row) {
  # Extract raw alignment coordinates
  qstart_raw <- as.integer(row$hsp_qstart)
  sstart_raw <- as.integer(row$hsp_sstart)
  send_raw   <- as.integer(row$hsp_send)
  
  # Strand orientation
  forward    <- ifelse(!is.na(sstart_raw) && !is.na(send_raw), sstart_raw < send_raw, TRUE)
  
  p_first <- as.integer(row$gap_position_start)
  p_last  <- as.integer(row$gap_position_end)
  
  # Recompute codon-anchored bounds (extract_codon_at_reference_position logic)
  if (forward) {
    nt_codon_start <- sstart_raw + (p_first - qstart_raw) * 3L
    nt_codon_end   <- sstart_raw + (p_last  - qstart_raw) * 3L + 2L
  } else {
    nt_codon_end   <- sstart_raw - (p_first - qstart_raw) * 3L
    nt_codon_start <- sstart_raw - (p_last  - qstart_raw) * 3L - 2L
  }
  
  # Gap-facing block bounds (get_gap_facing_nt logic)
  gf_start <- as.integer(row$gap_facing_nt_start)
  gf_end   <- as.integer(row$gap_facing_nt_end)
  
  # Coordinate relationship checks
  is_nested <- !is.na(gf_start) && !is.na(gf_end) && 
               (nt_codon_start >= gf_start) && (nt_codon_end <= gf_end)
  
  overlap_len <- ifelse(
    is.na(gf_start) || is.na(gf_end), 0L,
    max(0L, min(nt_codon_end, gf_end) - max(nt_codon_start, gf_start) + 1L)
  )
  
  # Phase offset relative to gap-facing sequence start
  phase_offset <- if (!is.na(gf_start) && !is.na(nt_codon_start)) {
    (nt_codon_start - gf_start) %% 3L
  } else {
    NA_integer_
  }

  tibble(
    pair_id          = row$pair_id,
    orientation_role = row$orientation_role,
    strand           = ifelse(forward, "+", "-"),
    gap_aa_range     = paste0(p_first, "-", p_last),
    codon_nt_window  = paste0(nt_codon_start, "-", nt_codon_end),
    gap_facing_window= paste0(gf_start, "-", gf_end),
    is_nested        = is_nested,
    overlap_bp       = overlap_len,
    phase_offset_mod3= phase_offset
  )
}

# ==============================================================
# 3. RUN AUDIT ON SAMPLE ROWS
# ==============================================================

# Select sample rows across N and C flanks
sample_rows <- tails_df %>%
  filter(!is.na(gap_facing_sequence), gap_facing_sequence != "") %>%
  head(10)

cat("--- Inspecting 10 Sample Tail Records ---\n\n")

results_list <- lapply(seq_len(nrow(sample_rows)), function(i) {
  audit_row_window(sample_rows[i, ])
})

results_df <- bind_rows(results_list)

print(results_df, width = Inf)

# Summary of Phase Shifts
cat("\n--- Phase Offset Distribution (Modulo 3) ---\n")
print(table(results_df$phase_offset_mod3, useNA = "ifany"))

cat("\n--- Window Nesting Summary ---\n")
print(table(results_df$is_nested, useNA = "ifany"))

cat("\n[COMPLETE] Targeted window and phase check finished.\n")