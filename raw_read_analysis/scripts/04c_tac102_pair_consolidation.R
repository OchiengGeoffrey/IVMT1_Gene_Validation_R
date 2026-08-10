# ==============================================================
# Step 4B.3c: TAC102 pair-level evidence consolidation
# ==============================================================

tails_new <- read.csv(
  "raw_read_analysis/reports/tac102_41pairs_unaligned_tails.csv",
  stringsAsFactors = FALSE
)

gap_new <- read.csv(
  "raw_read_analysis/reports/tac102_41pairs_gap_evidence.csv",
  stringsAsFactors = FALSE
)

cat("==============================================================\n")
cat("Step 4B.3c: TAC102 pair-level evidence consolidation\n")
cat("==============================================================\n\n")

cat("Input rows:", nrow(gap_new), "\n")
cat("Unique pairs:", length(unique(gap_new$pair_id)), "\n\n")


# --------------------------------------------------------------
# 1. Confirm exactly one N and one C observation per pair
# --------------------------------------------------------------

role_table <- table(gap_new$pair_id, gap_new$orientation_role)

cat("=== N/C ROLE STRUCTURE ===\n")
print(table(rowSums(role_table)))

bad_pairs <- names(which(
  apply(role_table, 1, function(x) {
    !identical(as.integer(x[c("N", "C")]), c(1L, 1L))
  })
))

cat("\nPairs with abnormal N/C structure:", length(bad_pairs), "\n")

if (length(bad_pairs) > 0) {
  print(bad_pairs)
}


# --------------------------------------------------------------
# 2. Reshape N and C observations to one row per pair
# --------------------------------------------------------------

N <- gap_new[
  gap_new$orientation_role == "N",
  c(
    "pair_id",
    "mate",
    "hsp_strand",
    "gap_facing_nt_length",
    "gap_position_start",
    "gap_position_end",
    "aligned_aa",
    "pident_aa",
    "evidence_class"
  )
]

C <- gap_new[
  gap_new$orientation_role == "C",
  c(
    "pair_id",
    "mate",
    "hsp_strand",
    "gap_facing_nt_length",
    "gap_position_start",
    "gap_position_end",
    "aligned_aa",
    "pident_aa",
    "evidence_class"
  )
]

names(N)[-1] <- paste0("N_", names(N)[-1])
names(C)[-1] <- paste0("C_", names(C)[-1])

pairs <- merge(N, C, by = "pair_id", all = TRUE)

cat("\nPair-level rows:", nrow(pairs), "\n")


# --------------------------------------------------------------
# 3. Pair-level evidence classification
# --------------------------------------------------------------

pairs$pair_class <- with(
  pairs,
  ifelse(
    N_evidence_class == "A" & C_evidence_class == "A",
    "A/A",
    ifelse(
      N_evidence_class == "A" | C_evidence_class == "A",
      "A/other",
      ifelse(
        N_evidence_class == "D" & C_evidence_class == "D",
        "D/D",
        "other"
      )
    )
  )
)


# --------------------------------------------------------------
# 4. Determine whether each flank reaches the defined TAC102 gap
# --------------------------------------------------------------

pairs$N_reaches_gap <- !is.na(pairs$N_gap_position_start) &
                       !is.na(pairs$N_gap_position_end)

pairs$C_reaches_gap <- !is.na(pairs$C_gap_position_start) &
                       !is.na(pairs$C_gap_position_end)

pairs$both_reach_gap <- pairs$N_reaches_gap &
                        pairs$C_reaches_gap


# --------------------------------------------------------------
# 5. Calculate coordinate span represented by each pair
# --------------------------------------------------------------

pairs$pair_gap_start <- ifelse(
  pairs$N_reaches_gap & pairs$C_reaches_gap,
  pmin(
    pairs$N_gap_position_start,
    pairs$C_gap_position_start
  ),
  NA
)

pairs$pair_gap_end <- ifelse(
  pairs$N_reaches_gap & pairs$C_reaches_gap,
  pmax(
    pairs$N_gap_position_end,
    pairs$C_gap_position_end
  ),
  NA
)

pairs$pair_gap_span_aa <- ifelse(
  !is.na(pairs$pair_gap_start) &
  !is.na(pairs$pair_gap_end),
  pairs$pair_gap_end - pairs$pair_gap_start + 1,
  NA
)


# --------------------------------------------------------------
# 6. Determine whether the two flanks collectively cover
#    the entire defined gap 649-732
# --------------------------------------------------------------

pairs$full_gap_649_732 <- pairs$both_reach_gap &
                          pairs$pair_gap_start <= 649 &
                          pairs$pair_gap_end >= 732


# --------------------------------------------------------------
# 7. Report pair-level evidence
# --------------------------------------------------------------

cat("\n=== PAIR-LEVEL EVIDENCE CLASSES ===\n")
print(table(pairs$pair_class, useNA = "ifany"))

cat("\n=== N-FLANK GAP REACH ===\n")
print(table(pairs$N_reaches_gap, useNA = "ifany"))

cat("\n=== C-FLANK GAP REACH ===\n")
print(table(pairs$C_reaches_gap, useNA = "ifany"))

cat("\n=== BOTH FLANKS REACH GAP ===\n")
print(table(pairs$both_reach_gap, useNA = "ifany"))

cat("\n=== FULL 649-732 GAP COVERAGE ===\n")
print(table(pairs$full_gap_649_732, useNA = "ifany"))


# --------------------------------------------------------------
# 8. Print the key pair-level table
# --------------------------------------------------------------

cat("\n=== PAIR-LEVEL SUMMARY ===\n")

print(
  pairs[
    order(
      pairs$full_gap_649_732,
      pairs$pair_class,
      decreasing = TRUE
    ),
    c(
      "pair_id",
      "N_mate",
      "C_mate",
      "N_evidence_class",
      "C_evidence_class",
      "N_gap_position_start",
      "N_gap_position_end",
      "C_gap_position_start",
      "C_gap_position_end",
      "both_reach_gap",
      "full_gap_649_732",
      "pair_class"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


# --------------------------------------------------------------
# 9. Save pair-level results
# --------------------------------------------------------------

out_file <- paste0(
  "raw_read_analysis/reports/",
  "tac102_41pairs_pair_level_evidence.csv"
)

write.csv(
  pairs,
  out_file,
  row.names = FALSE
)

cat("\nSaved:\n", out_file, "\n")


# --------------------------------------------------------------
# 10. Compact final summary
# --------------------------------------------------------------

cat("\n==============================================================\n")
cat("4B.3c SUMMARY\n")
cat("==============================================================\n")

cat("Both-flank pairs: ", nrow(pairs), "\n", sep = "")
cat("A/A pairs: ", sum(pairs$pair_class == "A/A"), "\n", sep = "")
cat("A/other pairs: ", sum(pairs$pair_class == "A/other"), "\n", sep = "")
cat("D/D pairs: ", sum(pairs$pair_class == "D/D"), "\n", sep = "")
cat(
  "Both flanks reach gap: ",
  sum(pairs$both_reach_gap, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Pairs spanning 649-732: ",
  sum(pairs$full_gap_649_732, na.rm = TRUE),
  "\n",
  sep = ""
)

cat("\n[COMPLETE] Step 4B.3c pair-level consolidation finished.\n")