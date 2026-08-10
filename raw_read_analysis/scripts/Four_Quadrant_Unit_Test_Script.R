library(Biostrings)

# Load updated get_gap_facing_nt definition
source("raw_read_analysis/scripts/04_tac102_gap_sequence_evidence.R")

# Mock 100 bp synthetic read sequence: "A" x 30, "C" x 40, "G" x 30
# Position 1..30 = A, 31..70 = C, 71..100 = G
test_seq <- DNAString(paste0(
  paste(rep("A", 30), collapse = ""),
  paste(rep("C", 40), collapse = ""),
  paste(rep("G", 30), collapse = "")
))

test_cases <- list(
  # Q1: N-flank, Forward (sstart < send)
  "N_Forward" = list(
    side = "N",
    hsp = list(sstart_raw = 10, send_raw = 40),
    exp_start = 41, exp_end = 100, exp_len = 60,
    exp_first_base = "C", exp_rc = FALSE
  ),
  # Q2: N-flank, Reverse (sstart > send)
  "N_Reverse" = list(
    side = "N",
    hsp = list(sstart_raw = 40, send_raw = 10),
    exp_start = 1, exp_end = 9, exp_len = 9,
    exp_first_base = "T", exp_rc = TRUE # "A" RC'd to "T"
  ),
  # Q3: C-flank, Forward (sstart < send)
  "C_Forward" = list(
    side = "C",
    hsp = list(sstart_raw = 60, send_raw = 90),
    exp_start = 1, exp_end = 59, exp_len = 59,
    exp_first_base = "A", exp_rc = FALSE
  ),
  # Q4: C-flank, Reverse (sstart > send)
  "C_Reverse" = list(
    side = "C",
    hsp = list(sstart_raw = 90, send_raw = 60),
    exp_start = 91, exp_end = 100, exp_len = 10,
    exp_first_base = "C", exp_rc = TRUE # "G" RC'd to "C"
  )
)

cat("==============================================================\n")
cat(" RUNNING 4-QUADRANT UNIT TEST: get_gap_facing_nt()\n")
cat("==============================================================\n")

all_passed <- TRUE
for (name in names(test_cases)) {
  tc <- test_cases[[name]]
  res <- get_gap_facing_nt(test_seq, tc$hsp, tc$side)
  
  obs_len <- length(res$sequence)
  obs_first_base <- if (obs_len > 0) as.character(res$sequence[1]) else ""
  
  pass_start <- isTRUE(res$nt_start == tc$exp_start)
  pass_end   <- isTRUE(res$nt_end == tc$exp_end)
  pass_len   <- isTRUE(obs_len == tc$exp_len)
  pass_base  <- isTRUE(obs_first_base == tc$exp_first_base)
  
  tc_pass <- pass_start && pass_end && pass_len && pass_base
  all_passed <- all_passed && tc_pass
  
  status <- if (tc_pass) "PASS" else "FAIL"
  cat(sprintf("[%s] %s: Start=%d (exp %d), End=%d (exp %d), Len=%d (exp %d), 1st Base=%s (exp %s)\n",
              status, name, res$nt_start, tc$exp_start, res$nt_end, tc$exp_end, 
              obs_len, tc$exp_len, obs_first_base, tc$exp_first_base))
}

if (all_passed) {
  cat("\n[SUCCESS] All 4 quadrants passed mathematical validation.\n")
} else {
  stop("\n[FAILURE] One or more quadrant tests failed.")
}