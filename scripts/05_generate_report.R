dir.create("results/report", recursive = TRUE, showWarnings = FALSE)

writeLines(
  capture.output(sessionInfo()),
  "reports/sessionInfo.txt"
)