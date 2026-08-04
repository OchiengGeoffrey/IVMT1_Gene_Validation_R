dir.create(
    PATHS$reports,
    recursive = TRUE,
    showWarnings = FALSE
)

writeLines(
    capture.output(sessionInfo()),
    file.path(PATHS$reports,
              "sessionInfo.txt")
)