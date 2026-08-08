# ==============================================================================
# Helper Functions for Phase II Raw Read Validation
# ==============================================================================

#' Log an info message to console and file
log_info <- function(msg, log_file = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- sprintf("[%s] INFO: %s", timestamp, msg)
  message(full_msg)
  if (!is.null(log_file) && file.exists(log_file)) {
    cat(full_msg, "\n", file = log_file, append = TRUE)
  }
}

#' Log an error message to console and file
log_error <- function(msg, log_file = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- sprintf("[%s] ERROR: %s", timestamp, msg)
  message(full_msg)
  if (!is.null(log_file) && file.exists(log_file)) {
    cat(full_msg, "\n", file = log_file, append = TRUE)
  }
}

#' Log a warning message to console and file
log_warn <- function(msg, log_file = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- sprintf("[%s] WARN: %s", timestamp, msg)
  warning(full_msg, call. = FALSE)
  if (!is.null(log_file) && file.exists(log_file)) {
    cat(full_msg, "\n", file = log_file, append = TRUE)
  }
}

#' Detect if FASTQ exists as .fastq or .fastq.gz
detect_fastq <- function(base_path, log_file = NULL) {
  candidates <- c(
    paste0(base_path, ".fastq"),
    paste0(base_path, ".fq"),
    paste0(base_path, ".fastq.gz"),
    paste0(base_path, ".fq.gz")
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    log_error(sprintf("FASTQ file not found: %s (tried %s)", base_path, paste(candidates, collapse=", ")), log_file)
    stop("Missing FASTQ file")
  }
  if (length(found) > 1) {
    log_warn(sprintf("Multiple FASTQ variants found for %s, using %s", base_path, found[1]), log_file)
  }
  return(found[1])
}

#' Verify write permissions on directories
verify_write_permissions <- function(paths, log_file = NULL) {
  for (name in names(paths)) {
    p <- paths[[name]]
    # PATHS also contains input FASTQ files.  They must be readable, but are
    # not output directories and should never be created as directories.
    if (file.exists(p) && !dir.exists(p)) {
      if (file.access(p, mode = 4) != 0) {
        log_error(sprintf("No read permission for input file: %s", p), log_file)
        stop(sprintf("Permission denied: %s", p))
      }
      next
    }
    if (!dir.exists(p)) {
      log_info(sprintf("Creating directory: %s", p), log_file)
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
    if (!file.access(p, mode = 2) == 0) {
      log_error(sprintf("No write permission for directory: %s", p), log_file)
      stop(sprintf("Permission denied: %s", p))
    }
  }
  log_info("All directories have write permission", log_file)
}

#' Verify that input FASTQ files exist and all reference FASTA files exist
verify_inputs_and_fastas <- function(paths, targets, log_file = NULL) {
  # Check FASTQs
  if (!file.exists(paths$fastq_r1)) {
    log_error("R1 FASTQ not found", log_file)
    stop("R1 FASTQ not found")
  }
  if (!file.exists(paths$fastq_r2)) {
    log_error("R2 FASTQ not found", log_file)
    stop("R2 FASTQ not found")
  }
  log_info(sprintf("FASTQ R1: %s", paths$fastq_r1), log_file)
  log_info(sprintf("FASTQ R2: %s", paths$fastq_r2), log_file)
  
  # Check each reference protein FASTA
  for (i in seq_len(nrow(targets))) {
    f <- targets$reference_file[i]
    if (!file.exists(f)) {
      log_error(sprintf("Reference FASTA missing: %s", f), log_file)
      stop(sprintf("Missing reference: %s", f))
    }
    log_info(sprintf("Reference FASTA OK: %s", basename(f)), log_file)
  }
}

#' Check if an external program is available (optional)
check_program <- function(prog, required = FALSE, log_file = NULL) {
  cmd <- Sys.which(prog)
  if (cmd == "") {
    msg <- sprintf("Program '%s' not found in PATH.", prog)
    if (required) {
      log_error(msg, log_file)
      stop(msg)
    } else {
      log_info(sprintf("%s (optional, skipping)", msg), log_file)
    }
  } else {
    log_info(sprintf("Program '%s' found at: %s", prog, cmd), log_file)
  }
}

#' Log system memory (Linux/Mac) – fallback for Windows
log_system_memory <- function(log_file = NULL) {
  if (Sys.info()["sysname"] == "Windows") {
    # WMIC is removed from current Windows releases. PowerShell CIM is its
    # supported replacement and is available on standard Windows installs.
    ps_command <- paste(
      "Get-CimInstance -ClassName Win32_OperatingSystem |",
      "Select-Object TotalVisibleMemorySize, FreePhysicalMemory |",
      "ConvertTo-Json -Compress"
    )
    mem <- tryCatch(
      system2(
        "powershell.exe",
        args = c("-NoProfile", "-NonInteractive", "-Command", shQuote(ps_command)),
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(mem) && length(mem) > 0 && !any(grepl("not recognized|error", mem, ignore.case = TRUE))) {
      log_info(paste("Memory info (Windows):", paste(mem, collapse = " ")), log_file)
    } else {
      log_warn("Unable to collect Windows memory information; continuing setup.", log_file)
    }
  } else {
    mem <- tryCatch(system2("free", "-h", stdout = TRUE, stderr = TRUE), error = function(e) NULL)
    if (!is.null(mem) && length(mem) > 0) {
      log_info(paste("Memory info:", paste(mem, collapse = " ")), log_file)
    } else {
      log_warn("Unable to collect system memory information; continuing setup.", log_file)
    }
  }
}
