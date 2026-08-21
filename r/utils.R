# Shared utility functions for the JABS nextflow postprocess pipeline.
# Source this file from other scripts: source("r/utils.R")

# Compare two vectors of NetworkFilenames; return a data frame of mismatches,
# or NULL if the sets are equal.
compare_NetworkFilenames <- function(x, y) {
  if (setequal(x, y)) return(NULL)

  x.notin.y <- tibble(
    NetworkFilename = setdiff(x, y),
    !!paste0(deparse(substitute(x)), "_not_in_", deparse(substitute(y))) := 1
  )
  y.notin.x <- tibble(
    NetworkFilename = setdiff(y, x),
    !!paste0(deparse(substitute(y)), "_not_in_", deparse(substitute(x))) := 1
  )

  out <- full_join(x.notin.y, y.notin.x, by = "NetworkFilename")
  out[is.na(out)] <- 0
  out
}

# Stop with a clear message if a required file does not exist.
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste0("Required file not found: '", file_path, "'"), call. = FALSE)
  }
}

# Stop with a clear message if a required directory does not exist.
check_dir_exists <- function(dir_path, hint = NULL) {
  if (!dir.exists(dir_path)) {
    msg <- paste0("Required directory not found: '", dir_path, "'")
    if (!is.null(hint)) msg <- paste0(msg, "\n", hint)
    stop(msg, call. = FALSE)
  }
}
