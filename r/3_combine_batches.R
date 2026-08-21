library(tidyverse)
library(data.table)
library(janitor)
options(error = NULL)

source("r/utils.R")

# ---- Configuration -----------------------------------------------------------

# Set working directory if not using an R project:
# working.directory <- "/path/to/project"
if (exists("working.directory")) setwd(working.directory)

final_feature_dir  <- file.path("Nextflow_Output", "final_nextflow_feature_data")
qc_missing_dup_dir <- file.path("qc", "missing_or_dup_data")
metadata_file      <- "metadata.csv"
exclude_file       <- "videos_to_exclude.txt"

# ---- Functions ---------------------------------------------------------------

validate_inputs <- function(feature_dir, qc_dir) {
  check_dir_exists(
    feature_dir,
    hint = "Run qc_check.R first to generate the feature files."
  )
  for (f in c("fecal_boli_final.csv", "gait_final.csv",
              "JABS_features_final.csv", "morphometrics_final.csv")) {
    check_file_exists(file.path(feature_dir, f))
  }
  check_dir_exists(
    qc_dir,
    hint = "Run qc_check.R first to generate the qc/ directory."
  )
}


load_feature_files <- function(feature_dir) {
  read_one <- function(pattern) {
    list.files(feature_dir, pattern = pattern, full.names = TRUE) |> read_csv()
  }
  list(
    gait       = read_one("gait_final"),
    morpho     = read_one("morphometrics_final") |> relocate(NetworkFilename),
    fecal_boli = read_one("fecal_boli_final")    |> relocate(NetworkFilename),
    jabs       = read_one("features_final")
  )
}


check_networkfilenames <- function(features, out_dir) {
  g <- features$gait$NetworkFilename
  m <- features$morpho$NetworkFilename
  f <- features$fecal_boli$NetworkFilename
  j <- features$jabs$NetworkFilename

  mismatches <- Filter(Negate(is.null), list(
    compare_NetworkFilenames(g, m),
    compare_NetworkFilenames(g, f),
    compare_NetworkFilenames(g, j),
    compare_NetworkFilenames(m, f),
    compare_NetworkFilenames(m, j),
    compare_NetworkFilenames(f, j)
  ))

  if (length(mismatches) > 0) {
    result <- Reduce(function(x, y) merge(x, y, by = "NetworkFilename", all = TRUE), mismatches)
  } else {
    result <- "NETWORKFILE NAMES MATCH PERFECTLY ACROSS DATAFRAMES"
  }

  write.csv(result, file.path(out_dir, "NetworkFilenames_missing_in_data.csv"), row.names = FALSE)
  result
}


merge_features <- function(features, exclude_file) {
  df <- features$gait |>
    merge(features$jabs,       by = "NetworkFilename") |>
    merge(features$morpho,     by = "NetworkFilename") |>
    merge(features$fecal_boli, by = "NetworkFilename")

  excluded_videos_present <- FALSE
  if (file.exists(exclude_file)) {
    to_exclude <- read.table(exclude_file, quote = "\"", comment.char = "")
    if (nrow(to_exclude) > 0) {
      df <- df[!df$NetworkFilename %in% to_exclude$V1, ]
      excluded_videos_present <- TRUE
    }
  }

  list(data = df, excluded_videos_present = excluded_videos_present)
}


merge_with_metadata <- function(df, metadata_file, out_dir) {
  check_file_exists(metadata_file)
  metadata <- read_csv(metadata_file, col_types = cols(MouseID = col_character()))

  df <- df |>
    mutate(
      FileName = str_split_i(NetworkFilename, "/", i = 4),
      FileName = gsub("_trimmed.avi", "", FileName),
      MouseID  = str_split_i(FileName, "_", i = 1)
    ) |>
    relocate(FileName)

  metadata_check <- compare_NetworkFilenames(df$MouseID, metadata$MouseID)
  if (is.null(metadata_check)) {
    metadata_check <- "NO MICE MISSING IN METADATA & VICE VERSA"
  }
  write.csv(metadata_check, file.path(out_dir, "mice_missing_in_metadata.csv"), row.names = FALSE)

  list(
    data           = merge(metadata, df, by = "MouseID"),
    metadata_check = metadata_check
  )
}


report_errors <- function(networkfilename_check, metadata_check, excluded_videos_present) {
  errors <- c(
    if (!is.character(networkfilename_check))
      "NetworkFilenames do not match perfectly across feature files — see NetworkFilenames_missing_in_data.csv",
    if (!is.character(metadata_check))
      "Mice in metadata and data do not match perfectly — see mice_missing_in_metadata.csv",
    if (!excluded_videos_present)
      "No videos_to_exclude.txt found (or file is empty) — is manual QC complete?"
  )

  if (length(errors) == 0) {
    message("No errors to report.")
  } else {
    message(paste(c("Warnings:", errors), collapse = "\n  "))
  }
}

# ---- Main --------------------------------------------------------------------

validate_inputs(final_feature_dir, qc_missing_dup_dir)

features <- load_feature_files(final_feature_dir)

networkfilename_check <- check_networkfilenames(features, qc_missing_dup_dir)

merged <- merge_features(features, exclude_file)

result <- merge_with_metadata(merged$data, metadata_file, qc_missing_dup_dir)

write_csv(result$data, file.path(final_feature_dir, "merged_nextflow_dataset.csv"))

report_errors(networkfilename_check, result$metadata_check, merged$excluded_videos_present)
