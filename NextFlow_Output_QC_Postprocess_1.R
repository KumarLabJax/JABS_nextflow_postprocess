#!/usr/bin/env Rscript
# A script to clean NextFlow outputs for final processing
# Developed by Dr. Jake Beierle (don't forget the Dr., it's important)

# ----Documentation----
# See comprehensive documentation on the github repository
# https://github.com/jacobbeierle/JABS_nextflow_postprocess/tree/main

# ======================
# Set up / Libraries / Options
# ======================

library(tidyverse)
library(writexl)

# =======================
# Argument Configuration
# ======================

# Set your paths here:
input.dir <- "~/kumar-group/SING-grant/NextflowOutput/"
output.dir <- "~/kumar-group/SING-grant/Nextflow_postprocess"

# Set your parameters here:
params <- list(
  expected_length = 60*60*30 + 5*30,  
  max_tracklet_per_hour = 6,
  max_missing_pose = 0.005,
  max_missing_segmentation = 0.2,
  max_missing_keypoint = 0.01,
  fecal_boli_quantile_plotting = 0.05
)

# Optional: Override with YAML file (uncomment if you want to use this)
# yaml_vals <- yaml::read_yaml("path/to/your/config.yaml")
# params <- modifyList(params, yaml_vals)

# Print final configuration
cat("=== QC CONFIGURATION ===\n")
cat("Input directory: ", input.dir, "\n")
cat("Output directory:", output.dir, "\n")
cat("QC Parameters:\n")
print(params)

# ==================
# Create output directories
# ==================

for (subdirectory in c("final_nextflow_feature_data",
                       "qc/nextflow_qc_logs",
                       "qc/missing_or_dup_data",
                       "qc/qc_figs")) {
  dir.path = file.path(output.dir, subdirectory)
  dir.create(dir.path, recursive = T, showWarnings = F)
}

# ================
# Process and Publish QC logs with success or failure annotated in a CSV
# ================
# Read QC files in NextFlow_Output directory
qc_log <- list.files(
              path = "~/kumar-group/SING-grant/NextflowOutput/",
              pattern = "qc_batch_",
              full.names = TRUE,
              recursive = TRUE
              ) %>%
          read_csv(id = "QC_file", show_col_types = FALSE)

# Record why QC failed for each video
qc_log <- qc_log %>%
  mutate(passed_duration_QC = video_duration == params$expected_length,
         passed_tracklet_QC = pose_tracklets < params$max_tracklet_per_hour * params$expected_length / 108000,
         passed_segmentation_QC = seg_counts > (1 - params$max_missing_segmentation) * params$expected_length,
         passed_pose_QC = pose_counts > (1 - params$max_missing_pose) * params$expected_length,
         passed_kp_QC = missing_keypoint_frames < params$max_missing_keypoint * params$expected_length)

# Apply thresholds defined above to create a separate 'failed QC' data frame
qc_log.failed <- qc_log %>%
  filter(if_any(starts_with("passed_"), ~ !.x))

# Write final Nextflow QC files for review by a human
write.csv(qc_log, file.path(output.dir, "qc/nextflow_qc_logs/qc_all.csv"), row.names = FALSE)
write.csv(qc_log.failed,  file.path(output.dir, "qc/nextflow_qc_logs/qc_failed.csv"), row.names = FALSE)

# List of expected videos from the QC log files
expected_videos <- qc_log$video_name |>
  gsub("_with_fecal_boli", "", x = _) |>
  gsub("_filtered", "", x = _) |>
  sub("^/",          "", x = _) |>
  unique()    # Remove duplicate

# ==================
# Helper functions for processing output data
# ==================
read_raw_data <- function(input_dir, pattern) {
  ##########################################################################
  # Read in multiple CSV files from a directory and harmonize the ID column
  #
  # Args:
  #   input_dir: path to the folder containing CSV files
  #   pattern: regex pattern to match file names
  #
  # Returns:
  #   A tibble with:
  #     - NetworkFilename as the first column
  #     - Cleaned NetworkFilename (no "_corrected", "_filtered", ".avi", leading "." or "/")
  #
  # Notes:
  #   - If NetworkFilename does not exist, the first column is used as ID
  #   - Useful for harmonizing output from different workflows before merging
  ##########################################################################
  
  # Read in data from multiple csv files of the same patterns
  raw_data <- list.files(path = input_dir, pattern = pattern,
                         recursive = T, full.names = T) |>
    read_csv(show_col_types = FALSE)
  
  # Use NetworkFilename as the ID column (move it to the first column if not already so)
  if ("NetworkFilename" %in% names(raw_data)) {
    raw_data <- relocate(raw_data, NetworkFilename, .before = 1)
  } else {
    colnames(raw_data)[1] <- "NetworkFilename"
  }
  
  # Clean and harmonize the NetworkFilename across data
  raw_data$NetworkFilename <- raw_data$NetworkFilename |>
    gsub("_corrected", "", x = _) |>
    gsub("_filtered",  "", x = _) |>
    gsub("\\.avi$",    "", x = _) |>
    sub("^\\.",        "", x = _) |>
    sub("^/",          "", x = _)
  
  return(raw_data)
}

check_missing_and_dup <- function(expected_videos, data_df, corr_thres = 0.99) {
  ##########################################################################
  # Check for missing videos and duplicated rows in a dataset
  #
  # Args:
  #   expected_videos: character vector of expected NetworkFilename values
  #   data_df: dataframe containing a NetworkFilename column and data
  #   corr_thres: a number for how much correlated 2 rows are to be flagged
  #
  # Returns:
  #   A list with three elements:
  #     - missing_qc: videos present in expected_videos but missing in output_file
  #     - missing_output: videos present in output_file but missing in expected_videos
  #     - dup_data: rows in output_file that are duplicated (ignoring NetworkFilename),
  #                 after rounding numeric columns to zero digits
  #
  # Notes:
  #   - Useful for QC of experimental datasets, e.g., fecal boli or gait data
  #   - Rounds numeric columns before checking for duplicates to account for minor differences
  ##########################################################################
  
  video_missing_output <- setdiff(expected_videos, data_df$NetworkFilename)
  video_missing_in_qc <- setdiff(data_df$NetworkFilename, expected_videos)
  
  # Check for rows with identical data (i.e. something went wrong in video recording)
  dup_idx <- which(duplicated(data_df[, -1]) | duplicated(data_df[, -1], fromLast = TRUE))
  
  # Check for rows with high correlation 
  corr_mat <- data_df %>% 
    dplyr::select(where(is.numeric)) %>%
    scale() %>% t() %>%
    cor(use = "pairwise.complete.obs")
  
  # correlated row pairs above threshold
  row_pairs <- which(corr_mat > corr_thres & row(corr_mat) != col(corr_mat), arr.ind = TRUE)
  cor_idx <- unique(c(row_pairs[,1], row_pairs[,2]))
  
  # UNION: rows that are duplicated OR highly correlated
  all_idx <- sort(unique(c(dup_idx, cor_idx)))
  
  duplicated_rows <- data_df[all_idx,] 
  
  return(list(missing_qc = video_missing_in_qc,
              missing_output = video_missing_output,
              dup_data = duplicated_rows))
}

# ==================
# Process fecal boli data
# ==================
# Concatenate all instances
fecal_boli.raw <- read_raw_data(input_dir = "~/kumar-group/SING-grant/NextflowOutput/",
                                pattern = "fecal_boli.csv")

# Check for missing and duplicated data
fecal_boli.summary <- check_missing_and_dup(expected_videos = expected_videos,
                                            data_df = fecal_boli.raw)
# Write out all raw, merged fecal boli counts
write.csv(fecal_boli.raw, file.path(output.dir, "final_nextflow_feature_data/fecal_boli_raw.csv"), row.names = FALSE)

# ================
# Fecal boli QC plots
# ================
# Pivot longer to facilitate plotting for QC
fecal_boli.plot <- fecal_boli.raw |> 
  pivot_longer(
    cols = !c(NetworkFilename, nextflow_version),
    names_to = "min", 
    values_to = "fecal_boli",
    values_drop_na = TRUE) |> 
  mutate(min = parse_number(min))

# Plot fecal boli QC measures 
outFileNamePDF <- file.path(output.dir, "qc/qc_figs/fecal_boli_qc_figs.pdf") 
pdf(outFileNamePDF, 6, 6)

# Growth curve for all mice
ggplot(fecal_boli.plot, aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
  geom_line() +
  labs(title = "Fecal boli growth, all mice") +
  theme(legend.position = "none")

# Plot mice with lowest fecal boli
fecal_boli.plot  |> 
  summarise(across(fecal_boli, max), .by = NetworkFilename) |> 
  slice_min(fecal_boli, prop = params$fecal_boli_quantile_plotting) |> 
  select(NetworkFilename) |> 
  merge(fecal_boli.plot, by.x = "NetworkFilename") |> 
  ggplot(aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
    geom_line() +
    labs(title = paste("Lowest ", params$fecal_boli_quantile_plotting*100, "% of fecal boli mice", sep = "")) +
    theme(legend.position = "none")

# Plot mice with highest fecal boli
fecal_boli.plot  |> 
  summarise(across(fecal_boli, max), .by = NetworkFilename) |> 
  slice_max(fecal_boli, prop = params$fecal_boli_quantile_plotting) |> 
  select(NetworkFilename) |> 
  merge(fecal_boli.plot, by.x = "NetworkFilename") |> 
  ggplot(aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
    geom_line() +
    labs(title = paste("Highest ", params$fecal_boli_quantile_plotting*100, "% of fecal boli mice", sep = "")) +
    theme(legend.position = "none")

# Histogram of final fecal boli count
fecal_boli.plot  |> 
  arrange(desc(min)) |> 
  distinct(NetworkFilename, .keep_all = TRUE) |> 
  ggplot(aes(fecal_boli))+
    geom_histogram(binwidth = 1, boundary = 0)+
    labs(title = "Fecal boli highest bin, all mice") +
    ylab("count")

invisible(dev.off())

# ===============
# Process Gait Data
# ===============
# Import Gait Data
gait.raw <- read_raw_data(input_dir = "~/kumar-group/SING-grant/NextflowOutput/",
                          pattern = "gait.csv")

video_level_metrics = c("Distance Traveled", "Body Length", "Speed", "Speed Variance", "nextflow_version")

gait.wide_format <- gait.raw %>%
  # Remove variance measures from speed bins with fewer than 3 strides
  mutate(across(.cols = contains("Variance") & !all_of(video_level_metrics),
                .fns = ~ ifelse(`Stride Count` < 3, NA, .x))) %>%
  # Convert to wide format
  pivot_wider(id_cols = c(NetworkFilename, all_of(video_level_metrics)),
              names_from = `Speed Bin`, 
              values_from = -c(NetworkFilename, all_of(video_level_metrics), `Speed Bin`),
              names_sep = ".") %>%
  mutate(across(.cols = c(`Stride Count.10`, `Stride Count.15`, `Stride Count.20`, `Stride Count.25`),
                .fns  = ~ replace_na(.x, 0)))

# Check for missing and duplicated data with gait
gait.summary <- check_missing_and_dup(expected_videos = expected_videos, 
                                      data_df = gait.wide_format)

#output to final CSV
write.csv(gait.wide_format, file.path(output.dir, "final_nextflow_feature_data/gait_final.csv"), row.names = FALSE)

# =================
# Process JABS Feature Data
# =================
JABS.features <- read_raw_data(input_dir = "~/kumar-group/SING-grant/NextflowOutput/",
                               pattern = "features.csv")

# Check for missing data in JABS.features
JABS.features.summary <- check_missing_and_dup(expected_videos = expected_videos,
                                               data_df = JABS.features)

# Write the final csv
write.csv(JABS.features, file.path(output.dir, "final_nextflow_feature_data/JABS_features_final.csv"), row.names = FALSE)

# ====================
# Process morphometrics feature data
# ====================
morpho.raw <- read_raw_data(input_dir = "~/kumar-group/SING-grant/NextflowOutput/",
                            pattern = "morphometrics.csv")

# Check for missing and duplicated data in morphometric outputs
morpho.summary <- check_missing_and_dup(expected_videos = expected_videos,
                                        data_df = morpho.raw)

# =================
# Report and output data for all warnings
# =================
# Videos in QC but not in output
all_missing_data <- list("fecal_boli" = fecal_boli.summary$missing_output,
                         "gait" = gait.summary$missing_output,
                         "JABS_features" = JABS.features.summary$missing_output,
                         "Morphometrics" = morpho.summary$missing_output) 

# Videos in output but not in QC
videos_not_in_qc_report <- list("fecal_boli" = fecal_boli.summary$missing_qc,
                                "gait" = gait.summary$missing_qc,
                                "JABS_features" = JABS.features.summary$missing_qc,
                                "morphometrics" = morpho.summary$missing_qc)

# Output data that is duplicated in the data frames
all_duplicated_data <- list("fecal_boli" = fecal_boli.summary$dup_data,
                            "gait" = gait.summary$dup_data,
                            "JABS_features" = JABS.features.summary$dup_data,
                            "morphometrics" = morpho.summary$dup_data)

no_missing_output <- all(sapply(all_missing_data, length) == 0)
no_missing_qc     <- all(sapply(videos_not_in_qc_report, length) == 0)
no_dups           <- all(sapply(all_duplicated_data, nrow) == 0)

# Summarize and write warnings
cat("=== FINAL ERROR REPORT ===\n")
if (no_missing_output && no_missing_qc && no_dups) {
  cat("No errors to report\n")
  # Create placeholder files
  write.csv("No data missing", file.path(output.dir, "qc/missing_or_dup_data/missing_data.csv"), row.names = FALSE)
  write.csv("No data missing", file.path(output.dir, "qc/missing_or_dup_data/videos_not_in_qc_report.csv"), row.names = FALSE)
  write_xlsx(as.data.frame("No duplicated data"), file.path(output.dir, "qc/missing_or_dup_data/duplicated_data.xlsx"), row.names = FALSE)
} else {

  # Check for missing output data
  if (no_missing_output) {
    write.csv("No data missing", file.path(output.dir, "qc/missing_or_dup_data/missing_data.csv"), row.names = FALSE)
  } else {
    all_missing_data <- all_missing_data %>% 
      enframe(., name = "outputType", value = "video_path") %>% 
      unnest(cols = video_path) %>% 
      mutate(missing = TRUE) %>% 
      pivot_wider(id_cols = video_path, names_from = outputType, values_from = missing)
    write.csv(all_missing_data, file.path(output.dir, "qc/missing_or_dup_data/missing_data.csv"), row.names = FALSE)    
    cat(paste("Missing", colnames(all_missing_data)[-1], "data"), sep = "\n")
  }
  
  # Check for missing QC data
  if (no_missing_qc) {
    write.csv("No data missing", file.path(output.dir, "qc/missing_or_dup_data/videos_not_in_qc_report.csv"), row.names = FALSE)
  } else {
    videos_not_in_qc_report <- videos_not_in_qc_report %>%
      enframe(., name = "outputType", value = "video_path") %>% 
      unnest(cols = video_path) %>% 
      mutate(missing = TRUE) %>% 
      pivot_wider(id_cols = video_path, names_from = outputType, values_from = missing)
    write.csv(videos_not_in_qc_report, file.path(output.dir, "qc/missing_or_dup_data/videos_not_in_qc_report.csv"), row.names = FALSE)    
    cat(paste("Missing video in QC for", colnames(videos_not_in_qc_report)[-1], "data"), sep = "\n")
  }
  
  # Check for duplicated data
  if (no_dups) {
    write_xlsx(as.data.frame("No duplicated data"), file.path(output.dir, "qc/missing_or_dup_data/duplicated_data.xlsx"), row.names = FALSE)
  } else {
    write_xlsx(all_duplicated_data, path = file.path(output.dir, "qc/missing_or_dup_data/duplicated_data.xlsx"))
    cat(paste("Duplicated data for", names(all_duplicated_data)[sapply(all_duplicated_data, nrow) != 0]), sep = "\n")
  }
}
