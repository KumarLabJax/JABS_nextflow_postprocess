#!/usr/bin/env Rscript
# A script to clean NextFlow outputs for final processing
# Developed by Dr. Jake Beierle (don't forget the Dr., it's important)

# ----Documentation----
# See comprehensive documentation on the github repository
# https://github.com/jacobbeierle/JABS_nextflow_postprocess/tree/main

library(optparse)

# ======================
# Set up / Libraries / Options
# ======================
library(tidyverse)

# =======================
# Set QC and other Values
# ======================

# Set the QC values you will use to screen the Nextflow QC files
expected.length<- 60*60*30 + 5*30     # Video clipping duration plus 5 seconds)
max.tracklets.per.hour <- 6
# Max frames (as percent of all frames) missing pose
max.percent.segmentation.missing <- 0.2
# Max percentage of Frames missing pose
max.percent.pose.missing <- 0.005
# Max percentage of Frames missing pose
max.percent.kp.missing <- 0.01
# Proportion of Highest and lowest fecal boli mice to plot separately for QC
fecal_boli_percent_threshold <- 0.05

# ==================
# Create output directories
# ==================
# Set output directory
output.dir <- "/projects/kumar-lab/USERS/nguyetu/SING-grant/Nextflow_postprocess"
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
          read_csv(id = "QC_file")

# Record why QC failed for each video
qc_log <- qc_log %>%
  mutate(passed_duration_QC = video_duration == expected.length,
         passed_tracklet_QC = pose_tracklets < max.tracklets.per.hour * expected.length / 108000,
         passed_segmentation_QC = seg_counts > (1 - max.percent.segmentation.missing) * expected.length,
         passed_pose_QC = pose_counts > (1 - max.percent.kp.missing) * expected.length,
         passed_kp_QC = missing_keypoint_frames < max.percent.kp.missing * expected.length)

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

check_missing_and_dup <- function(expected_videos, data_df) {
  ##########################################################################
  # Check for missing videos and duplicated rows in a dataset
  #
  # Args:
  #   expected_videos: character vector of expected NetworkFilename values
  #   data_df: dataframe containing a NetworkFilename column and data
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
  
  video_missing_in_qc <- setdiff(expected_videos, data_df$NetworkFilename)
  video_missing_output <- setdiff(data_df$NetworkFilename, expected_videos)
  
  # Check for rows with identical data (i.e. something went wrong in video recording)
  # I remove the network file name col because data may have been mislabeled
  duplicated_rows <- data_df %>%
    mutate(across(where(is.numeric), ~ round(.x, digits = 0))) %>%
    filter(duplicated(across(-1)) | duplicated(across(-1), fromLast = T))

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
    cols = !NetworkFilename,
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

# Plot 10% mice with lowest fecal boli
fecal_boli.plot  |> 
  summarise(across(fecal_boli, max), .by = NetworkFilename) |> 
  slice_min(fecal_boli, prop = fecal_boli_percent_threshold) |> 
  select(NetworkFilename) |> 
  merge(fecal_boli.plot, by.x = "NetworkFilename") |> 
  ggplot(aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
    geom_line() +
    labs(title = paste("Lowest ", fecal_boli_percent_threshold*100, "% of fecal boli mice", sep = "")) +
    theme(legend.position = "none")

# Plot 10% mice with highest fecal boli
fecal_boli.plot  |> 
  summarise(across(fecal_boli, max), .by = NetworkFilename) |> 
  slice_max(fecal_boli, prop = fecal_boli_percent_threshold) |> 
  select(NetworkFilename) |> 
  merge(fecal_boli.plot, by.x = "NetworkFilename") |> 
  ggplot(aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
    geom_line() +
    labs(title = paste("Highest ", fecal_boli_percent_threshold*100, "% of fecal boli mice", sep = "")) +
    theme(legend.position = "none")

# Histogram of final fecal boli count
fecal_boli.plot  |> 
  arrange(desc(min)) |> 
  distinct(NetworkFilename, .keep_all = TRUE) |> 
  ggplot(aes(fecal_boli, ifelse(after_stat(count) > 0, after_stat(count), NA)))+
    geom_histogram(binwidth = 1, boundary = 0)+
    labs(title = "Fecal boli highest bin, all mice") +
    ylab("count")

dev.off()

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
#Process morphometrics feature data
# ====================
morpho.raw <- read_raw_data(input_dir = "~/kumar-group/SING-grant/NextflowOutput/",
                            pattern = "morphometrics.csv")

# Check for missing and duplicated data in morphometric outputs
morpho.summary <- check_missing_and_dup(expected_videos = expected_videos,
                                        data_df = morpho.raw)

# =================
# Merge and output all missing data
# =================
# Missing output data (but video is in QC)
all_missing_data <- list("fecal_boli" = fecal_boli.summary$missing_output,
                         "gait" = gait.summary$missing_output,
                         "JABS_features" = JABS.features.summary$missing_output,
                         "morphometrics" = morpho.summary$missing_output) %>%
  enframe(., name = "outputType", value = "video") %>% 
  unnest(cols = video)
# Select the dfs with actual missing data represented, dropping the rest
publish_missing_data <- NULL
for(i in seq_along(all_missing_data)){
  if(length(all_missing_data[[i]]) > 0){
    if(length(publish_missing_data) == 0){
      publish_missing_data <- all_missing_data[[i]]
    } else {publish_missing_data <- full_join(publish_missing_data, all_missing_data[[i]])}
  }
}

# Fill the csv with something if all QC passes
if(length(publish_missing_data) == 0){
  publish_missing_data <- "NO DATA MISSING"
}

# Write the csv
write.csv(publish_missing_data, file.path(output.dir, "qc/missing_or_dup_data/missing_data.csv"), row.names = FALSE)

# Publish all vids in some feature csv, but not in the qc dataframe
# Combine into a single list
videos_not_in_qc_report <- list(
  "fecal_boli_videos_missing_in_qc" = fecal_boli_videos_missing_in_qc,
  "gait_videos_missing_in_qc" = gait_videos_missing_in_qc,
  "JABS_features_videos_missing_in_qc" = JABS_features_videos_missing_in_qc,
  "morphometrics_videos_missing_in_qc" = morphometrics_videos_missing_in_qc
)

# Select the dfs with actual missing data represented, dropping the rest
publish_videos_not_in_qc_report <- NULL
for(i in seq_along(videos_not_in_qc_report)){
  if(length(videos_not_in_qc_report[[i]]) > 0){
    if(length(publish_videos_not_in_qc_report) == 0){
      publish_missing_data <- videos_not_in_qc_report[[i]]
    }else{publish_missing_data <- full_join(publish_missing_data, videos_not_in_qc_report[[i]])}
  }else{}
}

# Fill the csv with something if all QC passes
if(length(publish_videos_not_in_qc_report) == 0){
  publish_videos_not_in_qc_report <- "NO DATA MISSING"
}

write.csv(publish_videos_not_in_qc_report, file.path(output.dir, "qc/missing_or_dup_data/videos_not_in_qc_report.csv"), row.names = FALSE)

# Output data that is duplicated in the data frames
duplicated_data <- list(
  "dup_gait" = duplicate_gait_rows,
  "dup_JABS" = duplicate_JABS_feature_rows,
  "dup_morpho" = duplicate_morphometrics_rows,
  "dup_fboli" = duplicate_fboli_rows
)

#If all objecst in the list have a length of 0 (i.e. empty), report no dupli
if(all(sapply(duplicated_data, function(x) nrow(x)==0))){
  duplicated_data <- "NO DUPLICATED DATA"
  write_xlsx(as.data.frame(duplicated_data), path = file.path(qc.missing_dup.dir, "duplicated_data.xlsx"))
}else{
  write_xlsx(duplicated_data, path = file.path(output.dir, "qc/missing_or_dup_data/duplicated_data.xlsx"))
}

# =========================
# Write warnings for failed QC
# =========================
error.reporting <- NULL

# Report to terminal if data is missing from QC log but in feature tables
if(!is.character(publish_videos_not_in_qc_report)){
  if(length(fecal_boli_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"FECAL BOLI DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
  if(length(gait_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"GAIT DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
  if(length(JABS_features_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"JABS FEATURE DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
  if(length(morphometrics_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"YOU ARE MISSING MORPHOMETRIC FEATURE DATA") }
}

# Report to terminal if data is missing feature tables but in QC log
if(!is.character(publish_missing_data)){
  if(length(videos_with_missing_fecal_boli)){error.reporting <- c(error.reporting,"YOU ARE MISSING FECAL BOLI DATA") }
  if(length(videos_with_JABS_features_missing)){ error.reporting <- c(error.reporting,"YOU ARE MISSING JABS FEATURE DATA") }
  if(length(videos_with_morphometrics_features_missing)){ error.reporting <- c(error.reporting,"MORPHOMETRIC DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
}

# Report to terminal if there is duplicated data
if(!is.character(duplicated_data)){ error.reporting <- c(error.reporting,"YOU HAVE DUPLICATED DATA!") }


#Print out all errors after code done running
if(length(error.reporting) == 0){
  print("FINAL ERROR REPORT: NO ERRORS TO REPORT")
}else{
  print("FINAL ERROR REPORT:", )
  paste(error.reporting, collapse = "\n")
}