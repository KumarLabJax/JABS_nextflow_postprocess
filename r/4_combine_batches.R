#A script to merge NextFlow outputs for final analysis
#Developed by Dr. Jake Beierle (don't forget the Dr., it's important)


#----Documentation----
#See comprehensive documentation on the github repository
#https://github.com/jacobbeierle/JABS_nextflow_postprocess/tree/main

#Load Libraries----------------------------------------------------------------
library(tidyverse)
library(data.table)
library(janitor)
options(error = NULL) #helps with error handling in functions checking for directories and filenames
#########################################################Define the working directory and variables------------------------------------------------------------------

##If you are not using an R project, set your working directory

#working.directory <- "C:\\Users\\beierj\\Desktop\\2025-04-09_NTG_C1-C5_Analysis"

#Create functions--------------------------------------------------------------

#A function to compares two vectors for the presence or absence of NetworkFilenames
compare_NetworkFilenames <- function(x, y){
  #Check if there are missing NetworkFilenames in the vectors x and y
  if(setequal(x, y) == FALSE){
    #Create data frames with missing data, create boolean col for failure reason
    x.notin.y <- tibble(
      NetworkFilename = setdiff(x, y),
      !!gsub("\\$NetworkFilename", "", paste0(deparse(substitute(x)), "_not_in_", deparse(substitute(y)))) := 1
    )
    
    y.notin.x <- tibble(
      NetworkFilename = setdiff(y, x),
      !!gsub("\\$NetworkFilename", "", paste0(deparse(substitute(y)), "_not_in_", deparse(substitute(x)))) := 1
    )
    
    #Full joining two dataframes containing missing values by a dummy value N
    out <- full_join(x.notin.y, y.notin.x)
    out[is.na(out)] <- 0
    out
  }else{NULL} #If nothing is missing, return a NULL vector
}

#A function to make sure the file path exists, and if not, stop the code with an error being output
check_files_exist <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste0(
      "YOU DO NOT HAVE A '", file_path, "' FILE\n",
      "YOU NEED A '", file_path,"' FILE TO CONTINUE"
    ), call. = FALSE)
  }
}


#Set working directory, and check for expected file names and directories---------------------------------

#Set the expected directory to final nextflow feature csv files
final.nextflow.feature.dir <- "Nextflow_Output/final_nextflow_feature_data"
#Set directories for QC directories
qc.missing_dup.dir <- file.path("qc", "missing_or_dup_data")
qc.figs.dir <- file.path("qc", "qc_figs")

#Set the expected metadata filename
metadata.filename <- "metadata.csv"
#Set the expected name of the index of videos to exclude
excluded.videos <- "videos_to_exclude.txt"

#If you have defined a working directory above, set it here
if(exists("working.directory")){
  setwd(working.directory)
}

#Check for 'Nextflow_Output/final_nextflow_feature_data/'
#This is where processed data from previous script should be published, report error if there is no directory

if(!dir.exists(final.nextflow.feature.dir)){
  stop("YOU DO NOT HAVE A 'Nextflow_Output/final_nextflow_feature_data' DIRECTORY IN YOUR WORKING DIRECTORY \n
       YOU NEED TO RUN 'NextFlow_Output_QC_Postprocess_1.R' TO PRODUCE THESE FINAL FILES") 
}else{
  #Check to ensure all files needed for this code are present
  check_files_exist(file.path(final.nextflow.feature.dir, "fecal_boli_final.csv"))
  check_files_exist(file.path(final.nextflow.feature.dir, "gait_final.csv"))
  check_files_exist(file.path(final.nextflow.feature.dir, "JABS_features_final.csv"))
  check_files_exist(file.path(final.nextflow.feature.dir, "morphometrics_final.csv"))
  }

#Check that qc directory for missing and duplicated data exists, report error if it does not
if(!dir.exists(qc.missing_dup.dir)){
  stop("YOU DO NOT HAVE A 'qc/missing_or_dup_data' DIRECTORY IN YOUR WORKING DIRECTORY\n
       YOU SHOULD HAVE THEM IF YOU RAN 'NextFlow_Output_QC_Postprocess_1.R' TO PRODUCE FINAL NEXTFLOW FEATURE FILES") 
}else{qc.missing_dup.dir <- file.path("qc", "missing_or_dup_data")}


#Read in the final nextflow features files from 'NextFlow_Output_QC_Postprocess_1.R'------

#Read in gait
gait.final <- list.files(
  path = final.nextflow.feature.dir,
  pattern = "gait_final",
  full.names = TRUE) |> 
  read_csv()

#Read in morphometrics
morpho.final <- list.files(
  path = final.nextflow.feature.dir,
  pattern = "morphometrics_final",
  full.names = TRUE) |> 
  read_csv() |> 
  relocate(NetworkFilename)

#Read in fecal boli
fecal_boli.final <- list.files(
  path = final.nextflow.feature.dir,
  pattern = "fecal_boli_final",
  full.names = TRUE) |> 
  read_csv() |> 
  relocate(NetworkFilename)

#Read in JABS features
JABS.final <- list.files(
  path = final.nextflow.feature.dir,
  pattern = "features_final",
  full.names = TRUE) |> 
  read_csv()

#Check NetworkFilenames, and merge it all together!---------------------------------------------------------

#Create a list of dataframes representing NetworkFilenames missing in pairwise comparison of dataframes
NetworkFilename.check <- list(
  compare_NetworkFilenames(gait.final$NetworkFilename, morpho.final$NetworkFilename),
  compare_NetworkFilenames(gait.final$NetworkFilename, fecal_boli.final$NetworkFilename),
  compare_NetworkFilenames(gait.final$NetworkFilename, JABS.final$NetworkFilename),
  compare_NetworkFilenames(morpho.final$NetworkFilename, fecal_boli.final$NetworkFilename),
  compare_NetworkFilenames(morpho.final$NetworkFilename, JABS.final$NetworkFilename),
  compare_NetworkFilenames(fecal_boli.final$NetworkFilename, JABS.final$NetworkFilename)
)
#Remove the empty elements, i.e. instances where no NetworkFilenames were missing
NetworkFilename.check <- Filter(Negate(is.null), NetworkFilename.check)
#Merge all items on the list into a single data frame
merged.NetworkFilename.check <- Reduce(function(x, y) merge(x, y, by = "NetworkFilename", all = TRUE), NetworkFilename.check)


#Publish error reporting data frame if the length is above 0
#Otherwise publish confirmation of no errors and merge all dataframes
if(length(merged.NetworkFilename.check)){
  #Publish error report
  write.csv(merged.NetworkFilename.check, file.path(qc.missing_dup.dir, "NetworkFilenames_missing_in_data.csv"), row.names = FALSE)
}else{
  #Report that no data is missing
  merged.NetworkFilename.check <- "NETWORKFILE NAMES MATCH PERFECTLY ACROSS DATAFRAMES"
  #Publish error free report
  write.csv(merged.NetworkFilename.check, file.path(qc.missing_dup.dir, "NetworkFilenames_missing_in_data.csv"), row.names = FALSE)
}

#Merge all dataframes
nextflow_dataset.metadata_not_merged <- gait.final |> 
  merge(JABS.final, by = "NetworkFilename") |> 
  merge(morpho.final, by = "NetworkFilename") |> 
  merge(fecal_boli.final, by = "NetworkFilename")


#Remove Videos that Failed QC, Requires manual screening------------------------

if(!file.exists(excluded.videos)){
  videos_to_exclude.exists <- FALSE
}else{
  videos_to_exclude <- read.table(excluded.videos, quote="\"", comment.char="")
  if(length(videos_to_exclude)){
    nextflow_dataset.metadata_not_merged <- nextflow_dataset.metadata_not_merged[!nextflow_dataset.metadata_not_merged$NetworkFilename %in% videos_to_exclude$V1, ]
    videos_to_exclude.exists <- TRUE
  }else(videos_to_exclude.exists <- FALSE)
}

#Merge with metadata in top directory of project folder-------------------------

#Read Metadata, or throw error if it does not exist
if(file.exists(metadata.filename)){
  metadata <- read_csv(metadata.filename,
                       col_types = cols(MouseID = col_character()) )
}else{
  stop(paste0("YOU DO NOT HAVE A '", metadata.filename, "' FILE IN YOUR WORKING DIRECTORY YOU NEED A '", metadata.filename, "' file TO RUN THIS SCRIPT")) 
}

#Prepare NetworkFilename to Filename, to allow merging of metadata
nextflow_dataset.metadata_not_merged$FileName <- str_split_i(nextflow_dataset.metadata_not_merged$NetworkFilename, "/", i=4)
nextflow_dataset.metadata_not_merged$FileName <- gsub("_trimmed.avi", "", nextflow_dataset.metadata_not_merged$FileName)
nextflow_dataset.metadata_not_merged <- relocate(nextflow_dataset.metadata_not_merged, FileName)

nextflow_dataset.metadata_not_merged$MouseID <- str_split_i(nextflow_dataset.metadata_not_merged$FileName, "_", i=1)

#Check that all mice in metadata are represented in dataset, and vice versa
metadata.qc.check <- compare_NetworkFilenames(nextflow_dataset.metadata_not_merged$MouseID, metadata$MouseID)

#Publish error report
if(length(metadata.qc.check)){
  #Publish error report
  write.csv(metadata.qc.check, file.path(qc.missing_dup.dir, "mice_missing_in_metadata.csv"), row.names = FALSE)
}else{
  #Publish error free report
  metadata.qc.check <- "NO MICE MISSING IN METADATA & VICE VERSA"
  write.csv(metadata.qc.check, file.path(qc.missing_dup.dir, "mice_missing_in_metadata.csv"), row.names = FALSE)
}

#Merge Metadata
nextflow_dataset.final <- merge(metadata, nextflow_dataset.metadata_not_merged, by = "MouseID")

#Publish CSV--------------------------------------------------------------------

write_csv(nextflow_dataset.final, file.path(final.nextflow.feature.dir, "merged_nextflow_dataset.csv"))

#Error reporting----------------------------------------------------------------
error.reporting <- NULL

#Report to terminal if data is missing from QC log but in feature tables'
#Because passed QC assigns this object to a string reporting the lack of failed QC
#We can use is.character() to determine if QC passed
if(!is.character(merged.NetworkFilename.check)){
  error.reporting <- c(error.reporting, "NETWORKFILE NAMES DO NOT MATCH PERFECTLY ACROSS FINAL NEXTFLOW DATAFRAMES")
}

#Report to terminal if mice are missing from metadata and vice versa
if(!is.character(metadata.qc.check)){
  error.reporting <- c(error.reporting, "MICE IN METADATA AND DATA DO NOT MATCH PERFECTLY")
}


#Check that videos_to_exclude.txt is present and has data in it
if(videos_to_exclude.exists == FALSE){
  error.reporting <- c(error.reporting, "THERE IS NO 'videos_to_exclude.txt' IN YOUR WORKING DIRECTORY OR THERE ARE NO VIDEOS LISTED WITHIN IT, ARE THERE REALLY NO VIDEOS TO EXCLUDE BASED ON MANUALLY SCREEDED QC?")
}


#Print out all errors after code done running
if(length(error.reporting) == 0){
  print("FINAL ERROR REPORT: NO ERRORS TO REPORT")
}else{
  print("FINAL ERROR REPORT:", )
  paste(error.reporting)
}
