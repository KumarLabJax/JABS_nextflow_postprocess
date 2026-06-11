#A script to manage merged NextFlow data for final analysis
#Developed by Dr. Jake Beierle (don't forget the Dr., it's important)
#10-26-2025

#Documentation----
#See comprehensive documentation on the github repository
#https://github.com/jacobbeierle/JABS_nextflow_postprocess/tree/main


#Load Libraries----
library(tidyverse)
library(janitor)
options(error = NULL) #helps with error handling in functions checking for directories and filenames

######################################## Define the working directory------------------------------------------------------------------

#If you are not using an R project, set your working directory
#working.directory <- "C:\\Users\\beierj\\Desktop\\2025-04-09_NTG_C1-C5_Analysis"

## Set some pithy project name to be appended to your curated dataset
#project_name <- "NTG_C1-5"


#Sometimes features from nextflow are not useful because they are not reliable in your experimental
#     context, are all equal to each other, etc. I remove these to make exploraroty data analysis and 
#     subsequent statistical analysis more reliable. Below are expmples
#These represent seizure features that have high false positives (Jake's opinion)
#Recording these features in a vector makes reporting what was removed easier at the end
#These strings will be used to filter out cols that contain them, so be specific in your choice of wording

#features.removed.manually <- c("Side_seizure", "Tail_jerk", "Wild_jumping")
#features.removed.manually <- NULL

#Set z-score threshold for outlier screening
zscore.threshold <- 6

#Set the number of phenotypes to print with the highest, and most common outliers
how.many.highest.zscore.outlier.plots <- 16
how.many.most.freq.outlier.plots <- 8


#Create functions--------------------------------------------------------------

#A function to make sure the file path exists, and if not, stop the code with an error being output
check_files_exist <- function(x) {
  if (!file.exists(x)) {
    stop(paste0(
      "YOU DO NOT HAVE A '", x, "' FILE\n",
      "YOU NEED A '", x,"' FILE TO CONTINUE"
    ), call. = FALSE)
  }
}
#add this line of code to help the error handling
options(error = NULL)

#A function that checks for directories
directory_check_creation <- function(x){
  x <- paste0(x)
  if(!dir.exists(x)){
    dir.create(x)
  }
  else{print(paste0("'/", substitute(x), "' already exists"))}
}

#This functions takes 2 numeric data cols, a ID col and a number, and plot a scatterplot
#with a linear best fit line, and highlights the n points orthoganally furthest from that line.
qc_plot_lm_outliers <- function(x, y, id, n.points = 5){
  temp.dataframe <- data.frame(
    x = x,
    y = y,
    id = id
  )
  
  # Fit the model: y ~ x
  temp.model <- lm(y ~ x, data = temp.dataframe)
  # Extract the slope and intercept
  slope <- coef(temp.model)[["x"]]
  intercept <- coef(temp.model)[["(Intercept)"]]
  
  # Calculate perpendicular distance from each point to the line of best fit
  temp.dataframe <- temp.dataframe |>
    mutate(
      perp_dist = abs(slope * x - y + intercept) / sqrt(slope^2 + 1)
    ) |>
    arrange(desc(perp_dist)) |>
    mutate(
      percentile_order = row_number(),
      plot_label = percentile_order <= n.points
    )
  
  #create labels
  lab.title <- paste0(
    gsub(".+\\$", "", deparse(substitute(x))),
    " vs. ",
    gsub(".+\\$", "", deparse(substitute(y)))
  )
  x.title <- gsub(".+\\$", "", deparse(substitute(x)))
  y.title <- gsub(".+\\$", "", deparse(substitute(y)))
  
  #create the plot
  ggplot(temp.dataframe, aes(x, y))+
    geom_point() +
    geom_point(data = temp.dataframe[temp.dataframe$plot_label == TRUE,], colour = "red3") +
    geom_smooth(method = "lm") +
    ggrepel::geom_text_repel(
      data = temp.dataframe[temp.dataframe$plot_label == TRUE,], # Use the filtered data frame
      aes(label = id),# Specify the column with car names for labels
      colour = "red3",
      box.padding = 0.5,
      point.padding = 0.5,
      min.segment.length = 0.2,
      force = 50,
      force_pull = 0,
      fontface = "bold"
    ) +
    labs(title = lab.title,
         x = x.title,
         y = y.title) +
    theme_minimal(base_size = 12, base_family = "") +
    theme(
      plot.background = element_rect(fill = "white", color = "black", size = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12 + 2, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "grey2", fill = NA)
    )
}


#Below are a family of functions for streamlined outlier figs and compilation

#Takes a dataframe with the first col being unique IDs and calculates the z-score for each col, 
#removing features with 0 variance, and retaining cols that return var(x) == NA
z_score_remove_novar <- function(df){
  #ID no var cols, excluding strings
  nonzero_variance_cols <- sapply(df, function(x) var(x, na.rm = TRUE)) == 0
  nonzero_variance_cols[is.na(nonzero_variance_cols)] <- FALSE
  #Remove no var cols
  df <- df[, !nonzero_variance_cols]
  
  #Calculate Z score of each video for each phenotype
  cbind(df[1], as.data.frame(apply(df[-1], 2, scale)))
  
}

#Takes a dataframe of z-scores who's first col is unique IDs, and returns all cols that
#have an abs(z_score) higher or lower than the set threshold
filter_by_zscore <- function(df, threshold, id_col = 1){
  #ID cols with z greater than 6
  temp <- df[-id_col]
  temp <- temp[, colSums(!is.na(temp)) >= 3]
  cols.with.sig.zscore <- apply(temp, 2, function(x) if(max(x, na.rm = TRUE) >= threshold |min(x, na.rm = TRUE) >= threshold){TRUE}else{FALSE})
  
  cbind(df[id_col], temp[cols.with.sig.zscore])
}

#Takes a df of z scores and values, with identical cols, lengthens them by id col
#Merges them
lengthen_and_merge <- function(z_score_df, value_df, id_col = 1){
  #Lengthen the dataframe to allow for easy faceting by phenotype, add interaction to merge
  #with the value dataframe
  id_col <- id_col
  z_score_df.long <- z_score_df |> 
    pivot_longer(
      cols = -all_of(id_col),
      names_to = c("phenotype"),
      values_to = "z_score",
    )|> 
    mutate(
      interaction = paste(!!sym(colnames(z_score_df)[id_col]), phenotype, sep = "-")
    )
  
  #Repeat the above with the value dataframe
  value_df.long <- value_df |> 
    pivot_longer(
      cols = -all_of(id_col),
      names_to = c("phenotype")
    )|> 
    mutate(
      interaction = paste(!!sym(colnames(value_df)[id_col]), phenotype, sep = "-")
    )
  
  #Remove duplicated data before merging
  value_df.long <- value_df.long |>
    select(interaction, value)
  
  #Merge by interaction of phenotype and unique identifier for video
  data.outlier.plot <- merge(z_score_df.long, value_df.long, by = "interaction")
  #Remove interaction col from output
  data.outlier.plot[colnames(data.outlier.plot)!="interaction"]
}

#A function that unifies all the indivitual functions above, and returns a result of
# lengthen and merge ready for plotting
unified_zscore_processing <- function(df, threshold){
  df.zscore <- z_score_remove_novar(df) |> 
    filter_by_zscore(threshold)
  if(ncol(df.zscore)>=2){
    df.values <- df[colnames(df.zscore)]
    lengthen_and_merge(df.zscore, df.values)
  }else{df.zscore <- NULL}
  
}

#Makes facet wrapped boxplots displaying z score of out outliers
#Expects input from lengthen_and_merge
plot_feature_outliers_zscore <- function(df, threshold, n_col = 4,id_col_num = 1){
  ggplot(df, aes(phenotype, z_score)) +
    geom_boxplot() + 
    geom_hline(aes(yintercept = threshold), color = "red") +
    
    geom_point(
      data = subset(df, abs(z_score) >= threshold),
      color = "red4",
      size = 2
    ) +
    ggrepel::geom_text_repel(
      data = subset(df, abs(z_score) >= threshold),
      aes(label = !!sym(colnames(df)[id_col_num])),
      size = 3,
      color = "red4",
      force = 50, 
      force_pull = 0
    ) +
    facet_wrap(vars(phenotype), ncol = n_col, scales = "free_x")
}

#Makes facet wrapped boxplots displaying values of out outliers
#Expects input from lengthen_and_merge
plot_feature_outliers_value <- function(df, threshold, n_col = 4, id_col_num = 1){
  ggplot(df, aes(phenotype, value)) +
    geom_boxplot() + 
    geom_point(
      data = subset(df, abs(z_score) >= threshold),
      color = "red4",
      size = 2
    ) +
    ggrepel::geom_text_repel(
      data = subset(df, abs(z_score) >= threshold),
      aes(label = !!sym(colnames(df)[id_col_num])),
      size = 3,
      color = "red4",
      force = 50, 
      force_pull = 0
    ) +
    facet_wrap(vars(phenotype), ncol = n_col, scales = "free")
}

#Pull out only outliers from data for report in csv
#Expects output from lengthen_and_merge
parse_outliers_from_all <- function(df, data_type_str, threshold){
  count_name <- paste0("count_", data_type_str)
  pheno_name <- paste0("phenotype_", data_type_str)
  
  df |> 
    filter(
      abs(z_score) >= threshold
    ) |> 
    select(!c(z_score, value)) |> 
    group_by(!!sym(colnames(df)[1])) |> 
    summarise(
      !!sym(count_name) := n(),
      !!sym(pheno_name) := str_c(phenotype, collapse = ", "))
}




#Set working directory, and check for expected file names and directories---------------------------------

#If you have defined a working directory above, set it here
if(exists("working.directory")){
  setwd(working.directory)
}

#Set the expected directory to final nextflow feature csv files
merged_nextflow_dataset.dir <- file.path("Nextflow_Output", "final_nextflow_feature_data")
qc.figs.dir <- file.path("qc", "qc_figs")


#Check that qc_figs directoryexists, report error if it does not
if(!dir.exists(qc.figs.dir)){
  stop("YOU DO NOT HAVE A 'qc/qc_figs' DIRECTORY IN YOUR WORKING DIRECTORY\n
       YOU SHOULD HAVE THEM IF YOU RAN 'NextFlow_Output_QC_Postprocess_1.R' TO PRODUCE FINAL NEXTFLOW FEATURE FILES") 
}else{
  directory_check_creation(file.path(qc.figs.dir, "zscore_boxplots"))
  zscore.boxplot.dir <- file.path(qc.figs.dir, "zscore_boxplots")
}


#Check for "merged_nextflow_dataset.csv" in its expected directory
#This is where processed data from previous script should be published
#Report error if there is no directory
check_files_exist(file.path(merged_nextflow_dataset.dir, "merged_nextflow_dataset.csv"))

#Load merged dataset and set variables for final data output----
data.Nextflow <- read_csv(file.path(merged_nextflow_dataset.dir, "merged_nextflow_dataset.csv"))


#Check that qc directory for missing and duplicated data exists, report error if it does not
if(!dir.exists(qc.figs.dir)){
  stop(
  paste0("YOU DO NOT HAVE A '", qc.figs.dir, "' DIRECTORY IN YOUR WORKING DIRECTORY
         YOU SHOULD HAVE THIS DIRECTORY IF YOU RAN 'NextFlow_Output_QC_Postprocess_1.R' TO PRODUCE FINAL NEXTFLOW FEATURE FILES") 
  )
}else{qc.missing_dup.dir <- file.path("qc", "missing_or_dup_data")}


#Remove features not useful for your data-------------

#Manually remove the features you defined at the top of this script
if(length(features.removed.manually)){
  data.Nextflow <- data.Nextflow |> 
    select(!contains(features.removed.manually))
}


#Some of the morphometric features are always identical, and I remove these now too.
# Identify columns where the variance is not zero
nonzero_variance_cols <- sapply(data.Nextflow, function(x) var(x, na.rm = TRUE)) == 0
#Replace the NAs with FALSE, as these represent strings etc
nonzero_variance_cols[is.na(nonzero_variance_cols)] <- FALSE

#Recording these features in a vector makes reporting what was removed easier at the end
features.removed.zerovar <- colnames(data.Nextflow)[nonzero_variance_cols]

# Subset the data frame to remove zerovar cols
data.Nextflow <- data.Nextflow[, !nonzero_variance_cols]


#Print a csv of the features removed with reasons
features.removed.all <- data.frame(
  feature_removed = c(features.removed.manually, features.removed.zerovar),
  reason = c( rep("manually_removed", length(features.removed.manually)),  rep("zero_variance", length(features.removed.zerovar)))
)
if(nrow(features.removed.all)){
  write.csv(features.removed.all, "features_removed_from_curated_dataset.csv", row.names = FALSE)
}else{
  features.removed.all <- "NO FEATURES REMOVED"
  write.csv(features.removed.all, "features_removed_from_curated_dataset.csv", row.names = FALSE)
}



######################################## Adjust your data frame to align with your analysis----------------------------

#Now  you need to adjust your data for you final analysis. Perhaps you have multiple days of testing, which
#       is encoded deep within the mess that is NetworkFilenames. Perhaps you have doses you wish to exclude
#       from the final analysis so you don't have to keep doing that over and over again. This will be the 
#       most variable part of your script.

#I also like to arrange my metadata in the order that makes sense to me at this point. For my own sake.
#I also take this moment to adjust the metadata names I do not like


##Provide string of day for analysis
#data.Nextflow$Day <- str_split_i(data.Nextflow$NetworkFilename, "/", i=3)
#data.Nextflow$Day <- gsub("D", "", data.Nextflow$Day)
#data.Nextflow$Day <- as.numeric(data.Nextflow$Day)

##I'm removing doses and days I don't really care about, i.e. 2.5 and 5 mg/kg
#data.Nextflow <- subset(data.Nextflow, data.Nextflow$Tx!= 5 & data.Nextflow$Tx != 2.5)
#data.Nextflow <- subset(data.Nextflow, data.Nextflow$Day != 21)

##Reorder the factors
#data.Nextflow <- relocate(data.Nextflow, NetworkFilename, FileName, MouseID, PenID, ExptNumber, Cohort, Sex = sex, Day, Tx, LL)

########################################Summary information reporting--------------------------------------------------
#You can also take this time to report useful information for your analysis

#Example generating summaries of number of samples/timepoint
#data.Nextflow |> 
#  group_by(Day, Tx) |> 
#  summarise(N = n()) |> 
#  arrange(desc(Tx), Day) |> 
#  write_csv("final_n_per_timepoint.csv")


#Publish the final, curated data set--------------------------------------------
write.csv(data.Nextflow, paste0(project_name, "_final_nextflow_dataset.csv"), row.names = FALSE)

#Before QC plotting, I use janitor to clean col names because Nextflow col names are rough to work with----
colnames(data.Nextflow) <- colnames(clean_names(data.Nextflow))

#Preliminary QC figures: outlines from linear model of 2 phenotypes-------------

#Publish figures of linear relationships between 2 variables, with the 5 most distant points labeled
#see doc in function section for details
pdf(file.path(qc.figs.dir, "scatter_plot_lm_figs.pdf"), 7,7)
print(
  qc_plot_lm_outliers(data.Nextflow$distance_traveled,
                      data.Nextflow$bin_avg_55_locomotion_distance_cm,
                      data.Nextflow$network_filename)
)

print(
  qc_plot_lm_outliers(data.Nextflow$bin_sum_55_jumping_bout_behavior,
                      data.Nextflow$bin_sum_55_escape_bout_behavior,
                      data.Nextflow$network_filename)
)

print(qc_plot_lm_outliers(data.Nextflow$bin_sum_55_freeze_bout_behavior,
                          data.Nextflow$bin_sum_55_freezing_bout_behavior,
                          data.Nextflow$network_filename)
)

print(qc_plot_lm_outliers(data.Nextflow$distance_traveled,
                          data.Nextflow$speed,
                          data.Nextflow$network_filename)
)

print(qc_plot_lm_outliers(data.Nextflow$bin_sum_55_in_periphery_time_secs,
                          data.Nextflow$bin_sum_55_in_corner_time_secs,
                          data.Nextflow$network_filename)
)
dev.off()


#Preliminary QC figures: checking for outliers using z-score-----
#Subset datasets and create Network_filename + phenotypes

#Read in gait colnames for subsetting outlier data plots
gait.cols <- list.files(
  path = merged_nextflow_dataset.dir,
  pattern = "gait_final",
  full.names = TRUE) |> 
  read_csv(n_max = 1) |> 
  clean_names() |> 
  colnames()

#Read in gait colnames for subsetting outlier data plots
morpho.cols <- list.files(
  path = merged_nextflow_dataset.dir,
  pattern = "morphometrics_final",
  full.names = TRUE) |> 
  read_csv(n_max = 1) |> 
  clean_names()|> 
  colnames()

#Read in gait colnames for subsetting outlier data plots
JABS.cols <- list.files(
  path = merged_nextflow_dataset.dir,
  pattern = "features_final",
  full.names = TRUE) |> 
  read_csv(n_max = 1) |> 
  clean_names()|> 
  colnames()

#Gait outlier screening
gait.outliers <- unified_zscore_processing(data.Nextflow[colnames(data.Nextflow) %in% gait.cols], zscore.threshold)

if(length(gait.outliers) != 0){
  pdf(file.path(zscore.boxplot.dir, "gait_boxplot_outlier_figs.pdf"), 16, 20)
  print(plot_feature_outliers_zscore(gait.outliers, zscore.threshold))
  print(plot_feature_outliers_value(gait.outliers, zscore.threshold))
  dev.off()
}

#Morphometric outlier screening
morpho.outliers <- unified_zscore_processing(data.Nextflow[colnames(data.Nextflow) %in% morpho.cols], zscore.threshold)

if(length(morpho.outliers) != 0){
  pdf(file.path(zscore.boxplot.dir, "morpho_boxplot_outlier_figs.pdf"), 16, 20)
  print(plot_feature_outliers_zscore(morpho.outliers, zscore.threshold, 3))
  print(plot_feature_outliers_value(morpho.outliers, zscore.threshold, 3))
  dev.off()
}

#JABSmetric outlier screening
#I should really subset these features
JABS.outliers <- unified_zscore_processing(data.Nextflow[colnames(data.Nextflow) %in% JABS.cols], zscore.threshold)

unique(JABS.outliers$phenotype)


#Subset JABS figures if there are too many cols
if(length(JABS.outliers) != 0){
  if(length(unique(JABS.outliers$phenotype)) >= 20){
    
    x <- length(unique(JABS.outliers$phenotype))
    y <- ceiling(x/20)
    z <- c(1, 0 + ( 20 * 1:(y-1) )+1, x+1)
    
  pdf(file.path(zscore.boxplot.dir, "JABS_boxplot_outlier_figs.pdf"), 20, 30)
  for(i in 1:(length(z)-1)){
    J1 <- JABS.outliers |> 
      filter(
        phenotype %in% unique(JABS.outliers$phenotype)[c(z[i] : z[i+1]-1)]
      )
    print(plot_feature_outliers_zscore(J1, zscore.threshold))
    print(plot_feature_outliers_value(J1, zscore.threshold))
    
  }
  dev.off()
  
  }else{
    pdf(file.path(zscore.boxplot.dir, "JABS_boxplot_outlier_figs.pdf"), 20, 40)
    print(plot_feature_outliers(JABS.outliers, 5))
    dev.off()
  }
}

#Create report compiling number of outlier phenotypes for each mouse, for each dataframe--------------

gait.outlier.mice <- parse_outliers_from_all(gait.outliers, "gait", zscore.threshold)
morpho.outlier.mice <- parse_outliers_from_all(morpho.outliers, "morpho", zscore.threshold)
JABS.outlier.mice <- parse_outliers_from_all(JABS.outliers, "JABS", zscore.threshold)

all.outlier.mice.summary <- gait.outlier.mice |> 
  full_join(morpho.outlier.mice) |> 
  full_join(JABS.outlier.mice) |>
  group_by(!!sym(colnames(gait.outlier.mice)[1])) |> 
  mutate(
    count_gait = replace_na(count_gait, 0),
    count_morpho = replace_na(count_morpho, 0),
    count_JABS = replace_na(count_JABS, 0),
    count_all = sum(count_gait, count_morpho, count_JABS)
  ) |> 
  relocate(!!sym(colnames(gait.outlier.mice)[1]), count_all, count_gait, count_morpho, count_JABS) |> 
  arrange(desc(count_all))

write.csv(all.outlier.mice.summary, file.path(zscore.boxplot.dir, "outlier_videos_summary.csv"), row.names = FALSE)


#output highest Zscore phenos and most frequent z score phenos
#Merge the all data sets together
all.outliers <- rbind(gait.outliers, 
                      morpho.outliers,
                      JABS.outliers)

if(nrow(all.outliers >0)){
  #Pull out and plot the 16 phenotypes with the highest abs(Z score)
  highest.zscore.vals <- all.outliers |> 
    group_by(phenotype) |> 
    summarize(
      max_z = max(z_score)
    ) |> 
    arrange(desc(abs(max_z))) |> 
    slice_max(max_z, n = how.many.highest.zscore.outlier.plots) |>   
    select(phenotype)
  
  
  pdf(file.path(zscore.boxplot.dir, "highest_z_boxplot_outlier_figs.pdf"), 20, 20)
  print(
    inner_join(highest.zscore.vals, all.outliers) |> 
      plot_feature_outliers_zscore(threshold = zscore.threshold ,id_col_num = 2)
  )
  print(
    inner_join(highest.zscore.vals, all.outliers) |> 
      plot_feature_outliers_zscore(threshold = zscore.threshold ,id_col_num = 2)
  )
  dev.off()
  
  
  #Pull out and plot the 16 phenotypes with the most outliers
  most.freq.outliers <- all.outliers |> 
    group_by(phenotype) |> 
    summarize(
      count = sum(abs(z_score)>=zscore.threshold )
    ) |> 
    arrange(desc(count)) |> 
    slice_max(count, n = how.many.most.freq.outlier.plots) |> 
    select(phenotype)
  
  pdf(file.path(zscore.boxplot.dir, "most_freq_outliers_boxplot_outlier_figs.pdf"), 20, 20)
  print(
    inner_join(most.freq.outliers, all.outliers) |> 
      plot_feature_outliers_zscore(threshold = zscore.threshold ,id_col_num = 2)
  )
  print(
    inner_join(most.freq.outliers, all.outliers) |> 
      plot_feature_outliers_value(threshold = zscore.threshold ,id_col_num = 2)
  )
  dev.off()
}