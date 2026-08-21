
#A script to manage merged NextFlow data for final analysis
#Developed by Dr. Jake Beierle (don't forget the Dr., it's important)
#10-26-2025

#Documentation----
#See comprehensive documentation on the github repository
#https://github.com/jacobbeierle/JABS_nextflow_postprocess/tree/main


#----Load Libraries----

library(correlation)
library(janitor)
library(ppcor)
library(ComplexHeatmap)
library(colorRamp2)
library(tidyverse)

options(error = NULL)

#########################################################Define the working directory------------------------------------------------------------------

##If you are not using an R project, set your working directory
#working.directory <- "C:\\Users\\beierj\\Desktop\\2025-04-09_NTG_C1-C5_Analysis\\"

##Metadata for partial and regular correlations
metadata.to.retain <- c("sex", "day", "tx")


#Create functions--------------------------------------------------------------
#Creates a directory if it does not already exist
directory_check_creation <- function(x){
  x <- paste0(x)
  if(!dir.exists(x)){
    dir.create(x)
  }
  else{print(paste0("'/", substitute(x), "' already exists"))}
}

#I will document these when I am done developing them, TUAN. YEAH BUD. IM CALLING YOU OUT.
cor.to.list <- function(phenotype, correlate) {
  temp_corr_df <- data.frame(correlate, phenotype)
  colnames(temp_corr_df) <- c("correlate", "phenotype")
  temp_corr_df <- temp_corr_df[!is.na(temp_corr_df[, 2]), ]
  
  correlation <- cor.test(temp_corr_df$correlate, 
                          temp_corr_df$phenotype)
  list(
    pval = setNames(correlation[["p.value"]], colnames(temp_corr_df)[2]),
    r = setNames(correlation[["estimate"]][["cor"]], colnames(temp_corr_df)[2])
  )
}

pcor.to.list <- function(phenotype, correlate, control) {
  temp_corr_df <- data.frame(correlate, control, phenotype)
  colnames(temp_corr_df) <- c("correlate", "control", "phenotype")
  temp_corr_df <- temp_corr_df[!is.na(temp_corr_df[, 3]), ]
  
  correlation <- pcor.test(temp_corr_df$correlate,
                           temp_corr_df[3],
                           temp_corr_df$control)
  list(
    pval = setNames(correlation$p.value, colnames(temp_corr_df)[3]),
    r = setNames(correlation$estimate, colnames(temp_corr_df)[3])
  )
}

rm_df_to_cor_output_list <- function(df, rm, corr_var){
  output <- NULL
  
  #Create a list of unique values for repeated measure
  repeated.measure <- sort(unique(df[,rm]))
  
  #Ensure the cols are ordered the way you'd like below
  df <- cbind(df[c(rm, corr_var)],
              df[!(names(df) %in% c(rm, corr_var))])
  
  colnames(df)[1:2] <- c(c("rm", "corr"))
  
  current_rm_value = 1
  for(current_rm_value in repeated.measure) {
    temp_rm_df <- df[df$rm == current_rm_value,]
    temp_rm_df <- temp_rm_df[colnames(temp_rm_df) != "rm"]
    
    #Remove cols with less than 3 observations/level
    #Remove cols with less than 3 observations/level
    counts <- temp_rm_df |> 
      group_by(corr) |> 
      summarise(across(everything(), ~sum(!is.na(.))), .groups = "drop")
    #Keep columns where all groups have ≥ 3 non-NA values
    valid_cols <- names(counts)[colSums(counts >= 3) == nrow(counts)]
    temp_rm_df <- temp_rm_df[c("corr", valid_cols)]
    
    #Remove cols with 0 var
    nonzero_variance_cols <- sapply(temp_rm_df, function(x) var(x, na.rm = TRUE)) == 0
    #Replace the NAs with FALSE, as these represent strings etc
    nonzero_variance_cols[is.na(nonzero_variance_cols)] <- FALSE
    sum(nonzero_variance_cols)
    # Subset the data frame to remove zerovar cols
    temp_rm_df <- temp_rm_df[, !nonzero_variance_cols]
    
    
    results <- lapply(temp_rm_df[-c(1)], function(x) cor.to.list(x, temp_rm_df$corr))
    
    results <- as.data.frame(do.call(rbind, lapply(results, as.data.frame))) |> 
      rownames_to_column("measure") |> 
      mutate(!!sym(rm) := current_rm_value)
    
    if(length(output)){
      output <- rbind(output, results)
    }else{output <- results}
    
  }
  
  output |> 
    group_by(day) |> 
    summarize(
      count = n()
    )
  
  output <- list(r = output[c("measure", "r", paste(rm))], pval = output[c("measure", "pval", paste(rm))] )
  
  #Creat FDR list
  output$FDR <- output$pval 
  output$FDR$pval <- p.adjust(output$pval$pval, method = "fdr")
  colnames(output$FDR)[colnames(output$FDR)=="pval"] <- "FDR"
  
  #Make lists dfs
  output$r <- output$r |> 
    pivot_wider(
      id_cols = measure,
      names_from = !!sym(rm),
      values_from = r
    )
  output$pval <- output$pval |> 
    pivot_wider(
      id_cols = measure,
      names_from = !!sym(rm),
      values_from = pval
    )
  output$FDR <- output$FDR |> 
    pivot_wider(
      id_cols = measure,
      names_from = !!sym(rm),
      values_from = FDR
    )
  #make measure col into rownames for complex heatmapts
  output <- lapply(output, function(df) {
    column_to_rownames(df, var = "measure")})
  output
}

rm_df_to_pcor_output_list <- function(df, rm, corr_var, cont_var){
  output <- NULL
  
  #Create a list of unique values for repeated measure
  repeated.measure <- sort(unique(df[,rm]))
  
  #Ensure the cols are ordered the way you'd like below
  df <- cbind(df[c(rm, corr_var, cont_var)],
              df[!(names(df) %in% c(rm, corr_var, cont_var))])
  
  colnames(df)[1:3] <- c(c("rm", "corr", "cont"))
  
  current_rm_value = 0
  for(current_rm_value in repeated.measure) {
    temp_rm_df <- df[df$rm == current_rm_value,]
    temp_rm_df <- temp_rm_df[colnames(temp_rm_df) != "rm"]
    
    #Remove cols with less than 3 observations/level
    counts <- temp_rm_df |> 
      group_by(corr, cont) |> 
      summarise(across(everything(), ~sum(!is.na(.))), .groups = "drop")
    #Keep columns where all groups have ≥ 3 non-NA values
    valid_cols <- names(counts)[colSums(counts >= 3) == nrow(counts)]
    temp_rm_df <- temp_rm_df[c("corr", "cont", valid_cols)]
    
    #Remove cols with 0 var
    nonzero_variance_cols <- sapply(temp_rm_df, function(x) var(x, na.rm = TRUE)) == 0
    #Replace the NAs with FALSE, as these represent strings etc
    nonzero_variance_cols[is.na(nonzero_variance_cols)] <- FALSE
    sum(nonzero_variance_cols)
    # Subset the data frame to remove zerovar cols
    temp_rm_df <- temp_rm_df[, !nonzero_variance_cols]
    
    
    results <- lapply(temp_rm_df[-c(1,2)], function(x) pcor.to.list(x, temp_rm_df$corr, temp_rm_df$cont))
    
    results <- as.data.frame(do.call(rbind, lapply(results, as.data.frame))) |> 
      rownames_to_column("measure") |> 
      mutate(!!sym(rm) := current_rm_value)
    
    if(length(output)){
      output <- rbind(output, results)
    }else{output <- results}
  }
  
  output <- list(r = output[c("measure", "r", paste(rm))], pval = output[c("measure", "pval", paste(rm))] )
  
  #Creat FDR list
  output$FDR <- output$pval 
  output$FDR$pval <- p.adjust(output$pval$pval, method = "fdr")
  colnames(output$FDR)[colnames(output$FDR)=="pval"] <- "FDR"
  
  #Make lists dfs
  output$r <- output$r |> 
    pivot_wider(
      id_cols = measure,
      names_from = !!sym(rm),
      values_from = r
    )
  output$pval <- output$pval |> 
    pivot_wider(
      id_cols = measure,
      names_from = !!sym(rm),
      values_from = pval
    )
  output$FDR <- output$FDR |> 
    pivot_wider(
      id_cols = measure,
      names_from = !!sym(rm),
      values_from = FDR
    )
  #make measure col into rownames for complex heatmapts
  output <- lapply(output, function(df) {
    column_to_rownames(df, var = "measure")})
  output
}

plot_heatmap <- function(list, threshold, adjusted = TRUE, min_sig = 1, min_sig_exclude_row_1 = FALSE,
                         height = 25, width = 7.5, col_titles){
  #Select only significant Pvals of FDRs
  if(adjusted){
    keep <- data.frame(list$FDR <= threshold)
    plot <- list$r * keep
  }else{
    keep <- data.frame(list$pval <= threshold)
    plot <- list$r * keep
  }
  
  #Remove rows without enough significant correlations
  if(min_sig_exclude_row_1){
    non_zero_counts <- as.vector(rowSums(plot[-1] != 0, na.rm = TRUE) >= min_sig)
    plot <- plot[non_zero_counts,]
  }else{
    non_zero_counts <- as.vector(rowSums(plot != 0, na.rm = TRUE) >= min_sig)
    plot <- plot[non_zero_counts,]
  }
  
  #Set color schemes based on max abs(r) value
  plot.val.extreme <- max(abs(plot), na.rm = TRUE)
  correlation.color = colorRamp2(c(-plot.val.extreme, 0, plot.val.extreme), 
                                 c("darkred", "white", "darkblue"))
  
  #Create the plot
  heatmap.out <-  Heatmap(as.matrix(plot), cluster_rows = TRUE, cluster_columns = FALSE, 
                          show_row_dend = FALSE, row_dend_reorder = TRUE,
                          width = unit(width, "cm"), height = unit(height, "cm"),
                          name = "r", col = correlation.color,
                          border = TRUE, rect_gp = gpar(col = "slateblue", lwd = 1),
                          row_title =  "", row_names_side = "left", row_names_gp = gpar(fontsize = 10),
                          column_title = col_titles, column_title_side = "bottom", column_names_rot = 90)
  heatmap.out
}
#Set working directory, and check for expected file names and directories---------------------------------

#If you have defined a working directory above, set it here
if(exists("working.directory")){
  setwd(working.directory)
}

directory_check_creation(file.path("exploratory_figs"))
directory_check_creation(file.path("exploratory_figs", "heatmaps"))
heatmap.pub.dir <- file.path("exploratory_figs", "heatmaps")

merged_nextflow_dataset.dir <- file.path("Nextflow_Output", "final_nextflow_feature_data")


#----Import and preprocess data-----

#Import
data.nextflow <- list.files(pattern = "_final_nextflow_dataset\\.csv") |> 
  read.csv() |> 
  clean_names()

#Select non-fecal boli data and metadata needed for this analysis

#Read in gait colnames for subsetting dataframe
gait.cols <- list.files(
  path = merged_nextflow_dataset.dir,
  pattern = "gait_final",
  full.names = TRUE) |> 
  read_csv(n_max = 1) |> 
  clean_names() |> 
  colnames()

#Read in gait colnames for subsetting dataframe
morpho.cols <- list.files(
  path = merged_nextflow_dataset.dir,
  pattern = "morphometrics_final",
  full.names = TRUE) |> 
  read_csv(n_max = 1) |> 
  clean_names()|> 
  colnames()

#Read in gait colnames for subsetting dataframe
JABS.cols <- list.files(
  path = merged_nextflow_dataset.dir,
  pattern = "features_final",
  full.names = TRUE) |> 
  read_csv(n_max = 1) |> 
  clean_names()|> 
  colnames()

data.nextflow <- data.nextflow[colnames(data.nextflow) %in% metadata.to.retain|
                               colnames(data.nextflow) %in% gait.cols[-1] |
                               colnames(data.nextflow) %in% morpho.cols[-1] |
                               colnames(data.nextflow) %in% JABS.cols[-1]]

#colnames(data.nextflow)[colnames(data.nextflow) == "day"] <- "test.day"
#Make sex a bianary col representing # of Y chromosomes
data.nextflow$sex <- if_else(data.nextflow$sex == "F", 0, 1)

data.nextflow.f <- subset(data.nextflow, sex == 0) |> 
  select(!sex)
  
data.nextflow.m <- subset(data.nextflow, sex == 1) |> 
  select(!sex)


#----Loop for heatmaps----
corr.all.tx <- rm_df_to_pcor_output_list(data.nextflow, "day", "tx", "sex")

corr.all.tx.f <- rm_df_to_cor_output_list(data.nextflow.f, "day", "tx")


#Make some heatmaps-----
ht.all.corr.05.noadj <- plot_heatmap(corr.all.tx, 0.05, adjusted = FALSE,
                                     min_sig=2, min_sig_exclude_row_1 = TRUE, col_titles = "Day",
                                     height = 7.5, width = 20)

pdf(file="exploratory_figs/Heatmaps/test.pdf", width = 7.5, height = 12.5)
#png(file="exploratory_figs/Heatmaps/All_mice_tx_corr_05_pval.png", width = 17.5, height = 30, units= "cm", res = 600)
draw(test, column_title = "NTG:Behavior Correlation, all mice", column_title_gp = gpar(fontsize = 25))
dev.off()



