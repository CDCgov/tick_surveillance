#!/usr/bin/env Rscript
# script for filtering the scientific name column from unassigned_sequence_report_all tab in sequencing_report by input string, if parameter is given
#
# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
#
# if running from Rscript (called from the pipeline)
if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  r_libdir                    = args[1]
  seq_report                  = args[2]
  non_ref_assignments         = args[3]
  filter_params               = args[4]
  output_dir                  = "./"
} else {
  # else if running via RStudio
  r_libdir                    = "NA"
  seq_abundance_table         = "../results/dada2/sequence_abundance_table.tsv"
  non_ref_assignments         = "../results/non_reference_sequence_assignments.xlsx"
  filter_params               = "NA"
  output_dir                  = "../results/"
}


# these libraries are part of the tidyverse, so will be availabile in the
# tidyverse singularity image we are using (or analogous conda env)
library(tidyverse)
library(readxl)
library(dplyr)

# openxlsx is not part of standard tidyverse, so may have to load it 
# from a specified path
# load openxlsx, either from pipeline's R lib dir or from R environment
if (r_libdir != "NA") {
  library(openxlsx, lib.loc=r_libdir)
  library(zip, lib.loc=r_libdir)
  
} else {
  library(openxlsx)
  library(zip)
}



# read in sequencing report wb
seq_report_wb2 <- loadWorkbook(seq_report)

# read in all_unassinged_seq tab from sequencing report
all_unassigned = read_excel(non_ref_assignments, sheet = "non_refseq_hits_with_metadata")


# read in target filter params
select_filter_pattern <- c(filter_params)
select_filter_pattern <- select_filter_pattern %>% str_replace_all("\\[|\\]", "")

# split string into list
split_pattern <- as.list(strsplit(select_filter_pattern,split = ",")[[1]])


#filter
filterd_unassigned_report <- filter(all_unassigned, grepl(paste(split_pattern, collapse ='|'), scientific_name, ignore.case = TRUE))

for (sheet_name in names(seq_report_wb2)) {
  if (sheet_name == "FIILTER_unassigned_seqs" ) {
    sheet_name <- TRUE
  }
}

if (sheet_name == TRUE) {
  removeWorksheet(seq_report_wb2, "FIILTER_unassigned_seqs")
  addWorksheet(seq_report_wb2, "FIILTER_unassigned_seqs")
}else {
  addWorksheet(seq_report_wb2, "FIILTER_unassigned_seqs")
  
}

# write to new tab in sequencing report
#addWorksheet(seq_report_wb2, "FIILTER_unassigned_seqs")

if(nrow(filterd_unassigned_report) == 0){
  none_identified <- paste(c("The filter params (" , select_filter_pattern, ") were not identified in the scientific_name column of the unassigned sequence report"), collapse = "")
  writeData(seq_report_wb2, "FIILTER_unassigned_seqs", none_identified)
}else{
  writeData(seq_report_wb2, "FIILTER_unassigned_seqs", filterd_unassigned_report)
}
saveWorkbook(seq_report_wb2, "sequencing_report.xlsx", overwrite = TRUE)