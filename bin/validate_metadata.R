#!/usr/bin/env Rscript
#
# This script performs validation on input metadata and sample IDs.  
#
# Mark Stenglein 6/16/2022
#

library(tidyverse)
library(readxl)

# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#

# if running from Rscript (called from the pipeline)
if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  sample_metadata_file = args[1]
  sample_ids_file = args[2]
  output_dir = "./"
} else {
  # if running via RStudio
  sample_metadata_file = "../input/AK_metadata.xlsx"
  sample_ids_file = "../results/sample_ids.txt"
  output_dir = "../results/"
}


# This script: 
# - validates that the metadata file contains required columns
# - validates that the metadata file contains values for certain columns
# - validates that sample names (fastq file names) are not repeated
# - validates that metadata sample names (in Index column) are not repeated
# - confirms that there is one-to-one correspondence between fastq sample IDs and metadata IDs


# nextflow already checks that this file exists
metadata <- read_excel(sample_metadata_file)

# read in sample IDs from file in pipeline
# these are derived from the actual fastq file names
sample_ids <- read.delim(sample_ids_file, header=F)
colnames(sample_ids) <- c("sample_id")

# check that sample IDs are unique
repeated_sample_ids <- sample_ids %>% 
  group_by(sample_id) %>% summarize(n=n()) %>% filter(n>1) 

# error if any sample IDs repeated
if (nrow(repeated_sample_ids) > 0) {
  message(paste0("Error: repeated sample IDs in fastq filenames.  Repeated IDs include:"))
  print(repeated_sample_ids, max=nrow(repeated_sample_ids))
  quit(save="no", status=1)
}


# check that required metadata columns exist.  Error and abort if not.
# required_metadata_columns <- c("Index", "Not_A_Column_")
required_metadata_columns <- c("Index", 
			       "Pathogen_Testing_ID", 
			       "CSID", 
			       "Morphological_Ectoparasite_Genus_Species", 
			       "Lifestage", 
			       "State")

for (required_col in required_metadata_columns) {
  if (!required_col %in% colnames(metadata)) {
    message(paste0("Error: required metadata column: ", 
                   required_col, 
                   " not present in metadata spreadsheet: ",
                   sample_metadata_file))
    quit(save="no", status = 1)
  }
}

# check that metadata Indexes are unique
repeated_indexes <- metadata %>% 
  group_by(Index) %>% summarize(n=n()) %>% filter(n>1) 

# error if any sample IDs repeated
if (nrow(repeated_indexes) > 0) {
  message(paste0("Error: repeated Indexes in metadata spreadsheet: ", 
                 sample_metadata_file, 
                 ".  Repeated Indexes include:"))
  print(repeated_indexes, max=nrow(repeated_indexes))
  quit(save="no", status=1)
}


# There be a column in the metadata spreadsheet named "Index" 
# and there should be a 1-to-1 correspondence between the values in this column 
# and the fastq files supplied to the pipeline, which are listed in the sample IDS
# file from nextflow

# check that there is a metadata row in spreadsheet for all fastq
# exit if missing metadata
for (sample_id in sample_ids$sample_id) {
  if (!sample_id %in% metadata$Index){
    message(paste0("Error: there is no metadata for sample ID: ", sample_id,
                   " in metadata spreadsheet: ", 
                   sample_metadata_file)) 
    quit(save="no", status=1)
  }
}

# issue warning if there are not fastq for all samples defined in metadata 
for (index in metadata$Index) {
  if (!index %in% sample_ids$sample_id){
    message(paste0("Warning: there is no fastq files for one or more samples defined in the metadata spreadsheet.\n", 
                   " Metadata spreadsheet: ",  sample_metadata_file,
                   " Metadata Index: ", index)) 
  }
}

message("metadata and sample ID check complete.")
