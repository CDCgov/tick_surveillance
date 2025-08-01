#!/usr/bin/env Rscript
# Script for organizing the non-reference sequence blast assignment report
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
  seq_abundance_table         = args[2]
  metadata_file               = args[3]
  non_ref_assignments         = args[4]
  output_dir                  = "./"
} else {
  # else if running via RStudio
  r_libdir                    = "NA"
  seq_abundance_table         = "../results/dada2/sequence_abundance_table.tsv"
  metadata_file               = "../../2022_2_24_datasets/2021_4_group_2_v3.xlsx"
  non_ref_assignments         = "../results/non_reference_sequence_assignments.xlsx"
  output_dir                  = "../results/"
}


# these libraries are part of the tidyverse, so will be availabile in the
# tidyverse singularity image we are using (or analogous conda env)
library(tidyverse)
library(readxl)

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


# read in all needed files
sequence_abundance <- read.delim(seq_abundance_table, sep="\t", header=T)
non_ref_seq <- read_excel(non_ref_assignments)
non_ref_seq_df <- as.data.frame(non_ref_seq)
metadata <- read_excel(metadata_file )
non_ref_seq_wb <- loadWorkbook(non_ref_assignments)


# with sequence abundance table, remove all columns with 0 entires for abundance
seq_abundance = filter(sequence_abundance, abundance > 0)


# rename column name in seq_abundance df from dataset to Index
names(seq_abundance)[names(seq_abundance) == "dataset"] <- "Index"


# rename column name in non_ref_seq df from unassigned_sequence_name to sequence_number
names(non_ref_seq_df)[names(non_ref_seq_df) == "unassigned_sequence_name"] = "sequence_number"



# create new df that merges non_ref_seq and seq_abundance by sequence_number. 
seq_num_merge <- merge(seq_abundance, non_ref_seq_df, by = "sequence_number")


# remove column "unassigned_seuqnece" because its duplicate of "sequence"
#seq_num_merge <- subset(seq_num_merge, select = -c(unassigned_sequence))
seq_num_merge %>% select(-any_of(c("unassigned_sequence")))

# create new df that contains only subset of metadata
metadata_trim <- c("Index", "Pathogen_Testing_ID", "CSID")
metadata_trim_df = metadata[metadata_trim]


# create new df that merges metadata_trim_df with seq_num_merge
all_unassigned_report <- merge(seq_num_merge, metadata_trim_df, by = "Index")


# reorder columns of all_unassigned_report so all metadata values are in the beginning
all_unassigned_report <- all_unassigned_report[, c("Index", "Pathogen_Testing_ID", "CSID",
                                                  "abundance", "sequence_number", "scientific_name",
                                                  "common_name", "blast_name", "kingdom", "accession", 
                                                  "percent_identity", "alignment_length", "evalue", "sequence")]


# create new column that contains the bp lenght of the sequence
all_unassigned_report$sequence_length <- nchar(all_unassigned_report$sequence)


# re-order columns
all_unassigned_report <- all_unassigned_report[, c("Index", "Pathogen_Testing_ID", "CSID",
                                                  "abundance", "sequence_number", "scientific_name",
                                                  "common_name", "blast_name", "kingdom", "accession", 
                                                  "percent_identity", "alignment_length", "evalue", 
                                                  "sequence_length", "sequence")]

# group data by CSID, then sequence_number

all_unassigned_report <- all_unassigned_report[order(all_unassigned_report$Index, all_unassigned_report$sequence_number),]

# check if worksheet name already exists (if the pipeline is re-started after this process completed, the worksheet will need to be deleted for it to complete)

for (sheet_name in names(non_ref_seq_wb)) {
  if (sheet_name == "non_refseq_hits_with_metadata" ) {
    sheet_name <- TRUE
  }
}

if (sheet_name == TRUE) {
  removeWorksheet(non_ref_seq_wb, "non_refseq_hits_with_metadata")
  addWorksheet(non_ref_seq_wb, "non_refseq_hits_with_metadata")
}else {
  addWorksheet(non_ref_seq_wb, "non_refseq_hits_with_metadata")
  
}

# write to non_reference_sequence_assingments file

writeData(non_ref_seq_wb, "non_refseq_hits_with_metadata", all_unassigned_report)
saveWorkbook(non_ref_seq_wb, "non_reference_sequence_assignments.xlsx", overwrite = TRUE)