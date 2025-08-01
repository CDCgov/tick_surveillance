#!/usr/bin/env Rscript 
#
# This script tidies dada2 output
#

library(tidyverse)
library(openssl)

#
# This code block sets up input arguments to either come from the command line
# (if not running in interactive mode) or to use expected default values 
# (if running in interactive mode, i.e. RStudio).  
#
# This is useful in case you want to modify this script
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  dada_seqtab=args[1]
  track_reads=args[2]
  outdir="./"
} else {
  # if running via RStudio
  trimmed_path = "../results/dada2/dada_seqtab.txt"
  track_reads= "../results/dada2/dada_read_clean_all.csv"
  outdir="../results/dada2"
}

# read in dada2 output
t <- read.delim(dada_seqtab, sep="\t", header=T)

# create a tidy formatted table with:
# unique, merged sequences and their abundances in each dataset
tidy_sequence_table <- pivot_longer(t, -dataset, names_to = "sequence", values_to = "abundance" )

# create a unique sequence id (#) for each sequence: for keeping track of in further processing steps
# this sequence ID will be a sha1 hash of the actual sequence
# this will facilitate comparison of sequences across datasets
sequences <- tidy_sequence_table %>% 
	       group_by(sequence) %>% 
	       summarize () %>% 
	       mutate(sequence_number = sha1(sequence))

# join these sequence #s back into the tidy table
tidy_sequence_table <- left_join(tidy_sequence_table, sequences, by="sequence")

# write out the tidy-formatted table
write.table(tidy_sequence_table, paste0(outdir, "sequence_abundance_table.tsv"), sep="\t", col.names=T, row.names=F)

# this writes fasta from the sequences data frame 
# see: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"sequence_number"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"sequence"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# write out the sequences in fasta format
writeFasta(sequences, paste0(outdir, "observed_sequences.fasta"))

# write out version info into versions.yml
writeLines("", "versions.yml")

# read in dada read tracking csv
track_reads <- read.csv("dada_read_clean_all.csv")

# Add column that output the percentage of reads that passed dada2 cleanup
track_reads$PercentPass <- (track_reads$nonchim/track_reads$input)*100

# write read tracking with % passing to new csv
write.csv(track_reads, paste0(outdir, "dada_read_clean_all_with_pct_pass.csv"), col.names=T)

# Create summary stats table
track_sum <- track_reads %>%
  summarise(across(where(is.numeric), .fns =
                     list(Min = min,
                          Median = median,
                          Mean = mean,
                          Max = max,
                          SD = sd))) %>%
  pivot_longer(everything(), names_sep = "_", names_to = c(".value", "Variable"))

# write read tracking to csv
write.csv(track_sum, paste0(outdir, "dada_read_clean_summary.csv"), col.names=T)
