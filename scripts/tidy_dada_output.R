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
  r_bindir=args[1]
  dada_seqtab=args[2]
  outdir="./"
} else {
  # if running via RStudio
  r_bindir = "."
  trimmed_path = "../results/dada2/dada_seqtab.txt"
  outdir="../results/"
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
	       # mutate(sequence_number = row_number())

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





