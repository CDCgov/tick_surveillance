#!/usr/bin/env Rscript
#
# This script assigns sequences that weren't assigned to any reference sequences based
# on blastn alignments vs. the nt database
#
# Mark Stenglein 1/13/2021 
#


#
# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#
if (!interactive()) {
  # if running from Rscript
  args = commandArgs(trailingOnly=TRUE)
  r_libdir                    = args[1]
  unassigned_sequences_fasta  = args[2]
  blast_output_path           = args[3]
  output_path                 = "./"
} else {
  # if running via RStudio
  r_libdir                    = "NA"
  unassigned_sequences_fasta  = "../test/results/blast/unassigned_sequences.fasta"
  blast_output_path           = "../test/results/blast/unassigned_sequences.fasta.bn_nt"
  output_path                 = "../test/results/"
}

library(tidyverse)
# load openxlsx, either from pipeline's R lib dir or from R environment
if (r_libdir != "NA") {
  library(openxlsx, lib.loc=r_libdir)
  library(zip,      lib.loc=r_libdir)
} else {
  library(openxlsx)
  library(zip)
}

# read in fasta-formatted unassigned sequences
unassigned_sequences = data.frame(query = character(),
                                  sequence = character())
seq_lines <- readLines(unassigned_sequences_fasta)
for (line in seq_lines){
  seq_id_match = str_match(line, ">(.*)")[,2]
  if (!is.na(seq_id_match)){
    seq_id = seq_id_match
  } else {
    sequence = line
    unassigned_sequences <- unassigned_sequences %>% add_row(query = seq_id, sequence = sequence)
  }
}

# read in the output file from blastn
blast_df_all <- read.delim(blast_output_path, sep="\t", header=T)

# query name type should be character not integer
blast_df_all <- blast_df_all %>% mutate(qaccver = as.character(qaccver))

# rename some of the oddly-named columns
blast_df_all <- blast_df_all %>% rename(query = qaccver,
                                        subject = saccver,
                                        percent_identity = pident,
                                        alignment_length = length)

# take the highest N scoring blast hits for each query
# this # can be modified
max_blast_hits_per_query <- 10
blast_df <- blast_df_all %>% 
  group_by(query) %>% 
  arrange(-bitscore, .by_group = T) %>% 
  filter(row_number() <= max_blast_hits_per_query)

# merge the sequence info and the blast output
blast_df <- left_join(blast_df, unassigned_sequences, by = "query")

# qaccver	saccver	pident	length	mismatch	gaps	qstart	qend	sstart	send	evalue	bitscore	staxid	ssciname	scomname	sblastname	sskingdom
# 12	MN086889.1	100.000	55	0	0	1	55	142	196	3.26e-18	102	373540	Borrelia lanei	Borrelia lanei	spirochetes	Bacteria
# 12	MK604329.1	100.000	55	0	0	1	55	349	403	3.26e-18	102	373540	Borrelia lanei	Borrelia lanei	spirochetes	Bacteria


# make a reporting dataset
# should contain these columns:
#
# *) sequence number
# *) sequence
# *) subject sci name
# *) subject common name
# *) subject blast name
# *) subject kingdom
# *) subject accession
# *) % alignment identity
# *) alignment length
# *) e-value

reporting_dataset <- blast_df %>% 
  select(query, 
         sequence, 
         ssciname, 
         scomname, 
         sblastname, 
         sskingdom, 
         subject,
	 percent_identity,
         alignment_length,
         evalue) %>%
  arrange(query)


# rename columns for better readability 
reporting_dataset <- reporting_dataset %>% 
  rename(
    unassigned_sequence_name = query,
    unassigned_sequence = sequence,
    scientific_name = ssciname,
    common_name = scomname,
    blast_name = sblastname,
    kingdom = sskingdom, 
    accession = subject)

# write data as plain-text
write.table(reporting_dataset, paste0(output_path, "non_reference_sequence_assignments.tsv"), quote=F, sep="\t", col.names=T, row.names=F)

# write as excel
wb <- createWorkbook("non_reference_sequence_assignments.xlsx")
addWorksheet(wb, "non_refseq_hits")
writeData(wb, "non_refseq_hits", reporting_dataset)
saveWorkbook(wb, paste0(output_path, "non_reference_sequence_assignments.xlsx"), overwrite = TRUE)

# write out version info into versions.yml
# TODO
writeLines("", "versions.yml")
