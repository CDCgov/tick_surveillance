#!/usr/bin/env Rscript 

# This script:
#
#   1. Assigns observed sequences to a set of reference sequences based
#      on blast alignments 
#   2. Populates a surveillance reporting table
#   3. Creates additional, more detailed output
#
# Mark Stenglein 05/19/2022
#


# This code block sets up input arguments to either come from the command line
# (if running from the pipeline, so not interactively) or to use expected default values 
# (if running in interactive mode, for example in RStudio for troubleshooting
# or development).  
#

# if running from Rscript (called from the pipeline)
if (!interactive()) {
  args = commandArgs(trailingOnly=TRUE)
  r_libdir                    = args[1]
  tidy_table_path             = args[2]
  blast_output_path           = args[3]
  sample_metadata_file        = args[4]
  targets_tsv_file            = args[5]
  surveillance_columns_file   = args[6]
  input_min_non_control_reads = args[7]
  output_dir                  = "./"
} else {
  # else if running via RStudio
  r_libdir                    = "../lib/R/"
  tidy_table_path             = "../results/dada2/sequence_abundance_table.tsv"
  blast_output_path           = "../results/blast/observed_sequences.fasta.bn_refseq"
  sample_metadata_file        = "../../2022_2_24_datasets/2021_4_group_2_v3.xlsx"
  targets_tsv_file            = "../refseq/targets.tsv"
  surveillance_columns_file   = "../refseq/surveillance_columns.txt"
  input_min_non_control_reads = 50
  output_dir                  = "../results/"
}

# these libraries are part of the tidyverse, so will be availabile in the
# tidyverse singularity image we are using (or analogous conda env)
library(tidyverse)
library(lubridate)
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

# a function to write out sequences in fasta format
# takes as input vectors of sequence names and the sequences themselves
# adapted from: https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html
writeFasta <- function(sequence_headers, sequences, filename){
  fastaLines = c()
  if (length(sequence_headers) != length(sequences)) {
    message (paste0("ERROR: was expecting equal numbers of sequences and sequence headesr"))
    quit (status = 1)
  }
  for (i in 1:length(sequence_headers)){
    fastaLines = c(fastaLines, as.character(paste0(">", sequence_headers[i])))
    fastaLines = c(fastaLines, as.character(sequences[i]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}


# ------------------
# read input files
# -----------------

# ----------------------------
# sequence abundance table
# ----------------------------
# read in the sequence_abundance_table.tsv file, which contains abundances (read counts) for
# each sequence in all the datasets
sequence_abundance_table <- read.delim(tidy_table_path, sep="\t", header=T)

# make sure the dataset names have character types, in case they have integer value names (1, 2, ...)
sequence_abundance_table$dataset <- as.character(sequence_abundance_table$dataset)

# make a dataframe containing the sequences and their sequence number (which matches the blast query column)
# and add a column listing their length
sequences <- sequence_abundance_table %>% 
  group_by(sequence, sequence_number) %>% 
  summarize(.groups = "drop") %>%
  mutate(sequence_length = str_length(sequence)) %>% 
  rename(observed_sequence = sequence) %>% 
  ungroup()

# -------------
# blast input
# -------------

# read in the output file from blastn
blast_df_all <- read.delim(blast_output_path, sep="\t", header=T)

# rename some of the oddly-named columns 
blast_df_all <- blast_df_all %>% rename(query = qaccver,
                                        subject = saccver,
                                        percent_identity = pident,
                                        alignment_length = length,
                                        query_length = qlen,
                                        subject_length = slen)

# take the highest scoring blast hit for each query
blast_df <- blast_df_all %>% group_by(query) %>% arrange(-bitscore, .by_group = T) %>% filter(row_number() == 1)

# possibly TODO: Deal with the case when a sequence hits multiple refseqs equally well 
# This would be the case when a sequence is equally close to 
# two of the reference sequences.  E.g. a Borrelia burgdorferi equally similar to B31 and N40
blast_df_equal_hits <- blast_df_all %>% 
  group_by(query) %>% 
  mutate(max_bitscore = max(bitscore)) %>% 
  filter(bitscore == max_bitscore) %>% 
  summarize(n=n()) %>% 
  filter(n>1)

# merge the sequence info and the blast output
blast_df <- left_join(blast_df, sequences, by=c("query" = "sequence_number")) 

# 
# what fraction of the query was present in the alignment?
# 
# Scenarios for insertions or deletions in the alignment:
# 
# 1) If observed sequences has no indels relative to the reference
#   - alignment_length == query length
#   - gaps = 0
#
# 2) If observed sequences has a deletion relative to the reference
#   - alignment_length == (query_length + gaps)
#   - gaps = # of bases in deletion
#
# 3) If observed sequences has an insertion relative to the reference
#   - alignment_length == query length
#   - gaps = # of bases in insertion
#
# Note that any indels relative to reference sequence have the effect of decreasing 
# the % identity of the alignment in blast output in addition to influencing 
# gaps and alignment length values as described above.  So, any gaps in the alignment
# will produce a lower-scoring alignment less likely to be above cutoff.
#  

blast_df <- blast_df %>% mutate(percent_query_aligned = 100 * alignment_length / query_length,
                                percent_of_alignment_gaps = 100 * gaps / alignment_length)

# ----------
# metadata 
# ----------

# read in the sample metadata file
metadata_df <- read_excel(sample_metadata_file)

# keep track of original metadata names
original_metadata_names <- names(metadata_df)

# fix the metadata names using tidyverse make.names function
# in case some metadata column names have disallowed characters like - or + 
# (disallowed for tidyverse column names)
names(metadata_df) <- make.names(names(metadata_df),unique = TRUE)

# Make sure the Index column is character data, in case Indexes have integer values (1, 2, etc.)
metadata_df$Index <- as.character(metadata_df$Index)

# assign datasets to batches based on an optional metadata column named batch
if ( ! "batch" %in% colnames(metadata_df)) {
  # if no batch column specified in the metadata,
  # then consider all of the datasets as belonging to a single batch
  message ('INFO: There is no "batch" column in the metadata file so will consider all datasets to belong to a single batch')
  metadata_df <- metadata_df %>% mutate(batch = "1")
}

# Metadata columns in Excel date format get output as in integer that is the # of days since Jan 1, 1900 
# We want to make sure to output them as human-readable dates
# which columns have POSIXct date format?
posix_date_columns <- which(sapply(metadata_df, is.POSIXct))
# This converts columns with POSIXct format into Date format
for (date_col in posix_date_columns) {
  metadata_df[,date_col] <- lapply (metadata_df[,date_col], as.Date, format="yyyy-mm-dd")
}

# --------------------------
# Read in target information
# --------------------------
targets_df <- read.delim(targets_tsv_file, sep="\t", header=T)
targets_df$internal_control <- as.logical(targets_df$internal_control)

# rename target sequences to be clear that they include primer sequences, 
# which are stripped off of the actual observed sequences 
# since primer-derived bases in sequence reads are unreliable 
# (they derive from the primers and not the template)
targets_df <- targets_df %>% rename(target_sequence_incl_primers = sequence)

# convert empty strings in reporting column to NA values
targets_df$reporting_columns <- na_if(targets_df$reporting_columns, "")

# join target info with blast info
# This is necessary because the targets file contains info about minimum percent
# identity for an alignment between an observed an expected sequence, etc. 
blast_df <- left_join(blast_df, targets_df, by=c("subject" = "ref_sequence_name"))

# ----------------------------
# Read in surveillance columns
# ----------------------------
surveillance_columns <- read.delim(surveillance_columns_file, sep="\t", header=T, stringsAsFactors = F)

# --------------------------------------------------------------------------------------
# This input parameter specifies how many reads will be required to make a positive call 
# in the surveillance table
# --------------------------------------------------------------------------------------
# make sure it's not a character type but a numeric type
input_min_non_control_reads <- as.numeric(input_min_non_control_reads)

# --------------------------------------------------------------------------
# logic for saying that one of the observed sequences is close enough 
# to one of the reference sequences to assign it to that reference sequence
# --------------------------------------------------------------------------
assign_blast_hits_to_refseqs <- function() {
  
  # <<- is the global assignment operator 
  # (it modifies the global blast_df variable instead of creating a new variable 
  # with the same name in this function scope)
  blast_df <<- blast_df %>%  mutate(
    
    # assigned_to_target will be TRUE if this particular sequence
    # aligns to the relevant reference sequence above the thresholds
    # Note that the min_percent_identity and min_percent_aligned values 
    # are columns in the dataframe: they come from the targets_tsv_file
    # so can be assigned separately for each target
    assigned_to_target = if_else( (  percent_identity >= min_percent_identity & 
                                     percent_query_aligned >= min_percent_aligned &
                                     percent_of_alignment_gaps <= max_percent_gaps),
                                  TRUE, 
                                  FALSE))
}  


# Call function to see if sequences are close enough to a target to be assigned to that target
assign_blast_hits_to_refseqs()

# plot hits: percent identity vs. fraction aligned
assigned_plot <- ggplot(blast_df) +
  geom_point(aes(x=percent_query_aligned, 
                 y=percent_identity, 
                 fill=assigned_to_target),
             shape=21, alpha=0.8, size=2) +
  scale_fill_manual(values=c("red", "darkslateblue")) +
  theme_bw() + 
  facet_wrap(~subject) +
  xlab("% of query aligned to reference sequence") +
  ylab("% identity between query and reference sequence") 

# render the plot and output to PDF
ggsave(paste0(output_dir, "/assigned_targets_plot.pdf"), plot = assigned_plot, width=10, height=7.5, units="in")

# --------------------------------------------------------------------------------------
# separate assigned and unassigned sequences and write out unassigned sequences in fasta
# --------------------------------------------------------------------------------------
# assigned sequences
assigned_sequences <- blast_df %>% filter(assigned_to_target)

# unassigned sequences are those that weren't assigned to any target reference sequence
unassigned_sequences <- filter(sequences, !(sequence_number %in% assigned_sequences$query))

# write out the unassigned sequences in fasta format
# these will be dealt with separately by downstream processes
writeFasta(unassigned_sequences$sequence_number, 
           unassigned_sequences$observed_sequence, 
           paste0(output_dir, "unassigned_sequences.fasta"))


# get rid of 0 counts in sequence abundance table
sparse_sat <- filter(sequence_abundance_table, abundance > 0) 

# join sparse sat with columns from blast_df indicating whether the sequences
# have been assigned to a refseq 
sat_with_assigned <- left_join(sparse_sat, 
                               select(blast_df, query, assigned_to_target), 
                               by=c("sequence_number" = "query"))  

# turn NA values in assigned_to_target, resulting from queries producing no blast hits at all, 
# into F values
sat_with_assigned$assigned_to_target <- replace(sat_with_assigned$assigned_to_target, 
                                                is.na(sat_with_assigned$assigned_to_target), 
                                                FALSE)

# calculate fraction of reads in each dataset that are assigned or not
assigned_unassigned_counts <- sat_with_assigned %>% 
  group_by(dataset, assigned_to_target) %>% 
  summarize(reads = sum(abundance), .groups="drop") %>%
  ungroup() %>%
  group_by(dataset) %>%
  mutate(total_reads = sum(reads),
         fraction = reads / total_reads)
  
# clean-up dataframe
fraction_assigned <- assigned_unassigned_counts %>% 
  filter(assigned_to_target == T) %>%
  select(-assigned_to_target) %>%
  rename(assigned_reads = reads)

# output dataframe
write.table(fraction_assigned, file=paste0(output_dir, "fraction_reads_assigned.txt"),  
            quote=F, sep="\t", row.names=F, col.names=T)

# calculate_stats about numbers of reads assigned or not
assigned_stats <- fraction_assigned %>% ungroup() %>% summarize(mean_fraction_assigned = mean(fraction),
                                                  median_fraction_assigned = median(fraction),
                                                  sd_fraction_assigned = sd(fraction),
                                                  min_fraction_assigned = min(fraction),
                                                  max_fraction_assigned = max(fraction))

# output dataframe with stats
write.table(assigned_stats, file=paste0(output_dir, "fraction_reads_assigned_summary_stats.txt"),  
            quote=F, sep="\t", row.names=F, col.names=T)

# create a histogram of fraction of reads assigned
assigned_histogram = ggplot(fraction_assigned) +
  geom_histogram(aes(x=fraction), bins=50, fill="darkslateblue", color="black", size=0.15) +
  theme_classic() +
  xlab("Fraction of reads assigned to a reference sequence") + 
  ylab("Datasets")

ggsave(paste0(output_dir, "/Fraction_of_reads_assigned_histogram.pdf"), plot = assigned_histogram, width=10, height=7.5, units="in")


# -----------------------
# consolidate dataframes
# -----------------------

# keep track of metadata rows
metadata_key <- metadata_df %>% select(Index, Pathogen_Testing_ID, batch)

# join sparse SAT with metadata
# the resulting dataframe, dataset_df, initially just has basic metadata and 
# sequence info.  We will join in a lot of more info and this will be the 
# main table keepting track of results
dataset_df <- full_join(metadata_key, sparse_sat, by=c("Index" = "dataset"))  

# join blast results into dataset_df
dataset_df <- left_join(dataset_df, blast_df, by=c("sequence_number" = "query")) %>% 
  filter(assigned_to_target == T)  %>%
  # get rid of one of two identical redundant sequence columns
  select(-observed_sequence) %>%
  rename(observed_sequence = sequence)

# ---------------------------------------------------------------------
# QC criterion: minimum # of reads mapping to internal pos. control
# this corresponds to the Acceptable DNA column in the reporting data
# ---------------------------------------------------------------------

# calculate the log10 # of tick internal-control mapping reads in datasets by batch
internal_control_batch_averages <- dataset_df %>%
  filter(internal_control == TRUE) %>%
  group_by(batch, Index) %>%
  mutate(per_sample_internal_ctrl_reads = log10(sum(abundance))) %>%
  ungroup() %>%
  group_by(batch) %>%
  summarize(mean_batch_internal_ctrl_reads = mean(per_sample_internal_ctrl_reads),
            sd_batch_internal_ctrl_reads   = sd  (per_sample_internal_ctrl_reads))

# -------------------------------------------------
# Minimum number of internal control mapping reads
# -------------------------------------------------

# for internal control targets, this is the mean of the batch log10 read # - 3 standard deviations
internal_control_batch_averages <- internal_control_batch_averages %>%
  mutate(minimum_internal_control_log_reads =
           mean_batch_internal_ctrl_reads - (3 * sd_batch_internal_ctrl_reads))

# join in info about internal control cutoffs
dataset_df <- left_join(dataset_df, internal_control_batch_averages, by="batch")

# create a plot of internal control (tick actin) reads in individual datasets
tick_reads_p <- dataset_df %>% 
  filter(internal_control) %>% 
  group_by(batch, Index, species) %>%
  summarize(abundance = sum(abundance), above_cutoff = log10(abundance) > minimum_internal_control_log_reads) %>%
  ggplot() +
  geom_violin(aes(x=batch, y=abundance), color="slateblue", fill="NA", size=0.5) +
  geom_jitter(aes(x=batch, y=abundance, fill=above_cutoff),
              shape=21, size=3, stroke=0.15, width=0.2) +
  scale_fill_manual(values=c("red", "darkslateblue")) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Reads mapping to tick internal control targets") +
  xlab("Batch") +
  ggtitle("Number of tick-mapping reads in individual tick or control datasets\n")

# plot output to PDF
ggsave(paste0(output_dir, "/tick_reads_per_dataset.pdf"), plot = tick_reads_p,  width=7, height=5, units="in")


# ----------------------------------------
# Generate and output surveillance table
# ----------------------------------------

# create map of targets to surveillance columns 
# this mapping will be used to populate values in surveillance report
targets_to_surv_df <- data.frame(target = character(),
                                 surveillance_column = character(),
                                 surveillance_column_type = character())
                                

# for each row in targets table
for (row in 1:nrow(targets_df)) {
  
  # what is the target name?
  target <- filter(targets_df, row_number() == row) %>% pull(ref_sequence_name)
  
  # what reporting columns are defined for this target?
  reporting_columns <- filter(targets_df, row_number() == row) %>% 
    pull(reporting_columns) %>% 
    str_split(pattern = ';', simplify = T)
  
  # collect targets for each reporting column into a named list
  for (col in reporting_columns) {
    
    # split the reporting column by colon: name:type
    # where type = "count" or "name"
    col_fields <- str_split(col, pattern = ":", simplify = T)
    col_name <- col_fields[,1]
    if (length(col_fields) > 1) {
      col_type <- col_fields[,2]
      if (!(col_type %in% c("count", "name"))) {
        message (paste0("ERROR: invalid column type for column ", col ," in target file (", targets_tsv_file , ")."))
        quit(status = 1)
      }
    } else {
      col_type <- NA_character_
    }

    if (!is.na(col_name)) {
      # double check that any surveillance column names defined in the targets file are also defined in the surveillance columns file
      if (!(col_name %in% surveillance_columns$column_name)) {
        message (paste0("ERROR: surveillance column ", col_name ," in target file (", targets_tsv_file , ")",
                        " not present in surveillance columns file (", surveillance_columns_file ,")"))
        quit(status = 1)
      }
      
      # map targets to surveillance columns
      targets_to_surv_df[nrow(targets_to_surv_df) + 1,] <- c(target, col_name, col_type)
    }
  }
}


# -------------------------------------------------------
# collapse unique sequences to level of targeted species
# this info will be output as its own table
# -------------------------------------------------------

# collapse all hits for same target species in each dataset 
# this will get output in the main excel output as a separate tab
dataset_by_spp <- dataset_df %>% 
  group_by(Index, species) %>%
  mutate(
    abundance = sum(abundance),   
    percent_identity = mean(percent_identity),
    percent_query_aligned = mean(percent_query_aligned),
    richness = n(),
  ) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  # get rid of unneeded columns
  select(Index, Pathogen_Testing_ID, batch, species, abundance, percent_identity, 
         percent_query_aligned, richness, internal_control, 
         minimum_internal_control_log_reads)

# ----------------------------------------------------------------------
# Collect information that will be used to populate surveillance report
# ----------------------------------------------------------------------

# join in surveillance column info into targets file
# because each target can be mapped to >1 surveillance column, this will result in rows with
# duplicate target info but multiple surveillance columns
dataset_with_surv_df <- left_join(dataset_df, targets_to_surv_df, by=c("subject" = "target"))

# create a column defining the minimum cutoff for non-internal control targets
# for non-internal control targets, this is a simple cutoff passed in as an argument 
# (e.g. 50 reads)
dataset_with_surv_df$minimum_non_control_reads = input_min_non_control_reads

# keep track of whether individual (dataset x target) combinations have sufficient reads to be
# called positive
dataset_with_surv_df <- dataset_with_surv_df %>% 
  mutate(
    individual_pos_neg_call = 
      case_when((internal_control & (log10(abundance) >= minimum_internal_control_log_reads)) ~ "Positive",
                (!internal_control & (abundance >= minimum_non_control_reads)) ~ "Positive",
                TRUE ~ "Negative"))

# make species-level calls
# first, calculate summed read counts by species
dataset_with_surv_df_species_level_calls <- dataset_with_surv_df %>% 
  filter(surveillance_column_type == "count") %>%
  group_by(Index, batch, species) %>%
  summarize(
    species_abundance             = sum(abundance),
    .groups="drop")

# join in species-level calls
dataset_with_surv_df <- left_join(dataset_with_surv_df, dataset_with_surv_df_species_level_calls)
    
# make species-level pos/neg calls using read counts summed by species
dataset_with_surv_df <- dataset_with_surv_df %>% 
  mutate(
    species_level_pos_neg_call = 
      case_when((internal_control & (log10(species_abundance) >= minimum_internal_control_log_reads)) ~ "Positive",
                (!internal_control & (species_abundance >= minimum_non_control_reads)) ~ "Positive",
                TRUE ~ "Negative"))


# ------------------
# name-type columns
# ------------------

# for each name-type surveillance column, keep track of which targets had sufficient #s 
# to be reported 
name_type_surveillance_values <- dataset_with_surv_df %>% 
  filter(surveillance_column_type == "name" & 
           species_level_pos_neg_call == "Positive") %>%
  group_by(Index, batch, surveillance_column) %>%
  summarize(
    contributing_target_names = str_trim(paste0(unique(species), sep = " ", collapse="")),
    .groups="drop"
  ) 

# -------------------
# count-type columns
# -------------------
# for each dataset and surveillance column, sum all read counts 
# for targets that contribute to that surveillance columns 
count_type_surveillance_values <- dataset_with_surv_df %>% 
  filter(surveillance_column_type == "count") %>%
  group_by(Index, batch, surveillance_column) %>%
  summarize(
    abundance             = sum(abundance),   
    percent_identity      = mean(percent_identity),
    percent_query_aligned = mean(percent_query_aligned),
    richness              = n(),
    internal_control      = any(internal_control),   
    minimum_non_control_reads = first(minimum_non_control_reads),
    .groups="drop"
  ) 

# join back in info about minimum # of internal control reads
count_type_surveillance_values <- left_join(count_type_surveillance_values, internal_control_batch_averages, by="batch")

# make positive negative calls for count-type surveillance columns
# based on summed abundances for all targets that contributed to that column
# and cutoffs for internal-control or non-control targets
count_type_surveillance_values <- count_type_surveillance_values %>% 
  mutate(pos_neg_call = 
           case_when(( internal_control & (log10(abundance) >= minimum_internal_control_log_reads)) ~ "Positive", 
                     (!internal_control & (abundance >= minimum_non_control_reads)) ~ "Positive", 
                     TRUE ~ "Negative"))


# keep only necessary columns
count_type_surveillance_columns <- count_type_surveillance_values %>% select(Index, batch, surveillance_column, abundance, pos_neg_call)

# ----------------------------
# populate surveillance table
# ----------------------------

# columns in this table are defined in the surveillance_columns.txt file (param --surveillance_columns)

# how big should the DF be?
num_surv_cols <- nrow(surveillance_columns)
num_surv_rows <- nrow(metadata_key)  

# create an empty DF
surv_df <- data.frame(mat = matrix(ncol = num_surv_cols, nrow = num_surv_rows), stringsAsFactors = F)

# name the columns
colnames(surv_df) <- t(surveillance_columns$column_name)

# confirm Index column present in surveillance report columns, since it's referred to in a hardcoded way below
if (!("Index" %in% colnames(surv_df))) {
  message (paste0("ERROR: \"Index\" column not present in surveillance columns file (", surveillance_columns_file ,")"))
  quit(status = 1)
}

# prepopulate the columns with any default values
#
# convert default values from surveillance column file to a named list for easy access to values
# see: https://stackoverflow.com/questions/33418288/how-to-convert-a-matrix-to-dictionary-like-a-list
surv_default_text <- as.list(t(surveillance_columns[, "default_text"]))
names(surv_default_text) <- t(surveillance_columns[, "column_name"])

# prepopulate the columns with any default value text
for (name in names(surv_default_text)){
  surv_df[name] <- surv_default_text[name]
}

# -----------------------------------------------
# populate the columns in the surveillance table
# -----------------------------------------------

# ------------------------------------------------------------
# metadata columns populated from metadata input to pipeline
# ------------------------------------------------------------

# Pull in metadata from the metadata dataframe (from excel input)
populate_surveillance_metadata <- function(surv_df, meta_df) {
  
  for (column_name in colnames(meta_df)) {
    if (!(column_name %in% colnames(surv_df))) {
      message (paste0("INFO: metadata column ", column_name, " not present in surveillance table. Will not populate column in surveillance report table."))
      next
    }
    surv_df[column_name] <- meta_df[column_name]
  }
  return (surv_df)
}

# call the function to populate surveillance table metadata
surv_df <- populate_surveillance_metadata(surv_df, metadata_df)

# make a copy of the surveillance table that will have abundance info instead of pos/neg calls
surv_df_abundances <- surv_df
# replace "Negative" values in abundances table with empty text
surv_df_abundances [surv_df_abundances == "Negative"] <- ""

# Populate surveillance tables with observed positive calls and abundances
# TODO: calling this function is slow.  Why?
populate_surveillance_calls <- function(surv_df, 
                                        surv_df_abundances, 
                                        name_type_surveillance_data, 
                                        count_type_surveillance_data, 
                                        observed_surveillance_columns) {
  
  # for each observed surveillance column
  for (i in 1:nrow(observed_surveillance_columns)) {
    
    column_name = observed_surveillance_columns$surveillance_column[i]
    column_type = observed_surveillance_columns$surveillance_column_type[i]
    
    # for each index defined in metadata
    for (each_index in surv_df$Index) {
      
      # count-type column
      if (column_type == "count") {
        
        # pull out data for this index/column
        this_call <- count_type_surveillance_data %>% 
          filter(Index == each_index & surveillance_column == column_name) %>% 
          select(abundance, pos_neg_call)
        
        # should only have a single data point
        num_hits <- nrow(this_call)
        if (num_hits > 1){
          message (paste0("ERROR: more than one hit for index", each_index, " and column: ", this_column))
          quit(status = 1)
        } else if (num_hits == 1) {
          
          # populate pos/neg calls
          surv_df[column_name][surv_df["Index"] == each_index] <- this_call %>% pull(pos_neg_call)
          
          # populate abundances 
          surv_df_abundances[column_name][surv_df["Index"] == each_index] <- this_call %>% pull(abundance)
        }
        # don't need to deal with case of 0 hits
        
      } else if (column_type == "name") {
        # name type column
        
        # pull out data for this index/column
        this_call <- name_type_surveillance_data %>% 
          filter(Index == each_index & surveillance_column == column_name) %>% 
          select(contributing_target_names)
        
        # should only have a single data point
        num_hits <- nrow(this_call)
        if (num_hits > 1){
          message (paste0("ERROR: more than one hit for index", each_index, " and column: ", this_column))
          quit(status = 1)
        } else if (num_hits == 1) {
          
          # populate names that contributed to this hit
          surv_df[column_name][surv_df["Index"] == each_index] <- this_call %>% pull(contributing_target_names)
          
        }
      }
    }
  }
  
  # return both modified tables in a list
  df_returns <- list("surv_df" = surv_df, "surv_df_abundances" = surv_df_abundances)
  return(df_returns)
}

# which surveillance targets were actually observed in the data?
observed_surveillance_targets <- tibble(dataset_with_surv_df %>% 
  filter(!is.na(surveillance_column)) %>% 
  group_by(surveillance_column) %>% 
  summarize())

# pull out surveillance column types, inferred from how they are defined in the targets file
surveillance_column_types <-  targets_to_surv_df %>% 
  select(-target) %>% 
  group_by(surveillance_column, surveillance_column_type) %>% 
  distinct() %>%
  summarize(n_types = n(), .groups="drop")

# check to make sure each column is only defined as one type
more_than_one_type <- surveillance_column_types %>% filter(n_types > 1)
if (nrow(more_than_one_type) > 0) {
  offending_columns <- paste(unique(more_than_one_type %>% pull(surveillance_column)), collapse=" and ")
  message (paste0("ERROR: more than one type (name/column) defined for surveillance column(s): ", offending_columns))
  quit(status = 1)
}

# get rid of no-longer needed column
surveillance_column_types <- surveillance_column_types %>% select(-n_types)

# join in type of surv columns to those that were observed
observed_surveillance_targets <- left_join(observed_surveillance_targets, surveillance_column_types)

# this weird structure because multiple return values from R functions must be in a list
returned_dfs <- populate_surveillance_calls(surv_df, 
                                            surv_df_abundances, 
                                            name_type_surveillance_values, 
                                            count_type_surveillance_values,
                                            observed_surveillance_targets)
surv_df <- returned_dfs$surv_df
surv_df_abundances <- returned_dfs$surv_df_abundances

# rename Acceptable_DNA column from Pos/Neg -> True/False 
surv_df$Acceptable_DNA <- recode(surv_df$Acceptable_DNA, Positive = "TRUE", Negative = "FALSE")

# -------------------
# Write out results
# -------------------

# reorder dataset_df (all_data table) so that Pathogen_Testing_Id is second column.
# see: https://github.com/stenglein-lab/tick_surveillance/issues/29
# need to make sure Pathogen_Testing_ID column exists first: since it is defined
# in metadata input, and this input is flexible.  But if it *is* there, move it.
if ("Pathogen_Testing_ID" %in% colnames(dataset_df)) {
  dataset_df <- dataset_df %>% relocate(Pathogen_Testing_ID, .after = 1)
}

# reorder dataset_df (all_data table) so that observed_sequence column is after mismatch column.
# see: https://github.com/stenglein-lab/tick_surveillance/issues/27
if (all(c("observed_sequence","mismatch") %in% colnames(dataset_df))) {
  dataset_df <- dataset_df %>% relocate(observed_sequence, .after = mismatch)
}

# write all data as csv plain-text file
write.table(dataset_df, paste0(output_dir, "all_data.csv"), quote=F, sep=",", col.names=T, row.names=F)

# create an all_data csv that includes metadata 
# see: https://github.com/stenglein-lab/tick_surveillance/issues/59
dataset_plus_metadata <- left_join(dataset_df, metadata_df, by= c("Index", "Pathogen_Testing_ID", "batch"))
if (nrow(dataset_df) != nrow(dataset_plus_metadata)) {
    message (paste0("ERROR: merging dataset and metadata resulted in an unexpectedly number of rows"))
    quit (status = 1)
}
write.table(dataset_plus_metadata, paste0(output_dir, "all_data_and_metadata.csv"), quote=F, sep=",", col.names=T, row.names=F)

# create excel output
wb <- createWorkbook(paste0(output_dir, "sequencing_report.xlsx"))
modifyBaseFont(wb, fontSize = 11, fontColour = "black", fontName = "Helvetica")

# a generic style for all cells
all_cell_style <- createStyle( 
  border = "TopBottomLeftRight",
  borderColour = getOption("openxlsx.borderColour", "grey"),
  borderStyle = getOption("openxlsx.borderStyle", "thin"),
  halign = "left",
  # numFmt = "0.00",
  wrapText = F
)

# integer number format
integer_cell_style <- createStyle( 
  numFmt = "0"
)

# date format
date_cell_style <- createStyle( 
  numFmt = "yyyy-mm-dd"
)

# text format
text_cell_style <- createStyle( 
  numFmt = "TEXT"
)

# column headers
col_header_style <- createStyle( 
  border = "TopBottomLeftRight",
  borderColour = getOption("openxlsx.borderColour", "grey"),
  borderStyle = getOption("openxlsx.borderStyle", "thin"),
  fontName = "Helvetica",
  fontSize = 11,
  textDecoration = "bold",
  wrapText = F
)


# this function will do some basic styling for worksheets
style_worksheet <- function (wb, sheetname, df) {
  
  # style all cells
  addStyle(wb=wb, sheet=sheetname,
           style=all_cell_style,
           cols = 1:ncol(df),
           rows = 1:nrow(df)+1,
           gridExpand = T
  )
  
  # style date cells
  date_columns <- which(sapply(df, is.Date))
  
  for (date_col in date_columns) {
    # style date cells
    addStyle(wb=wb, sheet=sheetname,
             style=date_cell_style,
             cols = date_col:date_col,
             rows = 1:nrow(df)+1,
             gridExpand = T,
             stack = T
    )
  }
  
  # style column headers
  addStyle(wb=wb, sheet=sheetname,
           style=col_header_style,
           cols = 1:ncol(df),
           rows = 1:1,
           gridExpand = T,
           stack = T
  )
  
  # set column widths to auto
  setColWidths(wb, sheet = sheetname, cols = 1:ncol(df), widths = "auto")
}

# populate the workbook with data

addWorksheet(wb, "TestingResults")
writeData(wb, "TestingResults", surv_df)
style_worksheet(wb, "TestingResults", surv_df)

addWorksheet(wb, "surveillance_counts")
writeData(wb, "surveillance_counts", surv_df_abundances)
style_worksheet(wb, "surveillance_counts", surv_df_abundances)

addWorksheet(wb, "data_by_species")
writeData(wb, "data_by_species", dataset_by_spp)
style_worksheet(wb, "data_by_species", dataset_by_spp)

addWorksheet(wb, "all_data")
writeData(wb, "all_data", dataset_df)
style_worksheet(wb, "all_data", dataset_df)

addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata_df)
style_worksheet(wb, "metadata", metadata_df)

addWorksheet(wb, "targets")
writeData(wb, "targets", targets_df)
style_worksheet(wb, "targets", targets_df)


# add conditional formatting for surveillance worksheet

# light green fill
light_green_fill <- createStyle( 
  bgFill = "#CCFF99"
)

# red fill
red_fill <- createStyle( 
  bgFill = "#FF9999"
)

# Color-coding for Acceptable DNA column 

# identify position of Acceptable DNA column
# This assumes that the column is actually named Acceptable_DNA in the surveillance_columns file
acceptable_DNA_column <- which(colnames(surv_df) == "Acceptable_DNA")
  
conditionalFormatting(wb=wb, sheet="TestingResults", 
                      "colourScale",
                      cols = acceptable_DNA_column:acceptable_DNA_column,
                      rows = 1:nrow(surv_df)+1,
                      style = light_green_fill,
                      rule = "TRUE",
                      type = "contains"
)

# Acceptable DNA column 
conditionalFormatting(wb=wb, sheet="TestingResults", 
                      "colourScale",
                      cols = acceptable_DNA_column:acceptable_DNA_column,
                      rows = 1:nrow(surv_df)+1,
                      style = red_fill,
                      rule = "FALSE",
                      type = "contains"
)

# Pos/Neg calls - red color for positive calls
# this applies to all columns after the Acceptable_DNA column
conditionalFormatting(wb=wb, sheet="TestingResults", 
                      "colourScale",
                      cols = acceptable_DNA_column:ncol(surv_df),
                      rows = 1:nrow(surv_df)+1,
                      style = red_fill,
                      rule = "Positive",
                      type = "contains"
                      
)

# add integer style to surveillance counts
addStyle(wb, "surveillance_counts", integer_cell_style, rows=1:nrow(surv_df_abundances)+1, cols=acceptable_DNA_column+1:ncol(surv_df_abundances)+1, gridExpand =T, stack = T)

# write out the workbook
saveWorkbook(wb, paste0(output_dir, "sequencing_report.xlsx"), overwrite = TRUE)

# write out version info into versions.yml
# TODO
writeLines("", "versions.yml")

