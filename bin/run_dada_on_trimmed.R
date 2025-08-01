#!/usr/bin/env Rscript 
#
# This script runs dada2 on a set of trimmed fastq files 
# 
# See: 
# https://benjjneb.github.io/dada2/tutorial.html
#

library(dada2)

#
# This code block sets up input arguments to either come from the command line
# (if not running in interactive mode) or to use expected default values 
# (if running in interactive mode, i.e. RStudio).  
#
# This is useful in case you want to modify this script
#
if (!interactive()) {
  # if running from Rscript
  args               = commandArgs(trailingOnly=TRUE)
  trimmed_path       = args[1]
  outdir             = "./"
  input_maxN         = as.numeric(args[2])
  input_maxEE        = as.numeric(args[3])
  input_truncQ       = as.numeric(args[4])
  input_trimRight    = as.numeric(args[5])
  input_min_reads    = as.numeric(args[6])
  input_min_overlap  = as.numeric(args[7])
  input_max_mismatch = as.numeric(args[8])
} else {
  # if running via RStudio (for development or troubleshooting)
  trimmed_path     = "../test/results/trimmed_fastq"
  outdir           = "../test/results/"
  input_maxN       = 0
  input_maxEE      = 2
  input_truncQ     = 2
  input_trimRight  = 0
  input_min_reads  = 10
  input_min_overlap  = 12
  input_max_mismatch = 0
}

# get lists of fastq files in trimmed directory 
fnFs <- sort(list.files(trimmed_path, pattern="fully_trimmed.R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(trimmed_path, pattern="fully_trimmed.R2.fastq.gz", full.names = TRUE))

# extract sample names from fastq file names
# so that sample names will match execpted values provided to sample sheet 
# when demultiplexing

# remove _R[12]_trimmed.fastq.gz from file names to make better sample names
sample.names <- gsub(".fully_trimmed.R[12].fastq.gz", "", basename(fnFs))

# keep a map of sample names -> input filenames
sample_name_file_name_map <- data.frame(sample_id = sample.names, filename = basename(fnFs))

# This will place dada-filtered files in dada_filtered/ subdirectory of the pwd
filtFs <- file.path("dada_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("dada_filtered", paste0(sample.names, "_R_filt.fastq.gz"))


names(filtFs) <- sample.names
names(filtRs) <- sample.names

#
# Run dada2 filtering and trimming

# filterAndTrim options:
#
# truncQ 	
# (Optional). Default 2. Truncate reads at the first instance of a quality score less 
# than or equal to truncQ.

# truncLen 	
# (Optional). Default 0 (no truncation). Truncate reads after truncLen bases. 
# Reads shorter than this are discarded.

# maxEE
# (Optional). Default Inf (no EE filtering). After truncation, reads with higher than 
# maxEE "expected errors" will be discarded. Expected errors are calculated from the 
# nominal definition of the quality score: EE = sum(10^(-Q/10))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     maxN        = input_maxN, 
		     maxEE       = c(input_maxEE, input_maxEE), 
		     truncQ      = input_truncQ, 
		     rm.phix     = TRUE,
		     trimRight   = input_trimRight,
                     compress    = TRUE, 
		     multithread = TRUE) # On Windows set multithread=FALSE


# only keep datasets with >0 reads to avoid dada2 throwing an error because of empty datasets
# empty datasets will not contribute to error learning and will have no assigned ASVs, as expected.
# see: https://github.com/benjjneb/dada2/issues/469
keep <- out[,"reads.out"] > input_min_reads 
filtFs <- filtFs[keep]
filtRs <- filtRs[keep]


# also only keep file pairs for which the filtered files exist
# see: https://github.com/benjjneb/dada2/issues/375
exists <- file.exists(filtFs) & file.exists(filtRs)
# use exists logical vector to subset filtered file names
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]


# Estimate error rates: actual error rate as a function of quality score for each type   
# of base substitution.
#
# From the dada2 documentation (https://benjjneb.github.io/dada2/tutorial.html): 
#
# "The DADA2 algorithm makes use of a parametric error model (err) and every amplicon 
# dataset has a different set of error rates. The learnErrors method learns this error 
# model from the data, by alternating estimation of the error rates and inference of 
# sample composition until they converge on a jointly consistent solution."
#
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Could save this as output...
# plotErrors(errF, nominalQ=TRUE)
# plotErrors(errR, nominalQ=TRUE)

# save the error profile for bases with a Q-score of 35
# q35_errors <- errF$err_out[,36]
# q35_errors
# write.table(q35_errors, "q35_errors.txt", sep="\t")

# Run dada function to identify ASVs
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, 
		      dadaRs, filtRs, 
                      minOverlap  = input_min_overlap,
                      maxMismatch = input_max_mismatch,
		      verbose     = TRUE)

# create a table of dada2 results
seqtab <- makeSequenceTable(mergers)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#
# Handle empty datasets.  These should be reported as datasets with abundances of 0 for all sequences.
#
# These could correspond to datasets that were truly empty after sequencing or had no reads after filtering.
#

# Add rows with all 0s for datasets that had no reads after filtering 
# this corresponds to the FALSE values in the exists vector
empty_sample_names <- sample.names[!exists]
# create a matrix with all 0s
new_rows_with_zeros <- matrix(0L, nrow = length(empty_sample_names), ncol = ncol(seqtab.nochim)) 
# name them according to their original names
row.names(new_rows_with_zeros) <- empty_sample_names
# concatenate this new all-zeros table onto the end of the data seqtab
seqtab.nochim <- rbind(seqtab.nochim, new_rows_with_zeros)

# convert to a data frame
t <- as.data.frame(seqtab.nochim)

# make a new column based on row names
t$dataset <- rownames(t)

# write the output to a table
write.table(t, paste0(outdir, "dada_seqtab.txt"), sep="\t", col.names=T, quote=F)

# write out version info into versions.yml
writeLines(c("DADA2:", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")),paste0("    ShortRead: ", packageVersion("ShortRead")) ), "versions.yml")

# turn out object into data frame with sample ID
out_df <- as.data.frame(out)
out_df$filename <- rownames(out)
# switch from filename to sample ID
# merge in sample IDs based on file names
out_df <- merge(out_df, sample_name_file_name_map)
out_df <- out_df[ , c("sample_id", "reads.in", "reads.out")]
colnames(out_df) <- c("sample_id", "input", "filtered")

# Track reads and create output file (https://benjjneb.github.io/dada2/tutorial.html)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
getN <- function(x) sum(getUniques(x))

# this function returns a data frame instead of a named
getN_to_df <- function(x, column_name) {
   Ns <- sapply(x, getN)
   new_df <- data.frame(sample_id = names(Ns), temp_new_name = Ns)
   colnames(new_df) <- c("sample_id", column_name)
   new_df
}

# turn these lists into data frames with reasonable names
dadaFs_N  <- getN_to_df(dadaFs, "denoisedF")
dadaRs_N  <- getN_to_df(dadaFs, "denoisedR")
mergers_N <- getN_to_df(mergers, "merged")

# non-chimeric read counts in seqtab
nochim_N <- data.frame(sample_id = rownames(seqtab.nochim), nonchim = rowSums(seqtab.nochim))

# combine all columns using sample_id as key
track <- merge(out_df, dadaFs_N, all.x = T)
track <- merge(track, dadaRs_N, all.x = T)
track <- merge(track, mergers_N, all.x = T)
track <- merge(track, nochim_N, all.x = T)

# convert NA values to 0s
track[is.na(track)] <- 0

# write read tracking info to csv file
write.csv(track, paste0(outdir, "dada_read_clean_all.csv"), row.names = F)
