/*
  This process runs dada2 on input files
*/
process DADA2 {      
  publishDir "${params.dada_outdir}", mode: 'link'
  publishDir "${params.QC_and_summary_stats_dir}", pattern: 'dada_read_clean_all.csv', mode: 'copy'
  tag        "all"

  label 'process_high'
  label 'error_retry'

  // if using conda
  conda "${moduleDir}/environment.yml"                                          

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0"
  } else {
      container "quay.io/biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0"
  }

  input:
  path (all_trimmed_reads)
  val  (maxN)
  val  (maxEE)
  val  (truncQ)
  val  (trimRight)
  val  (min_reads)
  val  (min_overlap)
  val  (max_mismatch)

  output:
  path("dada_seqtab.txt")                    , emit: seqtab
  path("dada_filtered")                      , emit: dada_filtered_dir
  path "versions.yml"                        , emit: versions 
  path ("dada_read_clean_all.csv")           , emit: dada_read_tracking_all                                        

  script:
  """
  # Run dada2 using trimmed fastq as input and create a tabular output of results
  run_dada_on_trimmed.R . $maxN $maxEE $truncQ $trimRight $min_reads $min_overlap $max_mismatch
  """
}

process TIDY_DADA_OUTPUT {      
  publishDir "${params.dada_outdir}", mode: 'link'
  publishDir "${params.QC_and_summary_stats_dir}", pattern: 'dada_read_clean_summary.csv', mode: 'copy'
  tag        "all"

  label 'process_low'

  // if using conda
  conda "$baseDir/conda/R_conda_environment.yaml"                               

  // if using singularity 
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.4.1"
  }

  input:
  path(dada_seqtab)
  path(dada_read_tracking_all)  

  output:
  path "observed_sequences.fasta"      , emit: observed_sequences
  path "sequence_abundance_table.tsv"  , emit: sequence_abundance_table
  path "dada_read_clean_summary.csv"   , emit: dada_read_tracking_summary
  path "dada_read_clean_all_with_pct_pass.csv"   
  path "versions.yml"                  , emit: versions                                         

  script:
  """
  # This R script creates an output named observed_sequences.fasta containing all of the
  # unique sequences observed in the amplicon dataset
  # and a sequence_abundance_table.tsv, which lists the abundances of these
  # sequences in each dataset
  tidy_dada_output.R $dada_seqtab $dada_read_tracking_all
  """
}

