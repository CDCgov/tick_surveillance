/*
 This process aligns all of the unique sequences reported by dada2
 to the set of expected reference sequences using blastn.
*/
process COMPARE_OBSERVED_SEQS {
  publishDir "${params.blast_outdir}", mode: 'link', pattern: "*bn_refseq"
  tag "all"

  label 'process_low'

  // if using conda
  conda "${moduleDir}/environment.yml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_2"
  } else {
      container "quay.io/biocontainers/blast:2.16.0--hc155240_2"
  }

  input:
  path (observed_sequences) 
  path (sequence_abundance_table) 
  path (refseq_blast_db)
  path (blast_db_dir) 

  output:
  path(sequence_abundance_table)            , emit: sequence_abundance_table
  path("${observed_sequences}.bn_refseq")   , emit: blast_output
  path "versions.yml"                       , emit: versions                                         

  script:
  // this is similar to the default blastn output except gaps replaces gapopens, because seems more useful!
  // also include other useful columns like query length, subject length, etc
  def blastn_columns = "qaccver saccver pident length mismatch gaps qlen slen bitscore"
  """
  # this finds the blast DB files in the pwd
  DB=`find -L ./ -name "*.ndb" | sed 's/\\.ndb\$//'`                          
  blastn -db \$DB -task blastn -evalue ${params.max_blast_refseq_evalue} -query $observed_sequences -outfmt "6 $blastn_columns" -out ${observed_sequences}.bn_refseq.no_header
  # prepend blast output with the column names so we don't have to manually name them later
  echo $blastn_columns > blast_header.no_perl
  echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header
  cat blast_header ${observed_sequences}.bn_refseq.no_header > ${observed_sequences}.bn_refseq

  cat <<-END_VERSIONS > versions.yml                                          
  "${task.process}":                                                          
      blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')      
  END_VERSIONS  
  """
}


/*
  This process takes the blast alignment information from the above process
  and decides whether those sequence are sufficiently similar to the reference
  sequences to be assigned to them.
*/
process ASSIGN_OBSERVED_SEQS {
  label 'process_low'
  tag "all"

  // if using conda
  conda "$baseDir/conda/R_conda_environment.yaml"

  // if using singularity
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.4.1"
  }

  input:
  path(metadata)
  path(abundance_table)
  path(blast_output)
  path(R_lib_dir_input) 
  path(surveillance_columns_file) 
  path(targets_file) 

  output:
  path("unassigned_sequences.fasta") , emit: unassigned_sequences
  path("sequencing_report.xlsx")     , emit: surveillance_report
  path("all_data*.csv")              , emit: all_data_csv
  path("*.txt")                      , emit: txt
  path("*.pdf")                      , emit: pdf
  path "versions.yml"                , emit: versions                                         

  // output channels for tree-building process
  // path("sequencing_report.xlsx") into report_tree_ch

  script:
  // only use R lib dir for singularity
  def r_lib_dir = workflow.containerEngine == 'singularity' ? "${R_lib_dir_input}" : "NA"
  """
  assign_observed_seqs_to_ref_seqs.R $r_lib_dir $abundance_table $blast_output $metadata ${targets_file} ${surveillance_columns_file} ${params.min_reads_for_positive_surveillance_call} 
  """
}

