/*
 This process modifies the unassigned_seuqnece report to include sample metadata 
*/
process FILTER_UNASSIGNED_SEQUENCES {
  publishDir "${params.outdir}", mode: 'link'
  label 'process_low'

  // if using conda
  conda "$baseDir/conda/R_conda_environment.yaml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.4.1"
  }

  input:
  path(surveillance_report) 
  path(unassigned_seq_report)
  val(filter)
  path(R_lib_dir_input)

  output:
  path('sequencing_report.xlsx')               , emit: sequences_report_filter

  script:
  // only use R lib dir for singularity
  def r_lib_dir = workflow.containerEngine == 'singularity' ? "${R_lib_dir_input}" : "NA"
  """
  filter_unassigned_seq_report.R $r_lib_dir $surveillance_report $unassigned_seq_report $filter
  """
}

