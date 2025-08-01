/*
 This process modifies the unassigned_seuqnece report to include sample metadata 
*/
process ORG_UNASSIGNED_SEQUENCES {
  publishDir "${params.outdir}", mode: 'link'
  label 'process_low'

  // if using conda
  conda "$baseDir/conda/R_conda_environment.yaml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.4.1"
  }

  input:
  path(sequence_abundance_table) 
  path(metadata)
  path(non_ref_assignments)
  path(R_lib_dir_input)

  output:
  path('non_reference_sequence_assignments.xlsx')               , emit: org_unassigned_sequences_report

  script:
  // only use R lib dir for singularity
  def r_lib_dir = workflow.containerEngine == 'singularity' ? "${R_lib_dir_input}" : "NA"
  """
  org_unassigned_seq_report.R $r_lib_dir $sequence_abundance_table $metadata $non_ref_assignments
  """
}

