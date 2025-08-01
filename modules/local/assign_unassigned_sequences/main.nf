/*                                                                              
   This process parses the blast output from unassigned sequences and           
   populates a spreadsheet with the info.                                       
*/                                                                              
                                                                                
process ASSIGN_UNASSIGNED_SEQUENCES {                                           
  label 'process_low'                                                           
  tag "all"                                                                     
                                                                                
  // if using conda                                                             
  conda "$baseDir/conda/R_conda_environment.yaml"                               
                                                                                
  // if using singularity                                                       
  if (workflow.containerEngine == 'singularity'){                               
      container "docker://rocker/tidyverse:4.4.1"
  }                                                                             
                                                                                
  input:                                                                        
  path(blast_output)                                                            
  path(unassigned_sequences)                                                    
  path(R_lib_dir_input)                                                         
                                                                                
  output:                                                                       
  path("non_reference_sequence_assignments.xlsx")  , emit: unassigned_sequences_report
  path "versions.yml"                              , emit: versions             
                                                                                
  script:                                                                       
  // only use R lib dir for singularity                                         
  def r_lib_dir = workflow.containerEngine == 'singularity' ? "${R_lib_dir_input}" : "NA"
  """                                                                           
  assign_non_ref_seqs.R $r_lib_dir $unassigned_sequences $blast_output          
  """                                                                           
}
