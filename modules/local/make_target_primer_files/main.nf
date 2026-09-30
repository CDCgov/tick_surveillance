params.make_targets = true

process MAKE_TARGET_PRIMER_FILES {
  label      'process_low'
  tag        "${venv_path}"
  publishDir "${params.outdir}", mode: 'link'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10.4"
  } else {
      container "quay.io/biocontainers/python:3.10.4"
  }


  input:
  path (venv_path)
  path(all_ref)
  val(primer_mix)
  val(study_type)
   
  output:
  
  path ("targets.tsv"), emit: targets
  path ("primers.tsv"), emit: primers 

  
  script:
    // only need to activate the venv for singularity
  def activate_venv_command = workflow.containerEngine == 'singularity' ? "source ${venv_path}/bin/activate" : ""
  """
  $activate_venv_command
  python ${params.script_dir}/make_targets.py $all_ref $primer_mix $study_type
  """
}
