
/*
   This sets up a python virtual environment (venv) containing the packages
   needed by python scripts in this pipeline.
   
   see: https://docs.python.org/3/library/venv.html

*/
process SETUP_PYTHON_ENVIRONMENT {
  label      'process_low'
  tag        "${venv_input_path}"
  publishDir "${params.outdir}"

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10.4"
  } else {
      container "quay.io/biocontainers/python:3.10.4"
  }

  when:
  params.make_trees 

  input:
  val(venv_input_path)
  path(requirements_path)
   
  output:
  // this output will be a signal that venv setup is complete
  path (venv_input_path) , emit: venv_path
  path "versions.yml"  , emit: versions

  
  script:

  if (workflow.containerEngine == 'singularity') {
  """
    python -m venv --clear $venv_input_path
    source ${venv_input_path}/bin/activate
    # install modules needed for tree-building scripts
    # see https://github.com/stenglein-lab/tick_surveillance/issues/76
    # see https://github.com/stenglein-lab/tick_surveillance/issues/77
    pip install -r ${requirements_path}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
  """
  } else {
  // only need to make a python venv when using singularity
  """
    # just make an empty directory to pass to processes expecting it.
    mkdir $venv_input_path

    touch versions.yml
  """
  }
}


