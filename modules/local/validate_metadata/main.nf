/*
  This process validates that the metadata file is in the approprate
  format and contains the appropriate information.  Logic in the R script.
*/
process VALIDATE_METADATA {
  tag        "${metadata}"
  publishDir "${params.outdir}", mode: "copy"

  label 'process_low'

  // if using conda
  conda "$baseDir/conda/R_conda_environment.yaml"                               

  // if using singularity 
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.4.1"
  }

  input:
  path (metadata)         // the metadata-containing excel spreadsheet
  path (sample_ids)       // all the sample IDS inferred from the fastq file names

  output:
  path(metadata)       , emit: validated_metadata
  path(sample_ids)     , emit: sample_ids
  path "versions.yml"  , emit: versions                                         


  script:
  """
    validate_metadata.R $metadata $sample_ids
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
  """
}

