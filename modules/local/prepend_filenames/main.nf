/*
  This process prepends main output files with a prefix (by default the date, can be overridden)
*/

process PREPEND_OUTPUT_FILENAMES {
  label 'process_low'

  input:
  path(output_file) 

  output:
  path ("${params.output_prefix}${output_file}"), emit: prepended_file

  script:
  """
  # -L forces the resulting link to be a hard link 
  # even if it was created via a symlink
  ln -L $output_file ${params.output_prefix}${output_file}
  """
}

