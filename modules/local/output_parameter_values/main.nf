/*
  This (sub)workflow outputs metadata information about the pipeline:
    - parameter values specified for this run
    - information from the workflow namespace
    - information from the nextflow namespace
    for more information, see: https://docs.seqera.io/nextflow/reference/stdlib-namespaces
 */

// helper function that avoids possible issue with 
// undefined object properties
// returns object property value or '<unavailable>'
def safeGet(obj, prop) {
  try {
   obj[prop]
  }
  catch (ignored) {
    '<unavailable>'
  }
}


/*
  This process prepends main output files with a prefix (by default the date, can be overridden)
*/

process OUTPUT_PARAMETER_VALUES {
  label 'process_low'

  output:
  path ("parameter_values.txt"), emit: parameter_values
  path ("workflow_values.txt"),  emit: workflow_values
  path ("nextflow_values.txt"),  emit: nextflow_values

  script:

  // output parameter values to a file in format:
  // parameter_name<tab>parameter_value                                            
  // sort alphabetically (ignoring case)
  def params_txt = params
    .toSorted { a, b -> a.key.toLowerCase() <=> b.key.toLowerCase() }
    .collect { key, value -> "${key}\t${value}" }                                                          
    .join('\n')                                                                    

  // collect information about workflow namespace (e.g. launchDir, projectDir)
  def workflow_properties = 
   ['commandLine', 
    'commitId',
    'configFiles',
    'containerEngine',
    'homeDir',
    'launchDir',
    'outputDir',
    'profile',
    'projectDir',
    'repository',
    'resume',
    'revision',
    'runName',
    'scriptFile',
    'scriptId',
    'scriptName',
    'sessionId',
    'start',
    'stubRun',
    'userName',
    'workDir'
      ].sort { a, b -> a.compareToIgnoreCase(b) }
  
  def workflow_txt = workflow_properties
   .collect { p -> "workflow.${p}\t${safeGet(workflow,p)}" }
   .join('\n')

  // collect information about nextflow namespace (e.g. launchDir, projectDir)
  def nextflow_properties = 
   ['build', 
    'timestamp',
    'version'
      ].sort { a, b -> a.compareToIgnoreCase(b) }

  def nextflow_txt = nextflow_properties
   .collect { p -> "nextflow.${p}\t${safeGet(nextflow,p)}" }
   .join('\n')

  """
  echo "${params_txt}"   > parameter_values.txt
  echo "${workflow_txt}" > workflow_values.txt
  echo "${nextflow_txt}" > nextflow_values.txt
  """
}

