/*
  This workflow ensures that the local nt blast database and blast taxonomy 
  database are valid.  These are used to taxonomically assign unassigned sequences.
*/
workflow SETUP_BLAST_DB_AND_TAX {

  take: 
  // no input 

  main:
  
  ch_versions = Channel.empty()                           

  //  a channel containing the path to an (optional) local copy of the nt blast database.
  blast_db_ch = Channel.empty()
  
  if (params.blast_unassigned_sequences && params.local_nt_database_dir) {
    blast_db_ch = Channel.fromPath( params.local_nt_database_dir )
  } else if (params.blast_unassigned_sequences && params.remote_blast_nt) {
    // if remote BLASTing, don't need to point to a directory containing a copy
    // of the local nt BLAST database, but need to provide a non-empty channel
    // so process check_local_blast_database will run.  
    // This is based on the pattern described here: 
    // https://nextflow-io.github.io/patterns/optional-input/
    blast_db_ch = file( "not_a_real_db_path_but_keeps_channel_from_being_empty" )
  }

  CHECK_LOCAL_BLAST_DATABASE(blast_db_ch)
  ch_versions = ch_versions.mix(CHECK_LOCAL_BLAST_DATABASE.out.versions)

  /*
    Create a channel containing the path to an (optional) local copy of the 
    NCBI blast/taxonomy db to make BLAST taxonomically aware
   */
  blast_tax_ch = Channel.empty()
  
  if (params.blast_unassigned_sequences) {
    // if this path was provided as a parameter, then create a channel
    // from this path and set a boolean to true to indicate it's an existing
    // directory
    if (params.blast_tax_dir) {
       blast_tax_ch = Channel.fromPath( params.blast_tax_dir )
                             .map { path -> [ path , true ] }  
    } else {
       // if this path was *not* provided as a parameter, then create a channel
       // from a bogus path "blast_tax_dir" and set a boolean to false 
       // to indicate it *doesn't* refer to an existing
       blast_tax_ch = Channel.fromPath( "blast_tax_dir" )
                             .map { path -> [ path , false ] }  
    }
  } 

  CHECK_BLAST_TAX(blast_tax_ch)

  emit:
  blast_db_dir  = CHECK_LOCAL_BLAST_DATABASE.out.checked_db_dir
  blast_tax_dir = CHECK_BLAST_TAX.out.checked_blast_tax_dir
  versions      = ch_versions
}

/*
  This process confirms that the local NT database is valid 
  (if applicable: if going to blast unassigned sequences and 
   if not doing a remote blast)
 */

process CHECK_LOCAL_BLAST_DATABASE {
  label 'process_low'
  tag "all"

  // if using conda
  conda "${moduleDir}/blast_environment.yml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  when: 
  params.blast_unassigned_sequences 

  input:
  path(local_nt_database_dir) 

  output:
  path(local_nt_database_dir)     , emit: checked_db_dir 
  path "versions.yml"             , emit: versions                                         

  script:                                                                       
  def local_nt_database = "${local_nt_database_dir}/${params.local_nt_database_name}"

  params.remote_blast_nt ? 
  """
    echo "Don't need to check local BLAST database when running with -remote option"
    rm $local_nt_database_dir
    touch $local_nt_database_dir
    
    # create empty versions.yml
    touch versions.yml
  """ : 
  """
    # check for expected .nal file: if not present, output a helpful warning message
    if [ ! -f "${local_nt_database_dir}/${params.local_nt_database_name}.nal" ]                                 
    then                                                                          
      echo "ERROR: it does not appear that ${local_nt_database} is a valid BLAST database."
    fi   
  
    # check validity of database with blastdbcmd.  If not valid, this will error 
    # and stop pipeline.
    blastdbcmd -db ${local_nt_database} -info 

    cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        blast: \$(blastdbcmd -version 2>&1 | sed 's/^.*blastdbcmd: //; s/ .*\$//')      
    END_VERSIONS       
  """
}

/* 
  This process confirms that the NCBI blast taxdb is installed locally.

  This db is required for BLAST to produce taxonomically aware results.

  This process will check the db if it is defined (that is, if 
  params.blast_tax_dir is defined), or will download the db if 
  params.blast_tax_dir is not defined
 */
process CHECK_BLAST_TAX {
  label 'process_low'
  tag "all"

  // if using conda 
  conda "${moduleDir}/curl_environment.yml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/curl:7.80.0"
  } else {
      container "quay.io/biocontainers/curl:7.80.0"
  }

  when: 
  params.blast_unassigned_sequences

  input:
  tuple path(blast_tax_dir), val(existing_db) 

  output:
  path (blast_tax_dir), emit: checked_blast_tax_dir 

  script:
  // if a local blast_tax_dir is specified, check that it contains the expected files
  if (existing_db) {
  """
    # check that the directory exists
    if [ ! -d "${blast_tax_dir}" ] ; then
      echo "ERROR: BLAST taxonomy directory ${blast_tax_dir} (--blast_tax_dir) does not exist."
      exit 1
    fi 
    # check that appropriate files exist
    if [ ! -f "${blast_tax_dir}/taxdb.btd" ] ; then
      echo "ERROR: required BLAST taxonomy file taxdb.btd not in directory ${blast_tax_dir} (--blast_tax_dir)."
      exit 1
    fi 
    if [ ! -f "${blast_tax_dir}/taxdb.bti" ] ; then
      echo "ERROR: required BLAST taxonomy file taxdb.bti not in directory ${blast_tax_dir} (--blast_tax_dir)."
      exit 1
    fi 
  
  """ 
  } else {
  // if tax db doesn't already exist : download the necessary files and keep track of directory 
  """
    # make a new local directory to contain the files
    # first remove broken link to non-existent file
    rm $blast_tax_dir
    mkdir $blast_tax_dir
    # download taxdb files
    curl -OL https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
    # unpack archive
    tar xvf taxdb.tar.gz
    # move files to blast_tax_dir
    mv taxdb.??? $blast_tax_dir
    # move sqlite3 file needed for restricting by taxid
    mv taxonomy4blast.sqlite3 $blast_tax_dir
    # get rid of archive
    rm taxdb.tar.gz
  """
}
}


