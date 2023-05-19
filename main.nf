#!/usr/bin/env nextflow

/*
    CDC Tick Surveillance Amplicon Sequencing Analysis Pipeline

    December 13, 2021

    Mark Stenglein
*/


// TODO: command line options and parameter checking

// running the pipeline with these params set will output usage info
// e.g. nextflow run main.nf --h
params.help = false
params.h = false

WorkflowMain.initialise(workflow, params, log)                                  


/*
   Pipeline usage output message.

   This will provide usage information in case of a parameter failure
*/
def usageMessage() {

  log.info """

        """
}

/* 
  Pipeline informational output message
*/
// This 
def infoMessage() {

  // tie this versioning to github tags/releases

  log.info """

    Pipeline version:             ${params.pipeline_version}
    Reference sequences version:  ${params.refseq_version}
    Primers version:              ${params.primers_version}

  """
}

/*
  Print usage information if necessary and exit.
*/
if (params.help || params.h) {
    usageMessage()
    exit 0
}

/* 
  Check input parameters 
*/
def check_params_and_input () {

  // This list includes a list of files or paths that are required 
  // to exist.  Check that they exist and fail if not.  
  checkPathParamList = [
    params.fastq_dir, 
    params.script_dir, 
    params.targets, 
    params.metadata, 
    params.primers 
  ]
  // log.info("Checking for required input paths and files...")
  for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

  check_blast_parameters()
  
}


/*
  Check that BLAST (of unassigned sequences) parameters are set in a way that makes sense
 */
def check_blast_parameters() {
  if (params.blast_unassigned_sequences && !(params.local_nt_database_dir || params.remote_blast_nt)) {
    log.error("Error: blast_unassigned_sequences is true but neither local_nt_database_dir nor remote_blast_nt is specified.")
    System.exit(1)
  }
  if (params.blast_unassigned_sequences && (params.local_nt_database_dir && params.remote_blast_nt)) {
    log.error(
    """
    Error: blast_unassigned_sequences is true and both local_nt_database_dir and remote_blast_nt are specified.
           Only one of local_nt_database_dir and remote_blast_nt must be specified.
    """)
    System.exit(1)
  }
}


/* 
  Check parameters and input
*/
check_params_and_input()


/* 
  Check for existence of metadata file
*/
Channel
    .fromPath("${params.metadata}", 
                   checkIfExists: true) 
    .set {post_metadata_check_ch}

                                                                                
/* 
  Read in the targets tsv-format file that describes the expected target sequences. 
*/
Channel
    .fromPath(params.targets, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { targets_ch }
                                                                                


/* 
  generate one fasta for each reference sequence
*/
process generate_refseq_fastas {                                                      

  label 'process_low'

  input:                                                                        
  val targets from targets_ch

  output:                                                                        
  path ("${targets.ref_sequence_name}.fasta") into refseq_fastas_ch

  // makes fasta formatted records for each target
  script:                                                                       
  """
  printf ">%s\n%s\n" $targets.ref_sequence_name $targets.sequence > ${targets.ref_sequence_name}.fasta
  """
}

/* 
  concatenate individual fasta into one big reference sequence file 
*/
process combine_refseq_fasta {                                                      
  publishDir "${params.refseq_fasta_dir}", mode: 'link'

  label 'process_low'

  input:                                                                        
  path (individual_fastas) from refseq_fastas_ch.collect()

  output:                                                                        
  path ("reference_sequences.fasta") into refseq_fasta_ch

  // makes fasta formatted records for each targets 
  script:                                                                       
  """
  cat $individual_fastas > reference_sequences.fasta
  """
}

/*
   Setup indexes and dictionaries needed by downstream processes.

   Only do this once at beginning.
*/
process setup_indexes {

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  input:
  path (refseq_fasta) from refseq_fasta_ch

  output:
  // this output will be a signal that indexes are setup and processes that
  // need them can proceed
  val ("indexes_complete") into post_index_setup_ch
  
  // blast database
  // have to pass paths to both the blast db prefix/name and paths to the additional 
  // files composing the blast db so they will be available to downstream process
  tuple path ("$refseq_fasta"), path ("${refseq_fasta}.*") into refseq_blast_db_ch


  script:
  """
  # Blast database of the reference sequences
  makeblastdb -dbtype nucl -in ${refseq_fasta} -out "${refseq_fasta}"
  """
}

/*
   This sets up a python virtual environment (venv) containing the packages
   needed by python scripts in this pipeline. 
   
   see: https://docs.python.org/3/library/venv.html

*/
process setup_python_venv {
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10.4" 
  } else {
      container "quay.io/biocontainers/python:3.10.4"
  }

  when:
  params.make_trees 


  output:
  // this output will be a signal that venv setup is complete
  val ("venv_complete") into post_venv_setup_ch
  
  script:

  if (workflow.containerEngine == 'singularity') {
  """
    python -m venv ${params.python_venv_path}
    source ${params.python_venv_path}/bin/activate
    # install modules needed for tree-building scripts
    # see https://github.com/stenglein-lab/tick_surveillance/issues/76
    # see https://github.com/stenglein-lab/tick_surveillance/issues/77
    pip install -r ${params.python_requirements}
  """
  } else {
  """
    echo "only need to make a venv when using singularity"
  """
  }
}


/*
   This installs a couple R packages that are not included in the
   base tidyverse singularity image we are using.
*/
process setup_R_dependencies {
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  output:
  // this output will be a signal that venv setup is complete
  val ("r_dependencies_OK") into post_r_dep_setup_ch
  
  script:

  if (workflow.containerEngine == 'singularity') {
  """
     mkdir -p ${params.R_lib_dir}
     Rscript ${params.script_dir}/install_R_packages.R ${params.R_tar_dir} ${params.R_lib_dir}
  """
  } else {
  """
    echo "setup not necessary for conda environment"
  """
  }
}


/*
 These fastq files represent the main input to this workflow
 
 Expecting files with _R1 or _R2 in their names corresponding to paired-end reads
*/


Channel
    .fromFilePairs("${params.fastq_dir}/${params.fastq_pattern}", 
                   size: 2, 
                   // checkIfExists: true, 
                   maxDepth: 1)
    // this map step first strips off any _L### text from the end of the sample ID
    // and then any _S## text from the end 
    // this is text added onto IDs from samplesheet that are added by Illumina bcl2fastq
    .map { untrimmed_sample_id, fastq  ->
           def sample_id = untrimmed_sample_id.replaceFirst( /_L\d{3}$/, "")
           sample_id = sample_id.replaceFirst( /_S\d+$/, "")
           [sample_id, fastq] }
    .into {samples_ch_qc; samples_ch_trim; sample_ids_ch}


/* 

 This channel generates the primer sequences that were 
 used to amplify surveillance targets.  These primer sequences will be trimmed 
 off of read pairs and used to identify legitimate PCR products, which will
 contain an expected F/R primer pair at the ends. 

*/
Channel
    .fromPath(params.primers, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { primers_ch }


/*
   this combinatorially mixes the fastq file pairs with the primer sequences
   to be trimmed, creating a new channel with all possible combinations of
   input fastq and primer pair
*/

primers_and_samples = primers_ch.combine(samples_ch_trim)

/*                                                                              
   A function to return the complement of a DNA sequence                        
   combine with reverse() to get the reverse complement                         
                                                                                
   from http://groovyconsole.appspot.com/script/29005                           
*/                                                                              
String.metaClass.complement = {                                                 
  def complements = [ A:'T', T:'A', U:'A', G:'C',                               
                      C:'G', Y:'R', R:'Y', S:'S',                                 
                      W:'W', K:'M', M:'K', B:'V',                                 
                      D:'H', H:'D', V:'B', N:'N' ]                                
  delegate.toUpperCase().collect { base ->                                      
    complements[ base ] ?: 'X'                                                  
  }.join()                                                                      
}                                                                               

/*
 Output sample IDs into a file
*/
process output_sample_ids {

  label 'process_low'

  input:
  tuple val(sample_id), path(initial_fastq) from sample_ids_ch

  output:
  val(sample_id) into tabulate_sample_ids_ch

  script:
  """
  # this doesn't actually do anything but the process wants a script block
  echo $sample_id
  """
}

/*
 Output sample IDs into a file
*/
tabulate_sample_ids_ch
  .collectFile(name: 'sample_ids.txt', newLine: true)
  .set{ sample_ids_file_ch }

/*
  This process validates that the metadata file is in the approprate
  format and contains the appropriate information.  Logic in the R script.
*/
process validate_metadata_file {
  publishDir "${params.log_outdir}", mode: 'link', pattern: "*.txt"

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     
  
  input:
  path(metadata) from post_metadata_check_ch
  path(sample_ids) from sample_ids_file_ch

  output:
  path(metadata) into validated_metadata_ch
  path(sample_ids)

  script:
  """
  Rscript ${params.script_dir}/validate_metadata.R ${params.script_dir} $metadata $sample_ids 
  """
}

/*
 Run fastqc on input fastq 
*/
process initial_qc {
  label 'many_forks'
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
  } else {
      container "quay.io/biocontainers/fastqc:0.11.9--0"
  }

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  path("*.zip") into post_initial_qc_ch

  script:
  """
  fastqc $initial_fastq 
  """
}

/*
 Use multiqc to merge initial fastqc reports

 TODO: will this report old, accumulated fastqc reports if the pipeline is re-run without cleaning up the work directory?
*/
process initial_multiqc {

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
  }

  input:
  path(all_zip) from post_initial_qc_ch.collect()

  output: 
  path("initial_qc_report.html") into initial_multiqc_output_ch

  script:
  """
  multiqc --interactive -n "initial_qc_report.html" -m fastqc $all_zip
  """
} 

/*
 Trim primer sequences.

 This process will look output read pairs that have one of the expected primer pairs at 
 the expected positions at the ends of R1 and R2.  These primers will be trimmed and 
 matching read pairs will be output to a new file.  This avoid amplification products that
 were formed by unexpected primer combinations (since this is a multiplex PCR assay).  

 This trimming uses cutadapt as described in:

 https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads

 This step also throws away read pairs shorter than ${params.post_trim_min_length} 

 This process will be run once per input fastq file pair per primer pair

*/

process trim_primer_seqs {                                                      
  label 'many_forks'
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
  }
                                                                              
  input:                                                                      
  val("indexes_complete") from post_index_setup_ch
  tuple val(primers), val(sample_id), path(initial_fastq) from primers_and_samples

  output:                                                                     
  tuple val(sample_id), path("*.R1_${primers.primer_name}.fastq.gz"), path("*.R2_${primers.primer_name}.fastq.gz") into primer_trimmed_ch_ungrouped
                                                                              
  script:                                                                     
  def primer_f_rc = primers.primer_f_seq.reverse().complement()                           
  def primer_r_rc = primers.primer_r_seq.reverse().complement()                           
  def f1 = initial_fastq[0]
  def f2 = initial_fastq[1]

  // setup arguments for trimming this particular pair of primers
  // see: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads
  def primer_args = "-a " + primers.primer_f_seq + "..." + primer_r_rc + " -A " + primers.primer_r_seq + "..."  + primer_f_rc
                                                                              
  """                                                                         
  # the --discard-untrimmed option means that the output of this command will contain
  # only those reads that match this particular expected primer pair because only
  # those matching read pairs will be trimmed, and non-trimmed pairs will get discarded.
  #
  # from the cutadapt manual:
  # "Use --discard-untrimmed to throw away all read pairs in which R1 
  # doesn’t start with FWDPRIMER or in which R2 does not start with REVPRIMER"
  #
  cutadapt $primer_args --discard-untrimmed  -e ${params.amplicon_primers_max_error_fraction} --minimum-length ${params.post_trim_min_length} $f1 $f2 -o ${sample_id}.R1_${primers.primer_name}.fastq.gz -p ${sample_id}.R2_${primers.primer_name}.fastq.gz 
  """                                                                         
}
                                                                                
/*
  This process concatenates the outputs from the individual cutadapt calls

  Cutadapt is run multiple times per fastq pair to create individual outputs 
  containing read pairs that contained an expected primer pair at their ends

  This step concatenates all those outputs into a single file per original
  fastq pair. 

  Because *all* cutadapt output flows into the primer_trimmed_ch_ungrouped channel,
  it is necessary to bin these outputs into groups corresponding to the original
  fastq read pairs.  

  The trick to doing this is the .groupTuple() operator, which
  groups tuples (sample_id + read1_trimmed_fastq + read2_trimmed_fastq) based on 
  a shared sample_id.

  This process also does one more round of trimming to remove stanard Illumina
  adapter sequences

*/
process collect_cutadapt_output {                                               
  publishDir "${params.trimmed_outdir}", mode:'link'                                    

  label 'process_low'
                                                                                
  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
  }

  input:
  // the groupTuple() operator here will consolidate tuples with a shared sample_id
  tuple val(sample_id), path(individual_r1), path(individual_r2) from primer_trimmed_ch_ungrouped.groupTuple()
                                                                                
  output:                                                                       
  val(sample_id) into post_trim_ch
  val(sample_id) into post_trim_sample_ch
  tuple val(sample_id), path("*_trimmed.fastq.gz") into post_trim_qc_ch
                                                                                
  script:                                                                       

  // cutadapt parameters to trim TruSeq-style adapters
  // see: https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage
  def truseq_cutadapt = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  
  // this shortened version of the Nextera adapter found in some amplicons
  def nextera_cutadapt = "-a CTGTCTCTTATACA -A CTGTCTCTTATACA"

  // these will be the names of the new concatenated outputs, before truseq trimming
  def f1 = "${sample_id}.R1_individual_primers.fastq.gz"
  def f2 = "${sample_id}.R2_individual_primers.fastq.gz"

  """                                                                           
  # concatenate the individually-trimmed files (one for each amplicon target)
  # note that cat will concatenate .gz files just fine
  cat $individual_r1 > $f1
  cat $individual_r2 > $f2

  # one last cutadapt command to remove adapter sequences and too-short read pairs
  cutadapt $nextera_cutadapt \
    -O ${params.adapters_min_overlap} \
    --minimum-length ${params.post_trim_min_length} \
    -o ${sample_id}_R1_trimmed.fastq.gz \
    -p ${sample_id}_R2_trimmed.fastq.gz \
    ${f1} ${f2}
  """                                                                           
}         


/*
 Use fastqc to do QC on post-trimmed fastq
*/
process post_trim_qc {
  label 'many_forks'
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"
  } else {
      container "quay.io/biocontainers/fastqc:0.11.9--0"
  }

  input:
  tuple val(sample_id), path(input_fastq) from post_trim_qc_ch

  output:
  path("*.zip") into post_trim_multiqc_ch

  script:

  """
  fastqc $input_fastq
  """
}

/*
 Use multiqc to merge post-trimming fastq reports

 This report should show, for instance, the absence of any remaining
 Illumina adapter sequences. 

 TODO: will this report old fastqc reports if the pipeline is re-run?

*/
process post_trim_multiqc {

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
  }

  input:
  path(all_zip) from post_trim_multiqc_ch.collect()

  output:
  path("post_trim_qc_report.html") into post_trim_multiqc_output_ch

  """
  multiqc --interactive -n "post_trim_qc_report.html" -m fastqc $all_zip
  """
}


/* 
  Run dada2 on trimmed read pairs to:

 - perform error-correction on misscalled bases
 - merge read pairs
 - cluster identical read pairs
 - remove chimeric sequences
 - tabulate the frequencies of particular unique sequences

*/


trimmed_fastq_ch = Channel.fromPath( params.trimmed_outdir )

process run_dada_on_trimmed {
  publishDir "${params.dada_outdir}", mode: 'link'

  label 'process_high'
  label 'error_retry'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0"
  } else {
      container "quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0"
  }

  input:
  val(all_sample_ids) from post_trim_ch.collect()
  path (trimmed_fastq_dir) from trimmed_fastq_ch

  output:
  path("dada_seqtab.txt") into post_dada_run_ch

  script:                                                                       
  """                                                                             
  # Run dada2 using trimmed fastq as input and create a tabular output of results
  Rscript ${params.script_dir}/run_dada_on_trimmed.R ${params.script_dir} $trimmed_fastq_dir
  """             
}


/* 
  Tidy up dada2 output and write out observed sequences (ASVs) in fasta format
*/

process tidy_dada_output {
  publishDir "${params.dada_outdir}", mode: 'link'

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  input:
  path(dada_output) from post_dada_run_ch

  output:
  tuple path("observed_sequences.fasta"), path("sequence_abundance_table.tsv") into post_dada_ch

  script:                                                                       
  """                                                                             
  # This R script creates an output named observed_sequences.fasta containing all of the
  # unique sequences observed in the amplicon dataset
  # and a sequence_abundance_table.tsv, which lists the abundances of these
  # sequences in each dataset
  Rscript ${params.script_dir}/tidy_dada_output.R ${params.script_dir} $dada_output
  """             
}

/* 
 This process aligns all of the unique sequences reported by dada2
 to the set of expected reference sequences using blastn.
*/
process compare_observed_sequences_to_ref_seqs {
  publishDir "${params.blast_outdir}", mode: 'link', pattern: "*bn_refseq"

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  input:
  tuple path(sequences), path(abundance_table) from post_dada_ch
  tuple path(refseq_blast_db), path (refseq_blast_db_files) from refseq_blast_db_ch

  output:
  tuple path(abundance_table), path("${sequences}.bn_refseq") into post_compare_ch

  script:                                                                       
  // this is similar to the default blastn output except gaps replaces gapopens, because seems more useful!
  // also include other useful columns like query length, subject length, etc
  def blastn_columns = "qaccver saccver pident length mismatch gaps qlen slen bitscore"
  """                                                                           
  blastn -db ${refseq_blast_db} -task blastn -evalue ${params.max_blast_refseq_evalue} -query $sequences -outfmt "6 $blastn_columns" -out ${sequences}.bn_refseq.no_header
  # prepend blast output with the column names so we don't have to manually name them later
  echo $blastn_columns > blast_header.no_perl
  echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header 
  cat blast_header ${sequences}.bn_refseq.no_header > ${sequences}.bn_refseq
  """             
}

surveillance_columns_ch =  Channel.fromPath(params.surveillance_columns, type: 'file', checkIfExists: true)
targets_ch =  Channel.fromPath(params.targets, type: 'file', checkIfExists: true)
targets_msa_ch =  Channel.fromPath(params.targets, type: 'file', checkIfExists: true)

/*
  This process takes the blast alignment information from the above process
  and decides whether those sequence are sufficiently similar to the reference 
  sequences to be assigned to them.
*/
process assign_observed_sequences_to_ref_seqs {

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  input:
  path(metadata) from validated_metadata_ch
  tuple path(abundance_table), path(blast_output) from post_compare_ch
  val(R_setup_OK) from post_r_dep_setup_ch
  path(surveillance_columns_file) from surveillance_columns_ch
  path(targets_file) from targets_ch

  output:
  path("unassigned_sequences.fasta") into post_assign_to_refseq_ch
  path("sequencing_report.xlsx") into report_output_ch
  path("all_data*.csv") into csv_output_ch
  path("*.txt") into assigned_txt_ch
  path("*.pdf") into assigned_pdf_ch

  // output channels for tree-building process
  path("sequencing_report.xlsx") into report_tree_ch

  script:                                                                       
  // only use R lib dir for singularity
  def r_lib_dir = workflow.containerEngine == 'singularity' ? "${params.R_lib_dir}" : "NA"
  """                                                                           
  Rscript ${params.script_dir}/assign_observed_seqs_to_ref_seqs.R ${params.script_dir} $r_lib_dir $abundance_table $blast_output $metadata ${targets_file} ${surveillance_columns_file} ${params.min_reads_for_positive_surveillance_call}
  """             
}

/*
   Split up assigned observed sequeces by target
   for making trees
*/
process create_fasta_for_trees {
  publishDir "${params.tree_outdir}", mode: 'link'

  // need to specify label?
  // label 'lowmem'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10.4" 
  } else {
      container "quay.io/biocontainers/python:3.10.4"
  }

  when:
  params.make_trees

  input:
  path(sequencing_report) from report_tree_ch
  val(venv_setup) from post_venv_setup_ch
  path(targets_file) from targets_msa_ch

  output:
  path("*_all.fasta") into fasta_tree_ch

  script:
  // only need to activate the venv for singularity
  def activate_venv_command = workflow.containerEngine == 'singularity' ? "source ${params.python_venv_path}/bin/activate" : ""
  """
  $activate_venv_command
  python3 ${params.script_dir}/MPAS_create_fasta.py $sequencing_report $targets_file
  """
}

/*
   Build multiple-sequencing alignments for each group of sequences using MAFFT. 

   MAFFT documentation : https://mafft.cbrc.jp/alignment/software/manual/manual.html
*/
process make_tree_alignment {
  publishDir "${params.tree_outdir}", mode: 'link'

  label 'process_medium'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/mafft:7.505--hec16e2b_0"
  } else {
      container "https://depot.galaxyproject.org/singularity/mafft:7.505--hec16e2b_0"
  }

  input:
  path(all_fasta) from fasta_tree_ch.flatten()
  
  output:
  path("mafft_${all_fasta}") into msa_tree_ch

  shell:
  """
  mafft --adjustdirection --quiet --auto --nuc "$all_fasta" > "mafft_${all_fasta}"
  """
}

/*
   Build maximum likelihood for each group of sequences using IQ-TREE. 

   IQ-TREE documentation: www.iqtree.org/doc/
*/
process make_ml_tree {

  label 'process_medium'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1"
  } else {
      container "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1"
  }
  
  input:
  path(all_alignment) from msa_tree_ch

  output:
  path("tree_${all_alignment.baseName}.treefile") into ml_tree_ch
  
  shell:
  """
  iqtree -s $all_alignment -st DNA -quiet -m MFP -pre tree_${all_alignment.baseName}   
  """ 
}

/*
   Creates pdf files of each ML tree using ToyTree. 

   ToyTree documentation: https://toytree.readthedocs.io/en/latest/
*/
process view_phylo_tree {
  publishDir "${params.tree_outdir}", mode: 'link'

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10.4" 
  } else {
      container "quay.io/biocontainers/python:3.10.4"
  }

  input:
  path(iqtree) from ml_tree_ch.flatten()

  output:
  path("${iqtree.baseName}.pdf") into pdf_tree_ch

  shell:
  // only need to activate the venv for singularity
  def activate_venv = workflow.containerEngine == 'singularity' ? "source ${params.python_venv_path}/bin/activate" : ""
  """
  # source ${params.python_venv_path}/bin/activate                                
  $activate_venv
  python3 ${params.script_dir}/MPAS_view_tree.py $iqtree
  """
}


/*
  Create a channel containing the path to an (optional) local copy of the nt
  blast database.
 */
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


/*
  This process confirms that the local NT database is valid 
  (if applicable: if going to blast unassigned sequences and 
   if not doing a remote blast)
 */

process check_local_blast_database {
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  when: 
  params.blast_unassigned_sequences 

  input:
  path(local_nt_database_dir) from blast_db_ch

  output:
  path(local_nt_database_dir) into post_blast_db_check_ch

  script:                                                                       
  def local_nt_database = "${local_nt_database_dir}/${params.local_nt_database_name}"

  params.remote_blast_nt ? 
  """
    echo "Don't need to check local BLAST database when running with -remote option"
    rm $local_nt_database_dir
    touch $local_nt_database_dir
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
  """
}

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


/* 
  Confirm that local blast will be taxonomically aware
 */
process check_blast_tax {
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/curl:7.80.0"
  } else {
      container "quay.io/biocontainers/curl:7.80.0"
  }

  when: 
  params.blast_unassigned_sequences

  input:
  tuple path(blast_tax_dir), val(existing_db) from blast_tax_ch

  output:
  path (blast_tax_dir) into post_blast_tax_check_ch

  script:
  // if a local blast_tax_dir is specified, check that it contains the expected files
  existing_db ? 
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
  
  """ :
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
    # get rid of archive
    rm taxdb.tar.gz
  """
}

/* 
 This process BLASTs all of the unassigned unique sequences 
 (amplicon sequence variants) from dada2 against the NCBI
 nt database to try to figure out what they are
*/

process blast_unassigned_sequences {
  publishDir "${params.blast_outdir}", mode: 'link'

  label 'process_high'
  label 'process_high_memory'
  label 'error_retry'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  when: 
  params.blast_unassigned_sequences

  input:
  path(unassigned_sequences) from post_assign_to_refseq_ch
  path(blast_tax_dir), stageAs: 'local_blast_tax_dir' from post_blast_tax_check_ch
  path(local_nt_database_dir), stageAs: 'local_nt_db_dir'  from post_blast_db_check_ch

  output:
  path("${unassigned_sequences}.bn_nt") into post_blast_unassigned_ch
  path("${unassigned_sequences}") into post_blast_unassigned_seq_ch

  script:                                                                       
  // if BLASTing locally: the location/name of the blast db
  def local_nt_database = "${local_nt_database_dir}/${params.local_nt_database_name}"

  // columns to include in blast output: standard ones plus taxonomy info
  def blastn_columns = "qaccver saccver pident length mismatch gaps qstart qend sstart send evalue bitscore staxid ssciname scomname sblastname sskingdom"
/*
   	    staxid means Subject Taxonomy ID
   	  ssciname means Subject Scientific Name
   	  scomname means Subject Common Name
   	sblastname means Subject Blast Name
   	 sskingdom means Subject Super Kingdom
*/

  def blast_db_params = ""
  if (params.remote_blast_nt) {
    // remote blastn: slower but doesn't require locally installed nt database
    blast_db_params = "-db nt -remote"
  }
  else {
    // local blastn: faster but requires locally installed nt database
    blast_db_params = "-db ${local_nt_database}"
  }


  """
  export BLASTDB="$blast_tax_dir"

  # run blastn
  blastn $blast_db_params -task megablast -perc_identity ${params.blast_perc_identity} -qcov_hsp_perc ${params.blast_qcov_hsp_perc} -evalue ${params.max_blast_nt_evalue} -query $unassigned_sequences -outfmt "6 $blastn_columns" -out ${unassigned_sequences}.bn_nt.no_header

  # prepend blast output with the column names so we don't have to manually name them later
  # the perl inline command here is to replace spaces with tabs
  echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header 

  # merge custom header line with actual blast output
  cat blast_header ${unassigned_sequences}.bn_nt.no_header > ${unassigned_sequences}.bn_nt

  """
}

/*
   This process parses the blast output from unassigned sequences and
   populates a spreadsheet with the info. 
*/

process assign_non_ref_seqs {

  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity'){
      container "docker://rocker/tidyverse:4.2.2"
  }     

  input:
  path(blast_output) from post_blast_unassigned_ch
  path(unassigned_sequences) from post_blast_unassigned_seq_ch

  output:
  path("non_reference_sequence_assignments.xlsx") into unassigned_blast_output_ch

  script:                                                                       
  // only use R lib dir for singularity
  def r_lib_dir = workflow.containerEngine == 'singularity' ? "${params.R_lib_dir}" : "NA"
  """                                                                           
  Rscript ${params.script_dir}/assign_non_ref_seqs.R $r_lib_dir $unassigned_sequences $blast_output
  """             
}

/* 
  This process prepends main output files with a prefix (by default the date, can be overridden)
*/

process prepend_output_filenames {
  publishDir "${params.outdir}", mode: 'link'

  label 'process_low'
   
  input:
  path(output_file) from initial_multiqc_output_ch.mix( post_trim_multiqc_output_ch,
                                                        report_output_ch,
                                                        csv_output_ch.flatten(),
                                                        assigned_txt_ch.flatten(),
                                                        assigned_pdf_ch.flatten(),
                                                        unassigned_blast_output_ch)
                                                        
  output:
  path ("${params.output_prefix}${output_file}")

  script:
  """
  cp $output_file ${params.output_prefix}${output_file}
  """
}

