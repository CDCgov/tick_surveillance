/* 
 This workflow performs all quality and primer trimming
 */
workflow TRIM_READS {

  take: 
  reads         // [meta, [fastq_path]]
  f_primer_file // [path]
  r_primer_file // [path]

  main:
  TRIM_PRIMERS_AND_LOW_Q(reads, f_primer_file, r_primer_file)
  TRIM_ADAPTERS(TRIM_PRIMERS_AND_LOW_Q.out.trimmed_reads)
  
  emit:
  trimmed_reads        = TRIM_ADAPTERS.out.trimmed_reads
  cutadapt_trim_report = TRIM_PRIMERS_AND_LOW_Q.out.cutadapt_trim_report
  versions             = TRIM_PRIMERS_AND_LOW_Q.out.versions
}


/*
 Trim primer sequences and low quality.

 This process will look output read pairs that have one of the expected primer pairs at
 the expected positions at the ends of R1 and R2.  These primers will be trimmed and
 matching read pairs will be output to a new file.  This avoid amplification products that
 were formed by unexpected primer combinations (since this is a multiplex PCR assay).

 This trimming uses cutadapt as described in:

 https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads

 This step also throws away read pairs shorter than ${params.post_trim_min_length}

 This process will be run once per input fastq file pair 
*/

process TRIM_PRIMERS_AND_LOW_Q {
  label 'process_low'
  tag   "${meta.id}"
  publishDir "${params.cutadapt_trim_reports}", pattern: '*_summary.txt', mode: 'copy'

  // if using conda
  conda "${moduleDir}/cutadapt_environment.yml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:4.9--py312hf67a6ed_2"
  } else {
    container "quay.io/biocontainers/cutadapt:4.9--py312hf67a6ed_2"
  }
                       
  input:
  tuple val(meta), path(initial_fastq) 
  path(f_primer_file) 
  path(r_primer_file) 
 
  output:                                                                     
  // output the R1 and R2s as seperate list elements rather than a combined sub-list
  // because we want to handle them separately later
  tuple val(meta), 
        path("${meta.id}.primer_trimmed.R1.fastq.gz"), 
        path("${meta.id}.primer_trimmed.R2.fastq.gz")      , emit: trimmed_reads
  path "versions.yml"                                      , emit: versions
  path("${meta.id}.primer_trimming_summary.txt")           , emit: cutadapt_trim_report                                         

                                             
  script:
  def f1 = initial_fastq[0]
  def f2 = initial_fastq[1]

  // setup arguments for trimming this particular pair of primers
  // see: https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-both-ends-of-paired-end-reads
  // def primer_args = "-a " + primers.primer_f_seq + "..." + primer_r_rc + " -A " + primers.primer_r_seq + "..."  + primer_f_rc
  def primer_args = "-a file:${f_primer_file} -A file:${r_primer_file}"


  """
  # the --discard-untrimmed option means that the output of this command will contain
  # only those reads that match this particular expected primer pair because only
  # those matching read pairs will be trimmed, and non-trimmed pairs will get discarded.
  #
  # from the cutadapt manual:
  # "Use --discard-untrimmed to throw away all read pairs in which R1 
  # doesn’t start with FWDPRIMER or in which R2 does not start with REVPRIMER"
  #
  # --pair-adapters forces the primers applied to read1 and read2 to be 
  # considered together.  Note that this means that the primers must 
  # be supplied in the same order in the forward and reverse files
  # (done with toSortedList in GENERATE_PRIMERS_FILE).
  #
  # --trim-n removes N bases from the ends of reads
  #
  # -q and -Q trim low quality bases from the ends of read1 and read2
  #
  cutadapt \
    $primer_args \
    --discard-untrimmed  \
    --pair-adapters \
    --trim-n \
    -q ${params.basecall_quality_limit} \
    -Q ${params.basecall_quality_limit} \
    -e ${params.amplicon_primers_max_error_fraction} \
    --minimum-length ${params.post_trim_min_length} \
    $f1 \
    $f2 \
    -o ${meta.id}.primer_trimmed.R1.fastq.gz \
    -p ${meta.id}.primer_trimmed.R2.fastq.gz \
    > ${meta.id}.primer_trimming_summary.txt

  cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
    cutadapt: \$(cutadapt --version)
  END_VERSIONS         
  """
}

/*
  This process does one more round of trimming to remove standard Illumina
  adapter sequences, should they remain
*/
process TRIM_ADAPTERS {      
  label      'process_low'
  tag        "${meta.id}"
  publishDir "${params.trimmed_outdir}", pattern: "*.fastq.gz", mode: "link"
                                                                   
  // if using conda
  conda "${moduleDir}/cutadapt_environment.yml"                                          

  // if using singularity
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:4.9--py312hf67a6ed_2"
  } else {                                                                      
    container "quay.io/biocontainers/cutadapt:4.9--py312hf67a6ed_2"             
  } 

  input:
  tuple val(meta), path(R1_fastq), path(R2_fastq)
                                                                                
  output:                                                                       
  tuple val(meta), path("*.fully_trimmed.R[12].fastq.gz"), emit: trimmed_reads 
  path "versions.yml"                                    , emit: versions                                         
                                                          
  script:                                                                       

  // cutadapt parameters to trim TruSeq-style adapters
  // see: https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage
  def truseq_cutadapt = "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  
  // this shortened version of the Nextera adapter found in some amplicons
  def nextera_cutadapt = "-a CTGTCTCTTATACA -A CTGTCTCTTATACA"

  """
  # one last cutadapt command to remove adapter sequences and too-short read pairs
  cutadapt \
    $nextera_cutadapt \
    -O ${params.adapters_min_overlap} \
    --minimum-length ${params.post_trim_min_length} \
    -o ${meta.id}.fully_trimmed.R1.fastq.gz \
    -p ${meta.id}.fully_trimmed.R2.fastq.gz \
    ${R1_fastq} ${R2_fastq}

  cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
    cutadapt: \$(cutadapt --version)
  END_VERSIONS         
  """
}

