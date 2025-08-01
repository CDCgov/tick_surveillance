/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VALIDATE INPUTS                                                             
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                

/*
   This uses the parameter schema defined in nextflow_schema.json
   to validate parameters.  It does things like make sure required
   parameters are defined and that parameters have appropriate values.
 */
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)            

// Validate input parameters                                                    
WorkflowMPAS.initialise(params, log)                                           

/*
  Read in the targets tsv-format file that describes the expected target sequences.
  and create channel
 */
Channel
    .fromPath(params.targets, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { targets_ch }

/*
  A channel containing the path of the tsv file specifying surveillance targets
  
  This channel differs from targets_ch in that targets_ch contains a groovy
  data structure from the parsed file.  This ch just contains the file's path.
 */
Channel
     .fromPath(params.targets, type: 'file', checkIfExists: true)
     .set{ targets_file_ch }

/*

 This channel contains the primer sequences that were
 used to amplify surveillance targets read from the primers.tsv file.  

 These primer sequences will be trimmed off of read pairs and used to 
 identify legitimate PCR products, which contain an expected F/R primer pair 
 at their ends.
 */
Channel
    .fromPath(params.primers, checkIfExists: true)
    .splitCsv(header:true, sep:"\t", strip:true)
    .set { primers_ch }

//  A channel containing the tsv file specifying surveillance report columns 
Channel
     .fromPath(params.surveillance_columns, type: 'file', checkIfExists: true)
     .set{surveillance_columns_ch}


/*
 These fastq files represent the main input to this workflow
 */

Channel
    .fromFilePairs("${params.fastq_dir}/${params.fastq_pattern}",
                   size: 2,
                   maxDepth: 1)
    .map { untrimmed_sample_id, fastq  ->
           // strip off any _L### text from the end of the sample ID
           // and then any _S## text from the end
           // this is text added onto IDs from samplesheet that are added by Illumina bcl2fastq
           def sample_id = untrimmed_sample_id.replaceFirst( /_L\d{3}$/, "")
           sample_id = sample_id.replaceFirst( /_S\d+$/, "")

           // define a new empty map named meta for each sample                   
           // and populate it with id and single_end values                       
           // for compatibility with nf-core module expected parameters           
           // reads are just the list of fastq                                    
           def meta        = [:]      

           // this map gets rid of any _S\d+ at the end of sample IDs but leaves fastq
           // names alone.  E.g. strip _S1 from the end of a sample ID..          
           // This is typically sample #s from Illumina basecalling.              
           // could cause an issue if sequenced the same sample with              
           // multiple barcodes so was repeated on a sample sheet.                
           meta.id         = sample_id
                                                                                  
           // this pipeline only designed to work with paired-end data
           meta.single_end =  false

           // this last statement in the map closure is the return value          
           [ meta, fastq ] }     
      
    .set {reads_ch}

/*
  collect sample IDs in a file
*/

reads_ch.map{meta, reads -> meta.id}
  .collectFile(name: 'sample_ids.txt', newLine: true)
  .set{ sample_ids_file_ch }

/*
   this combinatorially mixes the fastq file pairs with the primer sequences
   to be trimmed, creating a new channel with all possible combinations of
   input fastq and primer pair
*/                 
    
// primers_and_reads_ch = primers_ch.combine(reads_ch)

// set channel for filter blast report process
Channel
    .value( params.filter_unassigned_seq )
    .collect()
    .set { filter_ch }

/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES                                                                
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
                                                                                
/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS                                           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              

                                                                                
// include modules and subworkflows

// setup indexes and R/python environments
include { SETUP_INDEXES                 } from '../subworkflows/setup_indexes'
include { SETUP_PYTHON_ENVIRONMENT      } from '../modules/local/setup_python_env/main'
include { SETUP_R_DEPENDENCIES          } from '../modules/local/setup_R_dependencies/main'

// make a fasta of all reference sequences
include { GENERATE_REFSEQ_FASTA         } from '../subworkflows/generate_refseq_fasta'

// validate metadata
include { VALIDATE_METADATA             } from '../modules/local/validate_metadata/main'

// primer, quality, and adapter trimming
include { GENERATE_PRIMER_FILES         } from '../subworkflows/generate_primers_file'
include { TRIM_READS                    } from '../subworkflows/trim_reads'

// dada2
include { DADA2                         } from '../modules/local/dada2/main'
include { TIDY_DADA_OUTPUT              } from '../modules/local/dada2/main'

// compare ASVs to expected sequences and make positive/negative calls
include { COMPARE_OBSERVED_SEQS         } from '../modules/local/assign_sequences/main'
include { ASSIGN_OBSERVED_SEQS          } from '../modules/local/assign_sequences/main'

// tree-building
include { GENERATE_TREES                } from '../subworkflows/generate_trees'

// BLAST unassigned sequences
include { SETUP_BLAST_DB_AND_TAX        } from '../subworkflows/setup_blast_db_and_tax'
include { CLASSIFY_UNASSIGNED_SEQUENCES } from '../subworkflows/classify_unassigned_sequences'

include { PREPEND_OUTPUT_FILENAMES      } from '../modules/local/prepend_filenames/main'

// include { ORG_UNASSIGNED_SEQUENCES      } from '../modules/local/org_unassigned_sequences/main'
// include { FILTER_UNASSIGNED_SEQUENCES   } from '../modules/local/filter_unassigned_sequences/main'

// nf-core modules 
include { FASTQC as FASTQC_PRE        } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST       } from '../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_PRE      } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_POST     } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_TRIM     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*                                                                              
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW                                                           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/                                                                              
                                                                                
// Info required for completion email and summary                               
def multiqc_report1 = []
def multiqc_report2 = []                                                        
def multiqc_report3 = []                                                        
                                                                                
workflow MPAS_WORKFLOW {                                                                
                                                                                
    // this will keep track of all software versions
    ch_versions = Channel.empty()                                               
                                                                                
    CUSTOM_DUMPSOFTWAREVERSIONS (                                               
        ch_versions.unique().collectFile(name: 'collated_versions.yml')         
    )                     

    // generate fasta file for all reference sequences defined in targets file
    GENERATE_REFSEQ_FASTA(targets_ch)

    // create a blast db of reference sequences
    SETUP_INDEXES (GENERATE_REFSEQ_FASTA.out.fasta)
    ch_versions = ch_versions.mix ( SETUP_INDEXES.out.versions )      

    // setup python venv to handle python dependencies for tree-building
    python_req_ch = Channel.fromPath(params.python_requirements, checkIfExists: true)
    SETUP_PYTHON_ENVIRONMENT(params.python_venv_path, python_req_ch)
    ch_versions = ch_versions.mix ( SETUP_PYTHON_ENVIRONMENT.out.versions )      

    // setup R packages not included in base environment 
    R_tar_dir_ch = Channel.fromPath(params.R_tar_dir, checkIfExists: true)
    SETUP_R_DEPENDENCIES(R_tar_dir_ch)
    ch_versions = ch_versions.mix ( SETUP_R_DEPENDENCIES.out.versions )      

    // Check for existence of metadata file and validate it
    metadata_ch = Channel.fromPath("${params.metadata}", checkIfExists: true)
    VALIDATE_METADATA(metadata_ch, sample_ids_file_ch)
    ch_versions = ch_versions.mix ( VALIDATE_METADATA.out.versions )      

    // run fastqc on input reads
    FASTQC_PRE ( reads_ch.map {meta, reads -> [meta, reads, "pre_trimming"] } )
    ch_versions = ch_versions.mix ( FASTQC_PRE.out.versions )      

    // Create files with primers to be trimmed.
    // These files will be input to cutadapt for trimming
    GENERATE_PRIMER_FILES(primers_ch)

    // trim reads of primers, adapters, and low quality bases 
    TRIM_READS ( reads_ch, GENERATE_PRIMER_FILES.out.f_primer_file, GENERATE_PRIMER_FILES.out.r_primer_file )
    ch_versions = ch_versions.mix ( TRIM_READS.out.versions )      
 
    // run fastqc on trimmed reads
    FASTQC_POST ( TRIM_READS.out.trimmed_reads.map {meta, reads -> [meta, reads, "post_trimming"] } )
    ch_versions = ch_versions.mix ( FASTQC_POST.out.versions )      

    // run dada2 on trimmed reads
    ch_all_trimmed_reads = TRIM_READS.out.trimmed_reads
                            .map { meta, reads -> [ reads ] }                                   
                            .flatten().collect()

    // params related to dada2 error correction, filtering, and filtering
    def maxN         = params.dada2_maxN
    def maxEE        = params.dada2_maxEE
    def truncQ       = params.dada2_truncQ
    def trimRight    = params.dada2_trimRight
    def min_reads    = params.dada2_min_reads
    def min_overlap  = params.dada2_min_overlap
    def max_mismatch = params.dada2_max_mismatch

    // run DADA2 on trimmed reads
    DADA2(ch_all_trimmed_reads, maxN, maxEE, truncQ, trimRight, min_reads, min_overlap, max_mismatch)
    ch_versions = ch_versions.mix ( DADA2.out.versions )      

    // tidy dada2 output 
    TIDY_DADA_OUTPUT(DADA2.out.seqtab, DADA2.out.dada_read_tracking_all)
    ch_versions = ch_versions.mix ( TIDY_DADA_OUTPUT.out.versions )      

    // compare (align) observed sequences to predefined reference sequences
    COMPARE_OBSERVED_SEQS(TIDY_DADA_OUTPUT.out.observed_sequences,
                          TIDY_DADA_OUTPUT.out.sequence_abundance_table,
                          GENERATE_REFSEQ_FASTA.out.fasta,
                          SETUP_INDEXES.out.db)
    ch_versions = ch_versions.mix ( COMPARE_OBSERVED_SEQS.out.versions )      

    // logic to assign observed sequences to reference sequences
    // and generate surveillance report
    ASSIGN_OBSERVED_SEQS (VALIDATE_METADATA.out.validated_metadata,
                          COMPARE_OBSERVED_SEQS.out.sequence_abundance_table,
                          COMPARE_OBSERVED_SEQS.out.blast_output,
                          SETUP_R_DEPENDENCIES.out.R_lib_dir,
                          surveillance_columns_ch,
                          targets_file_ch)
    ch_versions = ch_versions.mix ( ASSIGN_OBSERVED_SEQS.out.versions )      

    // generate phylogenetic trees
    GENERATE_TREES (ASSIGN_OBSERVED_SEQS.out.surveillance_report,
                    SETUP_PYTHON_ENVIRONMENT.out.venv_path,
                    targets_file_ch)      
    ch_versions = ch_versions.mix ( GENERATE_TREES.out.versions )      

    // setup BLAST db for classification of unassigned sequences
    SETUP_BLAST_DB_AND_TAX()
    ch_versions = ch_versions.mix ( SETUP_BLAST_DB_AND_TAX.out.versions )      

    // BLAST (vs. NCBI nt database) then classify unassigned sequences
    CLASSIFY_UNASSIGNED_SEQUENCES(ASSIGN_OBSERVED_SEQS.out.unassigned_sequences,
                                  SETUP_BLAST_DB_AND_TAX.out.blast_db_dir,
                                  SETUP_BLAST_DB_AND_TAX.out.blast_tax_dir,
                                  SETUP_R_DEPENDENCIES.out.R_lib_dir,
                                  COMPARE_OBSERVED_SEQS.out.sequence_abundance_table,
                                  VALIDATE_METADATA.out.validated_metadata,
                                  ASSIGN_OBSERVED_SEQS.out.surveillance_report,
                                  filter_ch)

    ch_versions = ch_versions.mix ( CLASSIFY_UNASSIGNED_SEQUENCES.out.versions )      

    /*
    // Update unassigned_sequences output file with metadata info
    ORG_UNASSIGNED_SEQUENCES(COMPARE_OBSERVED_SEQS.out.sequence_abundance_table,
                             VALIDATE_METADATA.out.validated_metadata,
                             CLASSIFY_UNASSIGNED_SEQUENCES.out.report,
                             SETUP_R_DEPENDENCIES.out.R_lib_dir)

    // Filter unassigned_sequence_report for target name. Run if parameter by user given
    FILTER_UNASSIGNED_SEQUENCES(ASSIGN_OBSERVED_SEQS.out.surveillance_report,
                                ORG_UNASSIGNED_SEQUENCES.out.org_unassigned_sequences_report,
                                filter_ch,
                                SETUP_R_DEPENDENCIES.out.R_lib_dir)
    */

    //                                                                          
    // MODULE: MultiQC                                                          
    //                                                                          
    workflow_summary    = WorkflowMPAS.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)                       
                                                                                
    methods_description    = WorkflowMPAS.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)                 
                                                                                
    ch_multiqc_files1 = Channel.empty()                                          
    ch_multiqc_files1 = ch_multiqc_files1.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files1 = ch_multiqc_files1.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files1 = ch_multiqc_files1.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files1 = ch_multiqc_files1.mix(FASTQC_PRE.out.zip.collect{it[1]}.ifEmpty([]))
    
    ch_multiqc_files2 = Channel.empty()                                          
    ch_multiqc_files2 = ch_multiqc_files2.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files2 = ch_multiqc_files2.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files2 = ch_multiqc_files2.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files2 = ch_multiqc_files2.mix(FASTQC_POST.out.zip.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files3 = Channel.empty()                                          
    ch_multiqc_files3 = ch_multiqc_files3.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files3 = ch_multiqc_files3.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files3 = ch_multiqc_files3.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files3= ch_multiqc_files3.mix(TRIM_READS.out.cutadapt_trim_report.collectFile(name: '*_summary.txt'))
                                                                                
    MULTIQC_PRE (                                                                   
        ch_multiqc_files1.collect(),                                             
        ch_multiqc_config.toList(),                                             
        ch_multiqc_custom_config.toList(),                                      
        ch_multiqc_logo.toList()                                                
    )
    MULTIQC_POST (                                                                   
        ch_multiqc_files2.collect(),                                             
        ch_multiqc_config.toList(),                                             
        ch_multiqc_custom_config.toList(),                                      
        ch_multiqc_logo.toList()                                                
    )
    MULTIQC_TRIM (                                                                   
        ch_multiqc_files3.collect(),                                             
        ch_multiqc_config.toList(),                                             
        ch_multiqc_custom_config.toList(),                                      
        ch_multiqc_logo.toList()                                                
    )

    multiqc_report1 = MULTIQC_PRE.out.report.collectFile(name: 'MultiQC_FastQC_PRE_TRIM.html', storeDir: "${params.multiQC_reports}")
    multiqc_report2 = MULTIQC_POST.out.report.collectFile(name: 'MultiQC_FastQC_POST_TRIM.html', storeDir: "${params.multiQC_reports}")
    multiqc_report3 = MULTIQC_TRIM.out.report.collectFile(name: 'MultiQC_Cutadapt_Report.html', storeDir: "${params.multiQC_reports}")                                    
                                                              

    // prepend main output file names
    ch_main_output_files = Channel.empty()

    ch_main_output_files = ch_main_output_files
                             .mix( 
                             ASSIGN_OBSERVED_SEQS.out.surveillance_report,
                             ASSIGN_OBSERVED_SEQS.out.all_data_csv,
                             ASSIGN_OBSERVED_SEQS.out.txt,
                             ASSIGN_OBSERVED_SEQS.out.pdf,
                             CLASSIFY_UNASSIGNED_SEQUENCES.out.sequences_report_filter,
                             CLASSIFY_UNASSIGNED_SEQUENCES.out.org_unassigned_sequences_report)
                             .flatten()
    PREPEND_OUTPUT_FILENAMES(ch_main_output_files)

}                                                                               

/*
 Upon workflow completion 
*/
workflow.onComplete {                                                           
  // print summary message 
  NfcoreTemplate.summary(workflow, params, log)                               
}             

