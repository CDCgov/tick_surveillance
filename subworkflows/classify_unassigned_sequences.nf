include { BLASTN_UNASSIGNED_SEQUENCES   } from '../modules/local/blast_unassigned_sequences/main'
include { ASSIGN_UNASSIGNED_SEQUENCES   } from '../modules/local/assign_unassigned_sequences/main'
include { ORG_UNASSIGNED_SEQUENCES      } from '../modules/local/org_unassigned_sequences/main'
include { FILTER_UNASSIGNED_SEQUENCES   } from '../modules/local/filter_unassigned_sequences/main'


/*
 This workflow taxonomically classifies unassigned sequences (those sequences
 that were observed but not assigned to any of the predefined reference seqs).
 */
workflow CLASSIFY_UNASSIGNED_SEQUENCES {

  take: 
  unassigned_sequences
  blast_db_dir
  blast_tax_dir
  R_lib_dir
  sequence_abundance_table
  metadata
  surveillance_report
  taxa_to_filter

  main:

  ch_versions = Channel.empty()                           

  // filter BLAST results by taxids of interest?
  // initialize param value to empty string
  ch_taxids_of_interest = Channel.value("")
  if (params.taxids_of_interest) {
     // if taxids of interest defined, create new value channel with string
     ch_taxids_of_interest = Channel.value(params.taxids_of_interest)
  }  

  BLASTN_UNASSIGNED_SEQUENCES(unassigned_sequences, blast_db_dir, blast_tax_dir, ch_taxids_of_interest) 
  ch_versions = ch_versions.mix(BLASTN_UNASSIGNED_SEQUENCES.out.versions)

  ASSIGN_UNASSIGNED_SEQUENCES(BLASTN_UNASSIGNED_SEQUENCES.out.blast_out, unassigned_sequences, R_lib_dir)
  ch_versions = ch_versions.mix(ASSIGN_UNASSIGNED_SEQUENCES.out.versions)

  // Update unassigned_sequences output file with metadata info               
  ORG_UNASSIGNED_SEQUENCES(sequence_abundance_table,
                           metadata,
                           ASSIGN_UNASSIGNED_SEQUENCES.out.unassigned_sequences_report,
                           R_lib_dir)             
                                                                                
  // Filter unassigned_sequence_report for target name. Run if parameter by user given
  FILTER_UNASSIGNED_SEQUENCES(surveillance_report,
                              ORG_UNASSIGNED_SEQUENCES.out.org_unassigned_sequences_report,
                              taxa_to_filter,                                      
                              R_lib_dir)      

  emit:
  report                          = ASSIGN_UNASSIGNED_SEQUENCES.out.unassigned_sequences_report
  org_unassigned_sequences_report = ORG_UNASSIGNED_SEQUENCES.out.org_unassigned_sequences_report
  sequences_report_filter         = FILTER_UNASSIGNED_SEQUENCES.out.sequences_report_filter

  versions   = ch_versions
}


