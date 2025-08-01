include { BLAST_MAKEBLASTDB } from '../modules/nf-core/blast/makeblastdb/main'                                                                                                                                                       
/*
   Setup indexes and dictionaries needed by downstream processes.
*/
workflow SETUP_INDEXES {

  take: 
  refseq_fasta

  main:
  BLAST_MAKEBLASTDB(refseq_fasta)
  
  emit:
  db       = BLAST_MAKEBLASTDB.out.db
  versions = BLAST_MAKEBLASTDB.out.versions
}
