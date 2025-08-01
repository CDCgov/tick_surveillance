/* 
 This workflow creates a fasta file containing all of the predefined reference
 sequences (defined in the targets file).
 */
workflow GENERATE_REFSEQ_FASTA {

  take: 
  targets

  main:
  GENERATE_REFSEQ_FASTAS(targets)
  COMBINE_REFSEQ_FASTAS(GENERATE_REFSEQ_FASTAS.out.fasta.collect())
  
  emit:
  fasta = COMBINE_REFSEQ_FASTAS.out.fasta
  // no specialized tool versions to emit
}

/*
  generate one fasta for each reference sequence
*/
process GENERATE_REFSEQ_FASTAS {
  label 'process_low'

  input:
  val targets 

  output:
  path ("*fasta"), emit: fasta

  // makes fasta formatted records for each target
  script:
  """
  printf ">%s\n%s\n" $targets.ref_sequence_name $targets.sequence > ${targets.ref_sequence_name}.fasta
  """
}

/*
  concatenate individual fasta into one big reference sequence file
*/
process COMBINE_REFSEQ_FASTAS {
  label 'process_low'
  publishDir "${params.outdir}", mode: "copy"

  input:
  path (individual_fastas) 

  output:
  path ("reference_sequences.fasta"), emit: fasta

  // makes fasta formatted records for each targets
  script:

  def num_fasta = individual_fastas.size()

  """
  cat $individual_fastas > reference_sequences.fasta
  """
}

