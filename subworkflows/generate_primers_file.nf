/* 
 This workflow creates a file containing all of the predefined primer
 sequences (defined in the primers file).

 This file will be input to cutadapt for primer trimming
 */
workflow GENERATE_PRIMER_FILES {

  take: 
  primers

  main:
  // this process takes one line of primers.tsv and outputs a file
  // in the format expected by cutadapt 
  GENERATE_ONE_PRIMER_FILE(primers)
 
  // This process combines all the individual primer files into one
  // forward-orientation file and one reverse-orientation file
  // The input is a sorted list so that the fwd and rev primers end up
  // in the same order in the 2 files.
  COMBINE_PRIMER_FILES(GENERATE_ONE_PRIMER_FILE.out.f_primer_file.toSortedList(), 
                       GENERATE_ONE_PRIMER_FILE.out.r_primer_file.toSortedList())
  
  emit:
  f_primer_file = COMBINE_PRIMER_FILES.out.f_primer_file
  r_primer_file = COMBINE_PRIMER_FILES.out.r_primer_file
}



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
  generate one primer file
*/
process GENERATE_ONE_PRIMER_FILE {
  label 'process_low'
  tag   "${primer.primer_name}"                                                       


  input:
  val primer  // a row from the primers.csv file, 
              // with names according to the header (first row) in that file

  output:
  path ("*.fwd.txt"), emit: f_primer_file
  path ("*.rev.txt"), emit: r_primer_file

  // makes 
  script:

  def primer_name = primer.primer_name
  def primer_f    = primer.primer_f_seq
  def primer_r    = primer.primer_r_seq
  def primer_f_rc = primer_f.reverse().complement()
  def primer_r_rc = primer_r.reverse().complement()
                                 
  """
  printf ">%s_fwd\n^%s...%s\n" $primer_name $primer_f $primer_r_rc > ${primer_name}.fwd.txt
  printf ">%s_rev\n^%s...%s\n" $primer_name $primer_r $primer_f_rc > ${primer_name}.rev.txt
  """
}

/*
  concatenate individual files into one file
*/
process COMBINE_PRIMER_FILES {
  label 'process_low'
  publishDir "${params.outdir}", mode: "copy"

  input:
  path (individual_f_files) 
  path (individual_r_files) 

  output:
  path ("fwd_orientation_primer_sequences.txt"), emit: f_primer_file
  path ("rev_orientation_primer_sequences.txt"), emit: r_primer_file

  // makes fasta formatted records for each targets
  script:

  """
  cat $individual_f_files > fwd_orientation_primer_sequences.txt
  cat $individual_r_files > rev_orientation_primer_sequences.txt
  """
}

