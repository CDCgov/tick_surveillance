#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $print_usage = 0;
my $error_file = undef;
my $quality_score = 35;
my $num_reads = 10;
my $read_length = 250;
my $paired_reads = 1;
my $fasta_in = undef;
my $fastq_prefix = undef;

my $usage = <<USAGE;

  This script will perform simple fastq simulation from a single fasta reference sequence.   

  It will create a certain # of fastq reads or read pairs that start at the end(s) of the 
  reference sequence.  It is possible to specify length, # of reads, and a simple error
  model.  

  The point of this tool is not to make as realistic as possible fastq data but to make
  fastq data that is sufficiently realistic for testing purposes and to serve as an 
  internal control. 

  Mark Stenglein,  10/21/2021

  Usage: $0 [-h] 

   [-e]          File describing expected error rates, in the format output by dada2, 
                 with lines for each possible substitution (or non-substitution) and
		 their expected frequencies, like this:

		 A2A     0.989457860532376
                 A2C     0.00194919206215395

		 If not provided, no errors will be introduced
		 default: no errors

   [-f]          File with a single reference sequence in FASTA format

   [-q]          Average quality score for simulated fatq
                 default: $quality_score

   [-n]          Number of reads to simulate
                 default: $num_reads

   [-l]          Maximum read length to simulate.  
                 The length of reads will be the maximum of this and the refseq length.
                 default: $read_length

   [-p]          Paired read simulation
                 default: $paired_reads

   [-pre]        Prefix for fastq output files, 
                 which will have _R1.fastq and _R2.fastq (for paired reads) suffixes
                 default: name of fasta input file [-f]

   [-h]          print this message

USAGE

if ((scalar @ARGV == 0) and -t STDIN) { print $usage and exit; }

GetOptions ("h" => \$print_usage,
            "e=s" => \$error_file,
            "q=i" => \$quality_score,
            "n=i" => \$num_reads,
            "l=i" => \$read_length,
            "f=s" => \$fasta_in,
            "p" => \$paired_reads,
            "pre=s" => \$fastq_prefix
    );



if (!$fasta_in)
{
   die ("error: must define a fasta input file\n$usage\n"); 
}

open (my $fasta_fh, "<", $fasta_in) or die ("Error: couldn't open error file: $error_file\n$usage\n");

my $refseq = "";
my $num_refseq = 0;
while (<$fasta_fh>)
{
   chomp;
   if (/^>/) 
   { 
      $num_refseq += 1;
      if ($num_refseq > 1)
      {
         die ("error: only a single reference sequence supported.  >1 sequence in $fasta_in\n");
      }
      next;
   } 
   $refseq .= $_;
}

my %error_cdf = ();

# -------------------------------------------
# read in error file and store error profile
# -------------------------------------------
if ($error_file)
{
   open(my $err_fh, "<", "$error_file") or die ("Error: couldn't open error file: $error_file\n$usage\n");

   my %subs = ();
   while (<$err_fh>)
   {
      chomp;
      my @fields = split;
      if (scalar @fields != 2) { die ("Error: unexepected format for error file on line: $_\n"); }
      my $sub = $fields[0];
      my $rate = $fields[1];
      if ($sub =~ /([ACGT])2([ACGT])/)
      {

         $subs{$1}{$2} = $rate;
      }
      else { die ("Error: invalid substitution format for error file on line: $_\n"); }
   }

   # record error CDF
   my @bases = qw (A C G T);
   foreach my $base (@bases)
   {
      my @outcome_bases = sort { $subs{$base}{$b} <=> $subs{$base}{$a} } keys %{$subs{$base}}; 
      my $cumulative_probability = 0;
      foreach my $outcome_base (@outcome_bases)
      {
	 $cumulative_probability += $subs{$base}{$outcome_base};
	 $error_cdf{$base}{$outcome_base} = $cumulative_probability;
	 # print "$base\t$outcome_base\t$subs{$base}{$outcome_base}\t$cumulative_probability\n";
      }
   }
}

if (!$fastq_prefix)
{
   $fastq_prefix = $fasta_in;
}

my $out_1_name = $fastq_prefix."_R1.fastq";
my $out_2_name = $fastq_prefix."_R2.fastq";

# open output files for writing out fastq
open (my $out_1, ">", $out_1_name) or die ("Error: couldn't open output file $out_1_name\n");
my $out_2 = undef;
if ($paired_reads)
{
   open ($out_2, ">", $out_2_name) or die ("Error: couldn't open output file $out_1_name\n");
}

# RC a DNA sequence
sub reverse_complement
{
   my $original_sequence = shift @_;
   my $rc = reverse $original_sequence;
   $rc =~ tr/ATGCatgc/TACGtacg/;
   return ($rc);
}

my $rc_refseq = reverse_complement($refseq);

# ----------------------
# create simulated reads
# ----------------------
for (my $i = 0; $i < $num_reads; $i++)
{
   # don't make reads longer than the refseq
   if ($read_length > length($refseq))
   {
      $read_length = length($refseq);
   }
   
   # 
   # reads start at ends of ref seq to simulate amplicon data
   # 
   # this strategy could be modified to simulate shotgun data
   # 
   my $fwd_seq_no_error = substr($refseq, 0, ($read_length - 1));
   warn "fwd_seq_no_error: $fwd_seq_no_error\n";
   my $fwd_seq = apply_errors($fwd_seq_no_error);
   warn "fwd_seq: $fwd_seq\n";
   output_fastq_record($out_1, ($i+1), $fwd_seq);

   # make paired read and output if doing that kind of thing
   my $rev_seq_no_error = substr($rc_refseq, 0, ($read_length - 1));
   warn "rev_seq_no_error: $rev_seq_no_error\n";
   my $rev_seq = apply_errors($rev_seq_no_error);
   warn "rev_seq: $rev_seq\n";
   if ($paired_reads)
   {
      output_fastq_record($out_2, ($i+1), $rev_seq);
   }
}

# -----------------------------
# output a single fastq record
# -----------------------------
sub output_fastq_record
{
   my $fh = shift @_;
   my $seq_name = shift @_;
   my $seq = shift @_;

   print $fh "\@$seq_name\n";
   print $fh "$seq\n";
   print $fh "+\n";

   # Quality scores
   # Phred+33
   # COULD_DO: make a more realistic Q-score distribution
   my $quality_char = chr($quality_score + 33);
   for (my $i = 0; $i < length ($seq); $i++)
   {
      print $fh "$quality_char";
   }
   print $fh "\n";
}

close($out_1);
close($out_2);

# --------------------------------------
# apply a simple error model if defined
# --------------------------------------
sub apply_errors
{
   my $seq = shift @_;
   if (!defined $error_file) { warn "seq: $seq\n"; return $seq; }

   my $new_seq = "";
   for (my $c = 0; $c < length($seq); $c++)
   {
      my $base = substr($seq, $c, 1); 
      my $p = rand;
      my @possible_bases = sort { $error_cdf{$base}{$a} <=> $error_cdf{$base}{$b} } keys %{$error_cdf{$base}}; 
      foreach my $possible_base (@possible_bases)
      {
         my $base_p = $error_cdf{$base}{$possible_base};
         if ($p < $base_p)
	 {
	    warn "base: $base possible_base: $possible_base p: $p base_p: $base_p \n";
	    $new_seq .= "$possible_base";
	    last;
	 }
      }
   }
   return $new_seq;
}

