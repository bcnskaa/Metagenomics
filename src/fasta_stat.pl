#!/usr/bin/perl -w 
use strict;
use Bio::SeqIO;

# Read in command arguments
my $infn = $ARGV[0];

#print "Reading from $infn";

my $total_seq_len = 0;
my $max_seq_len = 0;
my $min_seq_len = 1000000000;
my $seq_n = 0;

my $in = Bio::SeqIO->new(-file=>$infn, -format=>'fasta');
while(my $seq = $in->next_seq() ){
  my $id = $seq->id();
  
  #print "$id\t",$seq->length,"\n";
  
  $total_seq_len += $seq->length;
  if($seq->length > $max_seq_len)
  {
  	$max_seq_len = $seq->length;
  }
  
  if($seq->length < $min_seq_len)
  {
  	$min_seq_len = $seq->length;
  }
  $seq_n += 1;
}

print "seq_n=", $seq_n, " mean_len=", $total_seq_len / $seq_n, " max_len=", $max_seq_len, " min_len=", $min_seq_len, "\n"; 
