#!/usr/bin/perl -w 
use strict;
use Bio::SeqIO;


if(@ARGV != 1)
{
	print "Missing argument\n";
	print "Usage: " . $ARGV[0] . " FASTA-FILES\n";
	exit;
}

# Read in command arguments
my $infn = $ARGV[0];

unless (-e $infn)
{
	print $infn . " does not exist.\n";
	exit;
}



#print "Reading from $infn";

my $total_seq_len = 0;
my $max_seq_len = 0;
my $min_seq_len = 1000000000;
my $seq_n = 0;

my $in = Bio::SeqIO->new(-file=>$infn, -format=>'fasta');
my @seq_lens; # For calculating N50 and N90 values
while(my $seq = $in->next_seq() ){
  my $id = $seq->id();
  
  #print "$id\t",$seq->length,"\n";
  
  $total_seq_len += $seq->length;
  
  push @seq_lens, $seq->length;
  
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

# Sort the sequence length
@seq_lens = sort{$b<=>$a} @seq_lens;

# Calculate N50 and N90 values
my ($count, $half, $n50, $n90)=(0,0,0,0);
for (my $j = 0;$j < @seq_lens;$j++){
    $count+=$seq_lens[$j];
    
    if (($count >= $total_seq_len / 2)&&($half==0)) {
        #print "N50: $seq_lens[$j]\n";
        $n50 = $seq_lens[$j];
        $half = $seq_lens[$j];
    } elsif ($count >= $total_seq_len*0.9) {
        $n90 = $seq_lens[$j];
    }
}



#print "seq_n=", $seq_n, " mean_len=", $total_seq_len / $seq_n, " max_len=", $max_seq_len, " min_len=", $min_seq_len, $n50, $n90"\n"; 
print "$infn\tseq_n=", $seq_n, " mean_len=", $total_seq_len / $seq_n, " max_len=", $max_seq_len, " min_len=", $min_seq_len, " N50=", $n50, " N90=", $n90, "\n"; 

