#!/usr/bin/perl -w

# 2014 SKWoolf bcnskaa AT gmail DOT com

use lib "/home/siukinng/tools/lib/PerlLib/lib/perl5/site_perl/5.8.8";
use lib "/home/siukinng/tools/lib/PerlLib";
use Bio::DB::EUtilities;
use Bio::SeqIO;
use POSIX;

my $fn = "CAZy_id.lst";

open(IN, $fn);


my %cazy_list;
while (<IN>) {
        chomp;
        if(length($_) < 4) {
                #print $_, "\n";
                next;
        }

        @items = split("\t", $_);
        #if(length($items[4]) > 1) {
        if(!exists $cazy_list{$items[4]}) {
                #print $cazy_list{$items[4]}, " existed: ", $_, "\n";
                $cazy_list{$items[4]} = $_;
        }
        #print join(":", split("\t", $_)), "\n";
        #print $items[4], "\n";
}
close (IN);



my @all_acc_n = keys %cazy_list;
$all_acc_n_size = scalar @all_acc_n;
print "Size=", $all_acc_n_size, "\n";

#
#$step = 10000;
##print ceil($all_acc_n_size / $step), "\n";
#my @gis;
#for(my $i=0; $i < ceil($all_acc_n_size / $step); $i++)
#{
#        print "Processing ", $i * $step, "-", ($i + 1) * $step;
#
#        my @acc_n = @all_acc_n[($i * $step) .. (($i + 1) * $step)];
#        #print join(",", @acc_n), "\n";
#        my $query = Bio::DB::EUtilities->new(-eutil => 'efetch', -db => 'protein', -id => \@acc_n, -email => 'mymail@foo.bar', -rettype => 'gi');
#
#        my @res = split(m{\n},$query->get_Response->content);
#        print ": ", scalar @res, "\n";
#        push @gis, @res;
#   
#}
#
#open(OUT, ">" . $fn . ".out");
#foreach $gi (@gis)
#{
#        print OUT $gi . "\n";
#}
#close(OUT);


############### Retrieve protein sequences

my $faa_outfn = $fn . ".retrieved.faa";
open(OUT, ">" . $faa_outfn);

# Character table for generating random string
my @chars = ("A".."Z", "a".."z", 0..9);

$step = 1000;
#print ceil($all_acc_n_size / $step), "\n";

for(my $i=0; $i < ceil($all_acc_n_size / $step); $i++)
{
		my $tmp_file;
		$tmp_file .= $chars[rand @chars] for 1..9;
		$tmp_file .= ".tmp";
		
        print "Processing ", $i * $step, "-", ($i + 1) * $step . "(" . $tmp_file . ")";

        my @acc_n = @all_acc_n[($i * $step) .. (($i + 1) * $step)];
        #print join(",", @acc_n), "\n";
        my $query = Bio::DB::EUtilities->new(-eutil => 'efetch', -db => 'protein', -id => \@acc_n, -email => 'mymail@foo.bar', -rettype => 'fasta');
		$query->get_Response(-file => $tmp_file);
		
		my $seqs = Bio::SeqIO->new(-file => $tmp_file,  -format => 'fasta');
		
		while(my $seq = $seqs->next_seq)
		{
        	print OUT ">" . $seq->id . "\n";
        	print OUT $seq->seq . "\n\n";
		}

        unlink $tmp_file;
}

close(OUT);



