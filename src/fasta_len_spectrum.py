import os.path
import sys
import getopt
from collections import defaultdict
from Bio import SeqIO
import numpy


# Main 
def main(argv):
    fasta_fn = argv[0]
    
    cap_val = 300000
    bin_size = 10000
    
    seq_lens = []
    with open(fasta_fn, "r") as infile:
        for seq in SeqIO.parse(infile, "fasta"):
            seq_lens.append(len(seq.seq));
        
    # Get the maximum length    
    max_seq_len = max(seq_lens)
    
#     # Go through all the length
#     for i in range(max_seq_len):
#         if(seq_lens.count(i) > 0):
#             print "%d\t%d" % (i, seq_lens.count(i))
#             
    bins = numpy.arange(0, cap_val, bin_size)
    if cap_val < max_seq_len:
        numpy.append(bins, max_seq_len)  
    freqs, bins = numpy.histogram(seq_lens, bins)
    
    for i, f in enumerate(freqs):
        print(str(bins[i]) + "\t" + str(bins[i+1]) + "\t" + str(f))



# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    
    
