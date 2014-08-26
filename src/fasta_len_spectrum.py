import os.path
import sys
import getopt
from collections import defaultdict
from Bio import SeqIO


# Main 
def main(argv):
    fasta_fn = argv[0]
    
    seq_lens = []
    with open(fasta_fn, "r") as infile:
        for seq in SeqIO.parse(infile, "fasta"):
            seq_lens.append(len(seq.seq));
        
    # Get the maximum length    
    max_seq_len = max(seq_lens)
    
    # Go through all the length
    for i in range(max_seq_len):
        if(seq_lens.count(i) > 0):
            print "%d\t%d" % (i, seq_lens.count(i))



# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    
    
