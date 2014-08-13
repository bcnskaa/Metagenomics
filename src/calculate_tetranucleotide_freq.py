from __future__ import division
import itertools
import sys
from itertools import islice
from random import randint

# Require biopython
# To run this script, path to biopython has to be included in PYTHONPATH
from Bio import SeqIO


# Input: A file containing one or more fasta sequences
# Output: A file "tetra_freq.table" which is in tabular format
def calculate_freq_from_fasta(fasta_fn, outtbl_fn="tetra_freq.table"):
    print "[ calculate_freq_from_fasta ] ", "Reading sequences from ", fasta_fn
    
    seqDB = list(SeqIO.parse(fasta_fn, "fasta"))
    
    print len(seqDB), "sequences found"

    freq_tables = {}
    for seq in seqDB:
        print "parsing", seq.id, "(", len(seq.seq) ,"bp)"
        freq_tables[seq.id] = normalize(calculate_freq_table(str(seq.seq)))
    
    #freq_tables = {seq.id: calculate_freq_table(str(seq.seq)) for seq in SeqIO.parse(fasta_fn, "fasta")} 
    export_freq_table(outtbl_fn, freq_tables)
    

# Normalize freq_table by its total sum
def normalize(freq_table):
    total = sum([v for v in freq_table.values()])
    normalized_freq_table = {k:v/total for k,v in freq_table.items()}
    
    return normalized_freq_table


# Export a frequency table into tabular format to an out file
def export_freq_table(outfn, freq_table_dict):
    patterns = list(freq_table_dict[list(freq_table_dict.keys())[0]].keys())
    patterns.sort()
    
    with open(outfn, "w") as OUT:
        # Header 
        OUT.write("\t".join(["ID", "\t".join(patterns), "\n"]))
        
        for key, table in freq_table_dict.items(): 
            values = "\t".join([str(table[k]) if k in table.keys() else str(0) for k in patterns])
            OUT.write("\t".join([key, values, "\n"]))
            #[tt[key] for key in sorted(tt.keys())]


# Given a string (symbol: A, C, G, T) and length of kmer, we calculate occurences of each kmer in
# the string and return a dictionary object enlisting all kmer and their count. One character sliding window is used.
def calculate_freq_table(sequence, pattern_list=['A', 'C', 'G', 'T'], kmer_len=4):
    # Generate a list of k-mer
    seq_tuple = [sequence[i:i + kmer_len] for i in range(1, len(sequence) - (kmer_len - 1))]
    
    # Generate a frequency table through an iterable and list comprehension
    #freq_table = {p : seq_tuple.count(p) for p in ["".join(l) for l in list(itertools.product(pattern_list, repeat=kmer_len))]}
    #patterns = ["".join(l) for l in list(itertools.product(pattern_list, repeat=kmer_len))]
    
    patterns = get_patterns(pattern_list, kmer_len)
    freq_table = {p:seq_tuple.count(p) for p in set(patterns)}
    
    return freq_table


# Construct a list of ALL possible kmers of length (default=4)
def get_patterns(pattern_list=['A', 'C', 'G', 'T'], kmer_len=4):  
    return ["".join(l) for l in list(itertools.product(pattern_list, repeat=kmer_len))] 


# Generate a random sequence of length (default=1000bp) for test-bench code
def generate_random_seq(length=10000, base = ["A", "T", "C", "G"]):
    base_len = len(base) - 1
    seq = "".join([base[randint(0, base_len)] for i in range(1, length)])
    
    return seq


# Test-bench codes
def testbench():
    seqs = {"".join([str(randint(0,9)) for i in range(1,10)]):generate_random_seq() for i in range(1, 10)}
    tbl = {k:calculate_freq_table(seq) for k, seq in seqs.items()}


# Main 
def main(argv):
    calculate_freq_from_fasta(argv[0])

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    