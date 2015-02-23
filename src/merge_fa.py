# 2014 SKWoolf bcnskaa AT gmail DOT com
import os.path
import sys
import getopt
from collections import defaultdict

# Require biopython
# To run this script, path to biopython libraries has to be included in PYTHONPATH
from Bio import SeqIO


# Main 
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)

    fasta_char_per_line = 70
    header = None
    fasta_ifn = None
    fasta_ofn = None

    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i"):
            fasta_ifn = arg
        elif opt in ("-o"):
            fasta_ofn = arg
        elif opt == "-n":
            header = arg
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)


    # Make sure everything is provided
    if(fasta_ifn is None or fasta_ofn is None or header is None):
        print("Missing argument")
        print_usage()
        sys.exit(0)
        
       
    if(fasta_ofn is None):
        fasta_ofn = fasta_ifn + ".merged.fa"
        

    # Call the Bio.SeqIO
    # Create a dict index on the fasta sequence file
    seqs = SeqIO.index(fasta_ifn, "fasta")
    merged_seq = ""
    exported_n = 0
    for seq_id in seqs.keys():
        merged_seq = merged_seq + str(seqs[seq_id].seq) 
        exported_n += 1
        
    seq_lines = [merged_seq[i:i+fasta_char_per_line] for i in range(0, len(merged_seq), fasta_char_per_line)]
    with open(fasta_ofn, "w") as OUT:
        OUT.write(">" + header + "\n")
        OUT.write("\n".join(seq_lines))
        OUT.write("\n")

    print("Number of sequence merged: " + str(exported_n) + " (" + str(len(merged_seq)) + "bp)")

    
    

# Print the usage of this script 
def print_usage():
    print("  python merge_fa.py -i FASTA-INFILE -o FASTA-OUTFILE -n HEADER")
    print("      -i STRING  Fasta infile")
    print("      -o STRING  Merged outfile")
    print("      -n STRING  Header of merged fasta sequence")
    print(" ")



    
# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    