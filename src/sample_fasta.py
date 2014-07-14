import os.path
import sys
import getopt
import numpy

# Main 
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:r:",["fasta_infile=","sample_outfile=", ""])
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i", "--fasta_infile"):
            infn = arg
        elif opt in ("-o", "--sample_outfile"):
            outfn = arg
        elif opt in ("-r", "--sample_outfile"):
            outfn = arg
        else:
            print_usage()
            sys.exit(0)
        
    fasta_objs = import_fasta(infn)

# Extract fasta     
def import_fasta(infile):
    fasta_objs = {}
    
    with open(infile, "r") as ins:
        processed_n = 0
        
        cur_header = None
        # read the file handler
        while line in ins:
            line.replace("\n", "")
            if line.startswith(">"):
                cur_header = line
                
                # Initialize a new fasta word
                fasta_objs[cur_header] = ""
                processed_n = processed_n + 1
            else:
                # Make sure that we meet a header before
                if cur_header in locals():
                    fasta_objs[cur_header] = fasta_objs[cur_header] + line
                
    return fasta_objs


def do_sample(fasta_objs, resample_n, resample_seq_size):
    resample_fasta_objs = {}
    
    

def print_usage():
    print("Usage:")
    print("python sample_fasta.py -i FASTA_INFILE -o RESAMPLE_OUTFILE -")
    

# Invoke the main function
if __name__ == "__main__":
    
    main(sys.argv[1:])
    
