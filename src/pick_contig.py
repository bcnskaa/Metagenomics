# 2014 SKWoolf bcnskaa AT gmail DOT com

import os.path
import sys
import getopt
from collections import defaultdict

# 
from utils.blastres import *

# Require biopython
# To run this script, path to biopython libraries has to be included in PYTHONPATH
from Bio import SeqIO


# Main 
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:d:b:")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)
        
    
    blastfn = None
    infn = None
    outfn = None
    fastafn = None

    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i"):
            infn = arg
        elif opt in ("-o"):
            outfn = arg
        elif opt == "-d":
            fastafn = arg
        elif opt == "-b":
            blastfn = arg
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)


    # Make sure everything is provided
    if((infn is None and blastfn is None) or outfn is None or fastafn is None):
        print("Missing argument")
        print_usage()
        sys.exit(0)   


    # 
    if infn is not None:
        # Import a list
        list = read_list(infn)
 
 
#     res = ""
#     with open(blastfn, "r") as IN:
#         while res in IN:
#             print "result:", res, "."

    if blastfn is not None:
        print "blast file=", blastfn
        res = ""
        with open(blastfn) as IN:
            while res in IN:
                print "result:", res, "."   
        
        blastres = BlastParser.parse(blastfn)
        print "blastres=", len(blastres)
        list = BlastParser.getSubjectIDs(blastres)
    
    
    print "Number of element in", infn, "=", len(list)
    
    # 
    if len(list) == 0:
        print("The list has no member.")
        sys.exit(0)
    
    # Picking
    pick_fasta(fastafn, outfn, list)
    
    

# 
def pick_fasta(fastafn, outfn, list):
    print "Picking fasta sequence(s) from", fastafn, "based on the list", len(list), "..."
        
    with open(outfn, "w") as outfile:
        for seq in SeqIO.parse(fastafn, "fasta"):
            if seq.id in list:
                print "extracting", seq.id, "..."
                outfile.write(">" + seq.id + "\n")
                outfile.write(seq.seq + "\n\n")
    
    
#     with open(fasta_fn, "r") as infile, open(outfn, "w") as outfile:
#         while line in ins:
#             line.replace("\n", "")
#             
#             # Check if this sequence is on the list
#             if line.startswith(">"):
#                 if line.replace(">", "") in list:
#                     outfile.write(">" + line + "\n")
#                     flag = True
#                 else:
#                     flag = False
#             else:
#                 if cflag:
#                     outfile.write(line + "\n")


# # Read a list of selected fasta id
# def read_qid_from_blast(blastfn):
#     id_list = []
#     
#     # export the gi and name to an outfile
#     with open(outfile, "r") as infile:
#         id_list = [line.replace("\n", "") for line in infile]
#     return id_list

    
# Read a list of selected fasta id
def read_list(listfn):
    id_list = []
    
    try: 
        # export the gi and name to an outfile
        with open(listfn, "r") as infile:
            id_list = [line.replace("\n", "") for line in infile]
    finally:
        infile.close()
    return id_list

  
    
# Print the usage of this script 
def print_usage():
    print("A simple python script for picking contigs based on blast result (generated with -m8 option)")
    print(" ")
    print("Usage:")
    print("  python pick_contig.py -i BLAST-RESULT-INFILE -o FILTER-OUTFILE -d FASTA-SEQ-DB")
    print("      -i STRING  A list contains contig names to be picked")
    print("      -d STRING  File, from which FASTA sequences enlisted in the list file will be picked")
    print("      -o STRING  Picked FASTA sequence will be output to")
    print(" ")
    print(" ") 
    print("Ver 0.00001")



# Print status
def print_status(msg):
    print("[]", msg)
    
    
# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    