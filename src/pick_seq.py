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
        opts, args = getopt.getopt(argv,"hki:o:d:")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)
    
    lstfn = None
    outfn = None
    fasta_fn = None
    search_key = False

    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i"):
            lstfn = arg
        elif opt in ("-o"):
            outfn = arg
        elif opt == "-d":
            fasta_fn = arg
        elif opt == "-k":
            search_key = True 
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)


    # Make sure everything is provided
    if(lstfn is None or fasta_fn is None):
        print("Missing argument")
        print_usage()
        sys.exit(0)
        
       
    if(outfn is None):
        outfn = lstfn + ".fa"
        
    list = import_list(lstfn)
    
    # Call the Bio.SeqIO
    # Create a dict index on the fasta sequence file
#     fasta_seqs_idx = SeqIO.index(fasta_fn, "fasta")
#     exported_n = 0
#     with open(outfn, "w") as OUT:
#         for id in list:
#             if id in fasta_seqs_idx.keys():
#                 exported_n += 1
#                 SeqIO.write(fasta_seqs_idx[id], OUT, "fasta")
    
    if search_key:
        exported_n = pick_seq_contain_key(list, fasta_fn, outfn)
    else:
        exported_n = pick_seq(list, fasta_fn, outfn)
    
    print "Number of sequence exported: " + str(exported_n)
   


def pick_seq(id_list, fasta_fn, outfn, is_append=False):
    exported_n = 0
    
    if not os.path.isfile(fasta_fn):
        print(fasta_fn + " is not found.")
        return exported_n
    
    fasta_seqs_idx = SeqIO.index(fasta_fn, "fasta")
    if is_append:
        OUT = open(outfn, "a")
    else:
        OUT = open(outfn, "w")

    for id in id_list:
        if id in fasta_seqs_idx.keys():
            exported_n += 1
            SeqIO.write(fasta_seqs_idx[id], OUT, "fasta")
    OUT.close()
    return exported_n




def pick_seq_contain_key(key_list, fasta_fn, outfn, is_append=False, search_descript=True):
    exported_n = 0
    
    if not os.path.isfile(fasta_fn):
        print(fasta_fn + " is not found.")
        return exported_n
    
    fasta_seqs_idx = SeqIO.index(fasta_fn, "fasta")
    if is_append:
        OUT = open(outfn, "a")
    else:
        OUT = open(outfn, "w")

    for seq_id in fasta_seqs_idx.keys():
        desc = fasta_seqs_idx[seq_id].description
        for k in key_list:
            if k in desc:
                exported_n += 1
                fasta_seqs_idx[id].id = fasta_seqs_idx[seq_id].description
                SeqIO.write(fasta_seqs_idx[id], OUT, "fasta")
                break
    OUT.close()
    return exported_n


'''
    Import a list of fasta sequence headers
'''
def import_list(list_fn):
    l = []
    with open(list_fn) as IN:
        l = IN.read().splitlines()
    
    l = [f for f in l if len(f) > 0]
    l = list(set(l))
    return l


# 
# def subsample():
#     import random
#     import os
#     
#     
#     
#     
    

# Print the usage of this script 
def print_usage():
    print("A simple python script for picking sequence based on blast result (generated with -m8 option)")
    print(" ")
    print("Usage:")
    print("  python pick_seq.py -i LIST_INFILE -o OUTFILE -d FASTA-FILE")
    print("      -i STRING  A list contains contig names to be picked")
    print("      -d STRING  File, from which FASTA sequences enlisted in the list file will be picked")
    print("      -o STRING  Picked FASTA sequence will be output to")
    print("      -k         Use the rows in the list file as keyword, and any sequence header containing any keyword will be picked")  
    print(" ")
    print(" ")
    print("Ver 0.00002")


      
# Print status
def print_status(msg):
    print("[]", msg)
    

    
# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    