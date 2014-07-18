import os.path
import sys
import getopt
import csv
import operator
from collections import defaultdict





# Main 
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)
        
    infn = None
    outfn = None
    lib_path = None

    
    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            infn = arg
        elif opt in ("-o", "--outfile"):
            outfn = arg
        elif opt == "-p":
            lib_path = arg
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)

    # Make sure infilename and outfilename are specified
    if (infn is None and outfn is None and lib_path is None):
        print("Names of input files or path to NCBI bacterial library are missing.\n")
        print_usage()
        sys.exit(0)
    
    
    
    ncbi_species = prepare_ncbi_bacterial_db(lib_path)
    
    
    
    
    sids = import_sids(infn)
    
    
    species_ids = [process_sid(id) for id in sids]
    
    
    


# Import subject ids from .sid file 
def import_sids(infn):
    with open(infn) as f:
        sids = f.read().splitlines()
    return sids

    
# Clear up sid
def process_sid(sid):
    sid.replace(" 16S ribosomal RNA, partial sequence", "")
    sid.replace(" 16S ribosomal RNA, complete sequence", "")
    
    return sid


# 
def prepare_ncbi_bacterial_db(path):
    print("Indexing ncbi bacterial db...")
    
    # Survey the name of available reference genomes
    #ncbi_species = os.popen("ls " + path).read().splitlines();
    # Survey the name of available reference genomes
    ncbi_species = os.popen("ls " + path + "/*/*.fna").read().splitlines();
    
    if len(ncbi_species) == 0:
        return None
    
    (species, seq_names) = 
    
    # tidy up ids
    # Only keep the information of genus and species, 
    ncbi_species = [" ".join(s.split("_")[0:2]) for s in ncbi_species if len(s.split("_")) > 1]




    print "Number of NCBI species =", len(ncbi_species)

    return ncbi_species
    


# This function checks if the given species id has a completed genome in NCBI Bacterial DB, 
# otherwise it suggests a closest sibling species. If no silbling species is available, None
# is returned.
def infer_organisms(species_id):  
    print("Inferring...")
    
    


# Print the usage of this script
def print_usage():
    print("This simple script generates a reference genomes library of prokaryotic species by referring to the subject ids.")
    print(" ")
    print("Usage:")
    print("  python prepare_ref_genomes.py -i SID-INFILE -o LIB-OUTFILE [-p PATH-TO-NCBI-BACTERIAL-LIB]")
    print("      -i STRING  A .sid file exported from filter_blast_res.py with -s option.")
    print("                 The subject ids or sids are expect to be")
    print("      -o STRING  Output library file")
    print("      -p STRING  Path to the NCBI Complete Bacterial Library [option]")
    print(" ")
    print(" ") 
    print("Version 0.1")


# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])