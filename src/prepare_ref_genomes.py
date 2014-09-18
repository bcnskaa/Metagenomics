#!/bin/python

# 2014 SKWoolf bcnskaa AT gmail DOT com

import os.path
import sys
import getopt
import csv
import operator
import time
from collections import defaultdict
from subprocess import call


ID_DELIM_CHR = "_"

#NCBI_BACTERIAL_DB_PATH="/home/jiapchen/sk/lib/BacteriaDB/all_fna"
NCBI_BACTERIAL_DB_PATH="/home/siukinng/db/BacteriaDB/all_fna"
DWGSIM_HOME="/home/siukinng/tools/dwgsim"
COVERAGE=10
READ_ERROR_RATE=0.001
READ_LENGTH=90

INFER_SLIM_MODE=False

# Main 
def main(argv):
    global COVERAGE
    global READ_ERROR_RATE
    global READ_LENGTH
    global INFER_SLIM_MODE
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:c:e:l:d:rs")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)
        
    infn = None
    outfn_prefix = None
    lib_path = NCBI_BACTERIAL_DB_PATH
    build_reference = False
    output_dir = None
    
    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            infn = arg
        elif opt in ("-o", "--outfile"):
            outfn_prefix = arg
        elif opt == "-p":
            lib_path = arg
        elif opt == "-c":
            COVERAGE = int(arg)
            #print("New coverage:", COVERAGE)
        elif opt == "-l":
            READ_LENGTH = int(arg)
        elif opt == "-e":
            READ_ERROR_RATE = float(arg)
        elif opt == "-r":
            build_reference = True    
        elif opt == "-d":
            output_dir = arg   
        elif opt == "-s":  
            INFER_SLIM_MODE = True  
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)

            
    if output_dir is None:
        output_dir = outfn_prefix


    # Make sure infilename and outfilename are specified
    #if (infn is None and outfn is None and lib_path is None):
    if (infn is None and outfn_prefix is None):
        print("Names of input files or path to NCBI bacterial library are missing.\n")
        print_usage()
        sys.exit(0)

    
    ncbi_bacterial_db = prepare_ncbi_bacterial_db(lib_path)
    
    sids = import_sids(infn)
    
    species_ids = [process_sid(id) for id in sids]
    
    # See if we can find the sequences matched with given species id
    bacterial_db_ids = []
    for species_id in species_ids:
        # Try to infer matched with or similar to the given species
        ids = infer_organisms(species_id, ncbi_bacterial_db, slim_mode=INFER_SLIM_MODE)
        if ids is not None:
            if len(ids) > 0:
                #print "Id found=", ids
                bacterial_db_ids.extend(ids)

    # length before making it unique 
    print "before", len(bacterial_db_ids)

    # Make a unique db_ids
    bacterial_db_ids = set(bacterial_db_ids)
    
    # length after making it unique 
    print "after", len(bacterial_db_ids)  
    
    for id in bacterial_db_ids:
        print id

    if build_reference:
        print "Constructing reference genome library..."
        contruct_referenece_genome_library(bacterial_db_ids, ncbi_bacterial_db, outfn_prefix)
    else:
        seq_n = 0
        
        # Create the output folder if it doesn't exist
        if not os.path.isdir(output_dir):
            print "Creating " + output_dir
            os.system("mkdir " + output_dir)
        
        # Copy sequences of reference genomes
        for id in bacterial_db_ids:  
            fn = output_dir + "/" + id + ".fna"
            paths_to_genome = ncbi_bacterial_db[id]
            
            print "Number of sequences to be processed for " + id + ":", len(paths_to_genome)
            
            for seq_fn in paths_to_genome:
                print "processing", seq_fn
                os.system("cat " + seq_fn + " >> " + fn)    
                seq_n += 1
                
        print "Sequence processed:", seq_n      
        
                


# Build read library from recruited reference genomes
def contruct_referenece_genome_library(selected_ids, db, output_prefix):
    
    # Temporary folders for storing reference genomes
    TMP_DIR = "./" + str(time.time()) + ".tmp"
    
    call(["mkdir", TMP_DIR])
    
    # Go through all ids
    for id in selected_ids:    
        fn = TMP_DIR + "/" + id + ".fna"
        outfn_prefix = TMP_DIR + "/" + id
        
        # Get paths to sequences 
        seq_paths = db[id]
        
        print "Number of sequences to be processed:", len(seq_paths)
        
        # Pool all sequences of a "species" together
        for seq_fn in seq_paths:
            #call(["cat", seq_fn, ">>", fn])
            print "processing", seq_fn
            os.system("cat " + seq_fn + " >> " + fn)
        
        # Construct a command
        cmd = DWGSIM_HOME + "/dwgsim -1 " + str(READ_LENGTH) + " -2 " + str(READ_LENGTH) + " -e " + str(READ_ERROR_RATE) + " -E " + str(READ_ERROR_RATE) + " -C " + str(COVERAGE) + " " + fn + " " + outfn_prefix 
        print "Executing command =", cmd
        os.system(cmd)
       
    # Mix all the reads into a single file
    print "Generating read1 data to ", output_prefix + "_1.fq", "..."
    cmd = "cat " + TMP_DIR + "/*.read1.fastq >> " + output_prefix + "_1.fq"
    os.system(cmd)
    print "Generating read2 data to ", output_prefix + "_2.fq", "..."
    cmd = "cat " + TMP_DIR + "/*.read2.fastq >> " + output_prefix + "_2.fq"
    os.system(cmd)
    
    # Mix all referernce genomes
    print "Concatenating referernce genomes to", output_prefix + "_all_refgenomes.fna"
    cmd = "cat " + TMP_DIR + "/*.fna >> " + output_prefix + "_all_refgenomes.fna"
    os.system(cmd)
    
    # Rename the temporary folder to output_prefix
    cmd = "mv " + TMP_DIR + " " + output_prefix
    os.system(cmd)



# Import subject ids from .sid file 
def import_sids(infn):
    with open(infn) as f:
        sids = f.read().splitlines()
    return sids

    
# Clear up sid
def process_sid(sid):
    sid = sid.replace(" 16S ribosomal RNA, partial sequence", "")
    sid = sid.replace(" 16S ribosomal RNA, complete sequence", "")
    
    return sid



# "Path" is a string indicating the location of the NCBI Bacterial Database
# Sequences should be stored under the folder with the species names
def prepare_ncbi_bacterial_db(path):
    print("Indexing ncbi bacterial db...")
    
    # Survey the name of available reference genomes
    #ncbi_species = os.popen("ls " + path).read().splitlines();
    
    # Survey the name of available reference genomes
    ncbi_species = os.popen("ls " + path + "/*/*.fna").read().splitlines();
    
    # Check if sequences exists
    if len(ncbi_species) == 0:
        return None
    
    # We need the species name and their genomic sequences. Split the name and sequences from the file path
    #(species, seq_names) = zip(*([s.split("/")[-2, -1] for s in ncbi_species]))
    
    
    # Construct the bacterial database
    ncbi_bacterial_db = defaultdict(list)
    for s in ncbi_species:
        # Make sure that the folder name is biologically meaningful, or skip otherwise
        if not len(s.split("_")) > 0:
            continue
        
        #(species, seq_names) = zip(([s.plit("/")[i] for i in [-2, -1]]))
        
        # Tidy up the id as we're only interested in the name of species, and not the strain info
        seq_path = s
        
        s = s.split("/")[-2]
        genus = s.split("_")[0]
        species = s.split("_")[1]
             
        #print species, "=", seq_name
        ncbi_bacterial_db[genus + ID_DELIM_CHR + species].append(seq_path)
        
    
    # tidy up ids
    # Only keep the information of genus and species, 
    #ncbi_species = [" ".join(s.split("_")[0:2]) for s in ncbi_species if len(s.split("_")) > 1]

    print "Number of NCBI species =", len(ncbi_species)

    #print "Number of sequences for Zymomonas mobilis =", len(ncbi_bacterial_db["Zymomonas_mobilis"])
    #for k, v in ncbi_bacterial_db.items():
    #    print k, "=", len(v)


    total = sum([len(v) for k, v in ncbi_bacterial_db.items()])
    print "Total number of sequences=", total
    
    return ncbi_bacterial_db
    


# This function checks if the given species id has a completed genome in the NCBI Bacterial DB, 
# otherwise it suggests a closest sibling species. If no silbling species is available, None
# is returned.
def infer_organisms(species_id, ncbi_bacterial_db, slim_mode=False):    
    #print "Matching", species_id, "to bacterial db"

    # Flagging our founding
    cflag_found = False
    found_keys = []

    # first check if the species_id can be used directly
    if species_id in ncbi_bacterial_db:
        cflag_found = True
        found_keys.append(species_id)
        return found_keys
    
    # Make sure the specie_id contains delimiting characters
    if not cflag_found and not len(species_id.split(" ")) > 1:
        print "No match is found for ", species_id
        return []
 
 
    # Extract the first two words of the species_id
    genus = species_id.split(" ")[0]
    species = species_id.split(" ")[1]
           
    # Genus + Species name
    if not cflag_found:
        print "Searching", genus + ID_DELIM_CHR + species
    
        if genus + ID_DELIM_CHR + species in ncbi_bacterial_db:
            cflag_found = True
            found_keys.append(genus + ID_DELIM_CHR + species)
            return found_keys
        
    # Genus only
    if not cflag_found:
        name = genus
        print "Searching", genus
        
        for key in ncbi_bacterial_db:
            if key.startswith(name + ID_DELIM_CHR):
                cflag_found = True
                found_keys.append(key)
                
                # If in slim mode, we export one similar item only
                if slim_mode:
                    return found_keys
                
            
    # Print the matched name
    #for k in found_keys:
    #    print "key=", k
    
    return found_keys

    


# Print the usage of this script
def print_usage():
    print("This simple script generates a reference genomes read library of prokaryotic species based on given subject ids.")
    print(" ")
    print("Usage:")
    print("  python prepare_ref_genomes.py -i SID-INFILE -o LIB-OUTFILE [-s] [-r] [-d output-dir] [-p PATH-TO-NCBI-BACTERIAL-LIB] [-l READ-LEN] [-e READ-ERROR-RATE] [-c COVERAGE]")
    print("      -i STRING  A .sid file exported from filter_blast_res.py with -s option.")
    print("                 The subject ids or sids are expect in 'genus+species+strain' format. (required)")
    print("      -o STRING  Prefix of an output library file. (required)")
    print("      -d STRING  Output directory [optional]")
    print("      -s         Slim selection, only one is picked from a list of strains of same species [optional]")
    print("      -r         Build a read library of reference genomes [optional, default=False]")
    print("      -p STRING  Path to the NCBI Complete Bacterial Library [optional, default=" + NCBI_BACTERIAL_DB_PATH + "]")
    print("      -l INT     Read length [optional, default=" + str(READ_LENGTH) + "]")
    print("      -e FLOAT   Read error rate [optional, default=" + str(READ_ERROR_RATE) + "]")
    print("      -c INT     Coverage [optional, default=" + str(COVERAGE) + "]")
   
    print(" ")
    print(" ") 
    print("Version 0.1")


# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])