import os.path
import sys
import getopt
import csv
import operator
from collections import defaultdict

import extract_gi_from_fna


#
# Given a blast result (generated with -m6 option), this script will return you the best hit for every query id
#
verbose = False
def main(argv):
    try:
        # Parse options from command lines
        cmd_opts, cmd_args = getopt.getopt(sys.argv[1:], "hig:v", ["help", "infile=", "gi_list="])
    except getopt.GetoptError as e:
        print str(e)
        print_usage()
        
        sys.exit(1)
    
    
    outfn = None    # output file name
    infn = None     # input file name
    gi_infn = None  # input file of gi list
    
    #verbose = False
    
    # Parse command line
    for opt, arg in cmd_opts:
        if opt == "-v":
            verbose = True
        elif opt in ("-h", "--help"):
            print_usage()
            sys.exit()
        elif opt in ("-o", "--outfile"):
            outfn = arg
        elif opt in ("-g", "--gi_list"):
            gi_infn = arg
            
            if(not os.path.isfile(gi_infn)):
                print(gi_infn + " is not found.")
                sys.exit()                 
        elif opt in ("-i", "--infile"):
            infn = arg
            
            if(not os.path.isfile(infn)):
                print(infn + " is not found.")
                sys.exit()       
        else:
            assert False, "unrecognized option"
        
    res = import_blast_m6_results(infn)
    gi_list = import_gi(gi_infi)
    
    print "Size of result =", len(res)


    # Retrieve best hits
    besthits = report_besthit(res)
    
    #for qid in besthits.keys():
    #    print qid, " = ", besthits[qid][1], " bitscore=", besthits[qid][11]
        
    #sorted_besthits = sorted(besthits.items(), key=lambda (k, v): v[1][11], reverse=True)
    #for besthit in sorted_besthits:
    #    print besthit[0], ":", besthit[1][1], " bitscore=", besthit[1][11]

    for k, v in sorted(besthits.items(), key=lambda e:e[1][11], reverse=True):
        print k, ":", v[1], " bitscore=", v[11]


def import_gi(gi_infilename):
    gi_list = {}
    with open(gi_infilename, "r") as infile:
        for line in infile:
            infos = (line.replace("\n","")).split("\t")
            gi_list[info[0]] = infos[1]
        infile.close()
    return gi_list

# Import the blast -m6 result into a dictionary of lists using query ids as keys
def import_blast_m6_results(infilename):
    if verbose:
        print("Checking if input file is a valid BLAST -m6 format.")
    
    blast_res = defaultdict(list)
    
    # infilename should be a tab-delimited csv file
    with open(infilename, "r") as infile:
        input_data = csv.reader(infile, delimiter='\t')
        
        for row in input_data:
            row[11] = float(row[11])
            blast_res[row[0]].append(row)
        infile.close()
            
    return blast_res


# Report the best hit for every query id 
def report_besthit(blast_m6_results):
    
    blast_besthits = {}
    for qid in blast_m6_results.keys():
        blast_besthits[qid] = sorted(blast_m6_results[qid], key=lambda x: x[11], reverse=True)[0]
    
    return blast_besthits
    

# Print the usage of this script
def print_usage():
    print("Usage:")
    print("python extract_blast_besthit.py -i INFILE -o OUTFILE")
    

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])