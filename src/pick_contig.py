import os.path
import sys
import getopt
import csv
import operator
from collections import defaultdict



# Main 
def main(argv):
    
    



# Print the usage of this script 
def print_usage():
    print("A simple python script for picking contigs based on blast result (generated with -m6 option)")
    print(" ")
    print("Usage:")
    print("  python pick_contig.py -i BLAST-RESULT-INFILE -o FILTER-OUTFILE [-q] [-b BITSCORE-CUTOFF] [-l ALIGNMENT-LENGTH-CUTOFF] [-p PERCENTAGE-IDENTITY]")
    print("      -i STRING  A list contains contig names to be picked")
    print("      -d STRING  File contains all contig sequences")
    print("      -o STRING  Output file contains selected contig sequences")
    print(" ")
    print(" ") 
    print("Ver 0.2")



# Print status
def print_status(msg):
    print("[]", msg)