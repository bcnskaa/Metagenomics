import os
import sys
import getopt
import csv
import operator
from collections import defaultdict
import inspect
import subprocess
from subprocess import call
import glob
import os.path
import time
from time import strftime
from datetime import datetime
import math
import numpy

# To run this script, path to biopython libraries has to be included in PYTHONPATH
from Bio import SeqIO
from Bio import SeqUtils
from Bio import GenBank

# Main 
def main(argv):
    
        
    

"""
    Using CDSs as anchors, to generate a list
"""
def construct_linking_map(gb_records):
    processed_n = 0
    #linking_map = {}
    linking_map = []
    for feature in gb_records.features:
        if feature.type == "CDS":
            if "gene" in feature.qualifiers.keys():
                gene_id = feature.qualifiers["gene"]
            else:
                gene_id = "-"
                
            locus_tag = feature.qualifiers["locus_tag"][0]
            spos = int(feature.location.start)
            epos = int(feature.location.end)
            protein_id = feature.qualifiers["protein_id"]
            product = feature.qualifiers["product"]
            translation = feature.qualifiers["translation"]
            gi = [xref.replace("GI:", "") for xref in feature.qualifiers["db_xref"] if "GI:" in xref]
            if len(gi) != 0:
                gi = gi[0]
            
            processed_n += 1
            #linking_map[locus_tag] = [spos, epos, gi, gene_id, protein_id, product, translation]
            linking_map.append([locus_tag, spos, epos, gi, gene_id, protein_id, product, translation])
            
       
    print "Processed item=" + str(processed_n) 
    #linking_map = sorted(linking_map, key=lambda v:v[0], reverse=True)
    linking_map = sorted(linking_map, key=lambda v:v[1], reverse=True)
        
    return linking_map 
        


"""
    Import data in Genbank format. The file expects to contain only one genbank object (one accession number)
"""
def import_genbank(gb_fn):
    
    gb_records = SeqIO.read(open(gb_fn,"r"), "genbank")
    
    #with open(gb_fn) as IN:
    #    gb_records = GenBank.read(IN)
        
    #IN.close()

    return gb_records
    
    
    
# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])