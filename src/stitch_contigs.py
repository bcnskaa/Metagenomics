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
import re

# To run this script, path to biopython libraries has to be included in PYTHONPATH
from Bio import SeqIO
from Bio import SeqUtils
from Bio import GenBank



# Main 
def main(argv):    
    gb_fns = [os.listdir(f) for f in os.listdir(argv[0]) if os.path.isdir(f)]

    for gb_fn in gb_fns:
        gb_records = import_genbank(gb_fn)
        linking_map = construct_linking_map(gb_records)
        
        

"""
 COG HMMsearch  
"""
def generate_locus_coords_for_reference_genome(gb_fn, dom_dir):
    gb_records = import_genbank(gb_fn)
    locus_map = construct_locus_map(gb_records)
    
    dom_list = import_doms(dom_dir)
    locus_coords = sort_map(dom_list, locus_map)
    
    return locus_coords


    
"""
 Contigs
"""
def generate_locus_coords_from_contigs(dom_dir, locus_coords):
    dom_list = import_doms(dom_dir)
    
    if len(dom_list) == 0:
        return []
    
    dom_list = {dom[1]:((dom[0][::-1]).split("_",1)[1])[::-1] for dom in dom_list}
    
    sorted_locus_coords = sorted(locus_coords.iteritems(), key=operator.itemgetter(1))
    
    stitched_contig_list = []
    for i, k in enumerate(sorted_locus_coords):
        COG_id = k[0]
        if COG_id in dom_list.keys():
            stitched_contig_list.append(dom_list[COG_id])
        else:
            print(COG_id + " is missing.")
   
    return locus_coords



"""
 
"""





"""
 Based on gb_records, we generate a list of locus tags
"""
def construct_locus_map(gb_records):
    processed_n = 0
    #linking_map = {}
    
    locus_map = []
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
            locus_map.append([locus_tag, spos, epos, gi, gene_id, protein_id, product, translation])
            
    print "Processed item=" + str(processed_n) 
    #linking_map = sorted(linking_map, key=lambda v:v[0], reverse=True)
    locus_map = sorted(linking_map, key=lambda v:v[1], reverse=True)
    
    return locus_map



# 
def get_map_by_gi(locus_map, gi):
    gi_map = {map[3]:map for map in locus_map}
    
    if gi in gi_map.keys():
        return(gi_map[gi])
    else:
        return None



# 
def sort_map(dom_list, locus_map):
    locus_coords = { dom[1]:get_map_by_gi(locus_map, dom[0])[1] for dom in dom_list if get_map_by_gi(locus_map, dom[0]) is not None }
    
    return locus_coords
  
    

"""
    Import data in Genbank format. The file expects to contain only one genbank object (one accession number)
"""
def import_genbank(gb_fn):
    gb_records = SeqIO.read(open(gb_fn,"r"), "genbank")
    
    #with open(gb_fn) as IN:
    #    gb_records = GenBank.read(IN)
        
    #IN.close()
    return gb_records
    
    

"""
  
"""
def extract_gi(query_id):
    items = query_id.split("|")
    return(items[1])
    
    
    
# import_doms(glob.glob("*")[0])
def import_doms(dom_dir, extract_gi=False):
    print("Importing from " + dom_dir)
    
    dom_fns = glob.glob(dom_dir + "/*.dom")
    dom_list = []
    if len(dom_fns) > 0:
        for i, dom_fn in enumerate(dom_fns):
            print("Processing " + dom_fn)
            doms = import_dom(dom_fn, extract_gi=extract_gi)
            if len(doms) > 0:
                dom_list.append(doms)
        
    return(dom_list)



# Import dom file
def import_dom(dom_fn, extract_gi=False):
    with open(dom_fn, "r") as IN:
        lines = [l for l in IN.read().splitlines() if not l.startswith("#")]
    
    #print(dom_fn + ": " + str(len(lines)))
    
    if len(lines) == 0:
        return []
    
    #print(lines[0])
    
    # Extract the query id
    #items = lines[0].split("\t")
    
    # Bug may exist
    items = re.split(r"\s+", lines[0])
    
    #print(len(items))
    
    if extract_gi:
        query_id = extract_gi(items[0])
    else:
        query_id = items[0]
         
    COG_id = items[3]
    print(query_id + "\t" + COG_id)  
     
    return [query_id, COG_id]
    
    
    
# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])