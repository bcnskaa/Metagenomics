from __future__ import print_function
from __future__ import division


import glob
import os
import numpy
import sys
from Bio import SeqIO



sys.path.append(os.path.abspath("/home/siukinng/tools/scripts"))


import mg_pipeline
import pick_seq
import merge_fa
import filter_blast_res



"""

"""



"""
Bin scaffolds
"""
def assign_scaffold_to_sid(scaffold_fa_fn, bla_fn):
    print("")
    blast_res = filter_blast_res.import_blast_res(bla_fn)
    
    
    
def estimate_bla_map(blast_res, cutoff_len=40000):
    print("")
    
    
    
    

"""
Calculate the coverage of all sequences in the coverage file
"""
def calculate_mean_coverage(coverage_fn):
    with open(coverage_fn) as IN:
        covs = IN.read().splitlines()
    
    cov_lst = {}    
    for cov in covs:
        
        try:
            [sid, pos, c] = cov.split("\t")
            c = int(c)
            if sid not in cov_lst.keys():
                # len_with_coverage, total_coverage, max_coverage
                cov_lst[sid] = [0, 0, 0]
                
            # len_with_coverage
            cov_lst[sid][0] += 1
            # total_coverage 
            cov_lst[sid][1] += c
            # max_coverage 
            if c > cov_lst[sid][2]:
                cov_lst[sid][2] = c
        except ValueError: 
            print(c)  
 
    return cov_lst
        
        


"""
Helper function for converting the following data:
OTU    Read
OTU_1    Read_1
    Read_2
    Read_3
OTU_2    Read_4
    Read_5
    Read_6
OTU_3    Read_7
OTU_4    Read_8
    Read_9
    Read_10
OTU_5    Read_11
    Read_12
OTU_6    Read_13
    Read_14
    Read_15
    Read_16
    Read_17
OTU_7    Read_18
OTU_8    Read_19
OTU_9    Read_20
OTU_10    Read_21


into:

OTU_1    Read_1    Read_2    Read_3        
OTU_2    Read_4    Read_5    Read_6        
OTU_3    Read_7                
OTU_4    Read_8    Read_9            
OTU_5    Read_11    Read_12            
OTU_6    Read_13    Read_14    Read_15    Read_16    Read_17
OTU_7    Read_18                
OTU_8    Read_19                
OTU_9    Read_20                
OTU_10    Read_21                

"""
def short(ifn, ofn=None):
    if ofn is None:
        ofn = ifn + ".tsv"
    
    with open(ifn) as IN:
        dat = IN.read().splitlines()
    
    cut_id = ""
    lst = {}
    for d in dat:
        id = d.split("\t")[0]
        v = d.split("\t")[1]
        if len(id) == 0:
            id = cut_id
        if id not in lst.keys():
            lst[id] = []
        lst[id].append(v)
        
    with open(ofn, "w") as OUT:
        for id in lst.keys():
            OUT.write(id + "\t" + "\t".join(lst[id]) + "\n")       
    return ofn
    
    
    
def print_msg(msg):
    mg_pipeline.print_status(msg)

