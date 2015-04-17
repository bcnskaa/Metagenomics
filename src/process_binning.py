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
        
        

        



    
    
    
def print_msg(msg):
    mg_pipeline.print_status(msg)

