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







def calculate_mean_coverage(coverage_fn):
    with open(coverage_fn) as IN:
        covs = IN.read().splitlines()
    
    cov_lst = {}
    for cov in covs:
        [sid, pos, c] = cov.split("\t")
        
        if sid not in cov_list.keys():
            



    
    
    
def print_msg(msg):
    mg_pipeline.print_status(msg)

