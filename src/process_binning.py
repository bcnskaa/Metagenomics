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
If sequence id is in format "XXXX.X", this function extracts the text before "."
"""
def wash_id(ids):
    washed_ids = [id.split(".")[0] for id in ids]
    return washed_ids
    
    
    
def print_msg(msg):
    mg_pipeline.print_status(msg)

