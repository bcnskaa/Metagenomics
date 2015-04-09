from __future__ import print_function
from __future__ import division

from Bio import Entrez
from urllib2 import HTTPError
import time

import glob
import os
import numpy
import sys
import random
from Bio import SeqIO
from collections import Counter 

sys.path.append(os.path.abspath("/home/siukinng/tools/scripts"))

import merge_fa


def process_gi(gi_fn, out_fn=None, resuming_from_fa_fn=None, chunk_size=200):
#def process_gi(gi_fn, out_fn=None, resuming_from_fa_fn=None, chunk_size=200):
    with open(gi_fn) as IN:
        gi_list = IN.read().splitlines()
    
    if resuming_from_fa_fn is not None:
        print("Resuming from " + resuming_from_fa_fn)
        seqs = SeqIO.index(resuming_from_fa_fn, "fasta")
        excluded_gi_list = list(seqs.keys())
        excluded_gi_list = [l.split("|")[1] for l in excluded_gi_list]
        
        tmp_gi_list = gi_list + excluded_gi_list
        gi_counts = Counter(tmp_gi_list)
        
        resume_gi_list = [gi for gi in gi_counts.keys() if gi_counts[gi] == 1]
        print(str(len(gi_list) - len(resume_gi_list)) + " are excluded.")
        
        gi_list = resume_gi_list
    
    
    gi_chunks = [gi_list[x:x+chunk_size] for x in xrange(0, len(gi_list), chunk_size)]
   
    if out_fn is None:
        out_fn = gi_fn + ".fa"
    OUT = open(out_fn, "w")   
    
    export_n = 0
    i = 0
    all_seqs = {}
    for gi_chunk in gi_chunks:
        i = i + 1
        print("Processing chunk " + str(i) + " out of " + str(len(gi_chunks)))
        seqs = extract(gi_chunk)
        if len(seqs) > 0:
            all_seqs.update(seqs)
            for k in seqs.keys():
                merge_fa.export_fasta(k, all_seqs[k], OUT=OUT)
                export_n = export_n + 1
        wait_time = random.randint(3, 20)
        print("Waiting " + str(wait_time) + " sec.")
        time.sleep(wait_time)


   
#     for k in all_seqs.keys():
#         merge_fa.export_fasta(k, all_seqs[k], OUT=OUT)
#         export_n = export_n + 1
    OUT.close()
    
    print(str(export_n) + " sequences exported.")
    
    return all_seqs
    
    

def extract(gi_list):
    Entrez.email = "ppr@na.edu"
    
    print("Extracting " + str(len(gi_list)) + " gi sequences from NCBI")
    
    request = Entrez.epost("Protein", id=",".join(gi_list))
    res = Entrez.read(request)
    we = res["WebEnv"]
    query = res["QueryKey"]
    
    hl = Entrez.efetch(db="protein", retmode="text", rettype="fasta", webenv=we, query_key=query)
    
    #hl = Entrez.efetch(db="Protein", retmode="xml", webenv=we, query_key=query)
    #h = Entrez.efetch(db="Protein", id=line, retmode="fasta")
    
    
    seqs = {}
    for record in SeqIO.parse(hl, "fasta"):
    #records = SeqIO.read(hl, "fasta")
        #print(record.description + ":" + str(len(record.seq)))
        seqs[record.description] = str(record.seq)
        
#     seqs = {}
#     for r in Entrez.parse(hl):
#         try:
#             gi = r
#             print(r)
#             #gi = [x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1]
#             #gi = [x for x in r['GBSeq_other-seqids'] if "gi" in x][0]
#         except ValueError:
#             gi = None
#         
#         if gi is not None:
#             seqs[gi] = r["GBSeq_sequence"]
    
    print(str(len(seqs)) + " sequences extracted.")
   
    return seqs

