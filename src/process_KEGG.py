##
# KEGG REST API:
# http://www.kegg.jp/kegg/docs/keggapi.html
#
##
from __future__ import print_function
from __future__ import division

import glob
from itertools import izip
import gc
import os
from collections import Counter
import operator



"""
Summary of the procedure:

1. Obtain the genome id at http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+genome:T00308

2. Specify the gene id: moth_0450

3. Obtain the information of the gene at http://www.genome.jp/dbget-bin/www_bget?mta:Moth_0450
        mta = genome abbreviation
        Moth_0450 = gene id
        
4. Parse the HTML for the tag: http://www.genome.jp/dbget-bin/www_bget?pfam:
   and extract the Pfam name.   
        
5. Import the Pfam HMM info from /home/siukinng/db/Markers/Pfam/Pfam.hmm.info

6. Map the extracted Pfam name to the HMM info

7. Export the data



Usage:

import process_KEGG

genes = ['Moth_1191', 'Moth_1192', 'Moth_1193', 'Moth_1194', 'Moth_1195', 'Moth_1196']
genes2pfam_tbl = process_KEGG.retrieve_gene_info(genes, "mta")
process_KEGG.export_gene2pfam_tbl(genes2pfam_tbl, "methylene-THF-reductase.info")

"""
def retrieve_gene_info(gene_list, genome_abbrev):
    print("Retrieving " + str(len(gene_list)) + " KEGG genes")
    
    pfam_info = import_pfam_info()
    
    gene2pfam_tbl = {}
    for gene in gene_list:
        print("Processing " + gene)
        cmd = "wget -o log http://www.genome.jp/dbget-bin/www_bget?" + genome_abbrev + ":" + gene + " -O " + gene
        os.system(cmd)
        
        # Parse the Pfam
        cmd = "grep -o \'pfam:[^\"]*\' " + gene + " | sed \"s/pfam://\""
        pfam_domains = os.popen(cmd).read().split("\n")
        pfam_domains = [p for p in pfam_domains if len(p) > 0]
        
        gene2pfam_tbl[gene] = []
        for pfam_domain in pfam_domains:
            gene2pfam_tbl[gene].append([pfam_domain] + pfam_info[pfam_domain])
    
    return gene2pfam_tbl
        




def import_pfam_info(pfam_info_fn="/home/siukinng/db/Markers/Pfam/Pfam.hmm.info"):
    print("Import Pfam info form " + pfam_info_fn)
    pfam_info = {}
    with open(pfam_info_fn) as IN:
        pfam_info = IN.read().splitlines()
    pfam_info = {p.split("\t")[0]: p.split("\t")[1].split("|") for p in pfam_info}
    
    print(str(len(pfam_info)) + " Pfam imported.")
    return pfam_info
    


def export_gene2pfam_tbl(gene2pfam_tbl, out_fn):
    OUT = open(out_fn, "w")
    
    for gene in gene2pfam_tbl.keys():
        OUT.write(gene + "\t" + "|".join([p[0] for p in gene2pfam_tbl[gene]]) + "\t" + "|".join([p[1] for p in gene2pfam_tbl[gene]]) + "\t" + "|".join([p[2] for p in gene2pfam_tbl[gene]]) + "\n") 

    OUT.close()
    