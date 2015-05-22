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
    
    

def export_coverage(cov, cov_summary_ofn):
    with open(cov_summary_ofn, "w") as OUT:
        OUT.write("id\tlength\tcoverage\n")
        for sid in cov.keys():
            OUT.write(sid+"\t"+str(cov[sid][0])+"\t"+str(cov[sid][3])+"\n")
    
    

"""
Calculate the coverage of all sequences in the coverage file
"""
def calculate_mean_coverage(coverage_fn, cov_summary_ofn=None, ref_fa_fn=None):
    from Bio import SeqIO
    
    print("Calculate coverage for " + coverage_fn)
    
    with open(coverage_fn) as IN:
        covs = IN.read().splitlines()
    
    ref_seqs_len = {}
    # Import reference sequences, and determine the percent of the ref genomes covered by reads
    if ref_fa_fn is not None:
        seqs = SeqIO.index(ref_fa_fn, "fasta")
        ref_seqs_len = {sid: len(seqs[sid].seq) for sid in seqs.keys()}
    
    
    nr_sids = list(set([cov.split("\t")[0] for cov in covs]))
    
    # len_with_coverage, total_coverage, max_coverage, mean_coverage, percent_genome_coverage
    cov_lst = {sid:[0, 0, 0, 0, 0] for sid in nr_sids}
        
    for cov in covs:     
        try:
            [sid, pos, c] = cov.split("\t")
            c = int(c)
            #if sid not in cov_lst.keys():
            #    # len_with_coverage, total_coverage, max_coverage, mean_coverage, percent_genome_coverage
            #    cov_lst[sid] = [0, 0, 0, 0, 0]

            # len_with_coverage
            cov_lst[sid][0] += 1
            # total_coverage 
            cov_lst[sid][1] += c
            # max_coverage 
            if c > cov_lst[sid][2]:
                cov_lst[sid][2] = c
        
        except ValueError:
            print(c)
    
    # Percent_genome_coverage
    if ref_fa_fn is not None:
        for sid in cov_lst.keys():
            if sid in ref_seqs_len.keys():
                cov_lst[sid][4] = cov_lst[sid][0] / ref_seqs_len[sid]
    
    for sid in cov_lst.keys():
        # Mean average
        cov_lst[sid][3] = cov_lst[sid][1] / cov_lst[sid][0]
    
    print("Number of sequences processed: " + str(len(cov_lst))) 
      
    if cov_summary_ofn is None:
        cov_summary_ofn = coverage_fn + ".summary"
    
    export_coverage(cov_lst, cov_summary_ofn)      
    
    return cov_lst


 
"""

""" 
def convert_vegan_biom_to_otu_table(biom_fn, otu_table_ofn=None, discard_eukaryota=True):
#def convert_vegan_biom_to_otu_table(biom_fn, otu_table_ofn=None, tax_map_ofn=None, discard_eukaryota=True):
    import json
    from itertools import cycle
    
    # k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Ruminococcaceae;g__;s__
    tax_ranks = ["k__", "p__", "c__", "o__","f__", "g__", "s__"]
    
    if otu_table_ofn is None:
        otu_table_ofn = biom_fn + ".otu"
    
    #if tax_map_ofn is None:
    #   tax_map_ofn = biom_fn + ".map"
    
    IN = open(biom_fn)
    biom = json.load(IN)
    
    # Determine the shape of OTU table
    nrow = biom['shape'][0]
    ncol = biom['shape'][1] 
    
    print("BIOM dimension: " + str(nrow) + ", " + str(ncol))
    
    row_ids = [r['id'] for r in biom['rows']]
    col_ids = [c['id'] for c in biom['columns']]
    OTU_tbl = [[0 for j in range(len(col_ids))] for i in range(len(row_ids))]
    
    #OUT = open(tax_map_ofn, "w")
    
    id2tax = {}
    
    # Process the tax map
    for i in range(nrow):
        id = biom['rows'][i]['id']
        tax = "None"
        if biom['rows'][i]['metadata'] is not None:
            if len(biom['rows'][i]['metadata']['Taxonomy']) > 2:
                tax = biom['rows'][i]['metadata']['Taxonomy'][2:]
                tax = tax + ["" for i in range(len(tax_ranks) - len(tax))]
                
                if "viruses" in tax[0]:
                    continue
                    
                if discard_eukaryota:
                    if tax[0] == "Eukaryota":
                        continue
                     
                tax = ";".join(z[0] + z[1] for z in zip(tax_ranks, tax))
                #tax = ";".join(biom['rows'][i]['metadata']['Taxonomy']).replace(" ", "_")
                id2tax[id] = clear_tax_label(tax)
                #OUT.write(id + "\t" + tax + "\n")
    #OUT.close()
    
    # Process OTU table
    OUT = open(otu_table_ofn, "w")
    mtx = biom['data']
    # Print column headers
    #OUT.write("\t".join(["Tax"]+col_ids) + "\n")
    OUT.write("# QIIME v1.2.1-dev OTU table\n")
    OUT.write("\t".join(["#OTU ID"]+col_ids) + "\tConsensus Lineage\n")
    for i in range(nrow):
        row_id = row_ids[i] 
        vals = "\t".join(map(str, mtx[i]))
        if row_id in id2tax.keys():
            OUT.write(row_id + "\t" + vals + "\t" + id2tax[row_id] + "\n")
    OUT.close()
    
    return otu_table_ofn



def clear_tax_label(tax_label):
    return tax_label.replace(" <phylum>", "").replace(" ", "_")
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





"""
# R script for doing CBS
library(lattice)
library(mixtools)
library(PSCBS)
library(DNAcopy)



#cov_fn = "SWH-Seed_Y0.scaffold_208.coverage"
cov_fn = "SWH-Cell55_Y2.sorted.bam.coverage"
cov <- read.table(cov_fn, sep="\t", header=F, col.names=c("chromosome", "x", "y"), stringsAsFactors=F)


scaffold_lens <- table(cov$chromosome)

len_cutoff <- 3000
filtered_scaffold_lens <- scaffold_lens[which(scaffold_lens > len_cutoff)]

i <- 2
selected_cov <- cov[which(cov$chromosome == names(filtered_scaffold_lens)[i]), ]




hist(selected_cov$y, breaks=100,prob=T); lines(density(selected_cov$y, adjust=2), lwd=1, col="red")

xyplot(selected_cov$y ~ selected_cov$x, pch=19, cex=0.2, xlim=c(min(selected_cov$x),max(selected_cov$x)) ,col='steelblue')



# http://www.r-bloggers.com/fitting-mixture-distributions-with-the-r-package-mixtools/
mixmdl <- normalmixEM(selected_cov$y)

plot(mixmdl,which=2)

cutoff = 0.95
if(length(mixmdl$lambda) > 1) {
    if(max(mixmdl$lambda) < cutoff) {
        mu <- mixmdl$mu[which.max(mixmdl$lambda)] 
        sigma <- mixmdl$sigma[which.max(mixmdl$lambda)] 
        normalized_cov <- cov
        normalized_cov$y <- (normalized_cov$y - mu) / (8 * sigma)

    } else {
        normalized_cov <- cov - mean(cov$y)
    }
} else {
    normalized_cov <- cov - mean(cov$y)
}

        
hist(normalized_cov$y, breaks=100,prob=T);
lines(density(normalized_cov$y, adjust=2), lwd=1, col="red")


segments <- segment(CNA(normalized_cov$y, normalized_cov$chromosome, normalized_cov$x))

segments <- segment(CNA(normalized_cov$y, normalized_cov$chromosome, normalized_cov$x), alpha=0.00001, undo.split="prune")

nrow(segments$segRows)



"""