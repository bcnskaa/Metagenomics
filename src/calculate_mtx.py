from __future__ import division
import os
import glob 
import numpy
from Bio import SeqIO


def main():
    # Calculate bin sequence len
    group_id = "combined"
    fasta_fn_ext = "fa"
    
    bin_fns = glob.glob(group_id + "_seqs/*." + fasta_fn_ext)
    
    bin_lens = {}
    for bin_fn in bin_fns:
        bin_id = (os.path.basename(bin_fn)).replace("." + fasta_fn_ext, "")
        print("Getting length from " + bin_id)
        seqs = SeqIO.index(bin_fn, "fasta")
        bin_lens[bin_id] = sum([len(str(seqs[sid].seq)) for sid in seqs])
    
    
    mtx = {}
    bla_fns = glob.glob("*.bla")
    for bla_fn in bla_fns:
    
        id = (os.path.basename(bla_fn)).replace(".bla", "")
        
        qid = id.split(".")[0]
        sid = id.split(".")[1]
        #qid = ".".join(id.split(".")[0:2])
        #sid = ".".join(id.split(".")[2:4])
        
        print("Processing " + id + "(sid=" + sid + " qid=" + qid + ")")
        with open(bla_fn) as IN:
            bla_res = IN.read().splitlines()
        if bla_res is None:
            mtx[id] = 0.0
        elif len(bla_res) > 0:
            #mtx[id] = numpy.mean([(int(b.split("\t")[3]) - int(int(b.split("\t")[5]) + int(b.split("\t")[4]))) / int(b.split("\t")[3]) for b in bla_res])
            mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / ((bin_lens[qid] + bin_lens[sid]) / 2)
            #mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / bin_lens[qid]
        else:
            mtx[id] = 0.0
        
    with open("bla.mtx", "w") as OUT:
        for k in mtx.keys():
            OUT.write(k + "\t" + k.replace(".", "\t") + "\t" + str(mtx[k]) + "\n")


"""
    bla_fn = "SWH-Seed_Y0.GZ-Seed_Y0.bla.b1000.0p80.0.bla"
"""
from Bio import SeqIO

def calculate_mtx_with_coverage(bla_fn, relative_abund=True, path_to_abund_dir="/home/siukinng/MG/scaffolds_5000", path_to_seq_dir="/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined", abund_dir_suffix="_5000"):
    query_metagenome_name = bla_fn.split(".")[0]
    subject_metagenome_name = bla_fn.split(".")[1]
    
    query_metagenome_abund_fn = path_to_abund_dir + "/" + query_metagenome_name + abund_dir_suffix + "/" + query_metagenome_name + ".abund"
    subject_metagenome_abund_fn = path_to_abund_dir + "/" + subject_metagenome_name + abund_dir_suffix + "/" + subject_metagenome_name + ".abund"
 
 
    score = {query_metagenome_name:0.0, subject_metagenome_name:0.0}

    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
    bla = [b.split("\t") for b in bla]
    
    if len(bla) == 0:
        return score
    
    
    with open(query_metagenome_abund_fn) as IN:
        query_metagenome_abund = IN.read().splitlines()
    if len(query_metagenome_abund) > 0:
        query_metagenome_abund = {a.split("\t")[0]:float(a.split("\t")[1]) for a in query_metagenome_abund}
    else:
        print("query metagenome " + query_metagenome_name + " does not have abundance information.")
        return score
    
    
    with open(subject_metagenome_abund_fn) as IN:
        subject_metagenome_abund = IN.read().splitlines()
    if len(subject_metagenome_abund) > 0:
        subject_metagenome_abund = {a.split("\t")[0]:float(a.split("\t")[1]) for a in subject_metagenome_abund}
    else:
        print("subject metagenome " + subject_metagenome_name + " does not have abundance information.")
        return score


    if relative_abund:
        query_metagenome_abund_total = sum([float(query_metagenome_abund[q]) for q in query_metagenome_abund.keys()])
        query_metagenome_abund = {q : float(query_metagenome_abund[q]) / query_metagenome_abund_total for q in query_metagenome_abund.keys()}
        subject_metagenome_abund_total = sum([float(subject_metagenome_abund[q]) for q in subject_metagenome_abund.keys()])
        subject_metagenome_abund = {q : float(subject_metagenome_abund[q]) / subject_metagenome_abund_total for q in subject_metagenome_abund.keys()}


    bin_lens = {}
    fa_fn = path_to_seq_dir + "/" + query_metagenome_name + ".fa"
    print("Getting length from " + fa_fn + " for " + query_metagenome_name)
    seqs = SeqIO.index(fa_fn, "fasta")
    #bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) * query_metagenome_abund[sid] for sid in seqs if sid in query_metagenome_abund.keys()]) 
    
    q_discard_n = sum([1 for sid in seqs if sid not in query_metagenome_abund.keys()])
    print("q_discard_n=" + str(q_discard_n))
    if relative_abund:
        bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) * query_metagenome_abund[sid] for sid in seqs if sid in query_metagenome_abund.keys()]) 
    else:
        bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if sid in query_metagenome_abund.keys()]) 
        #bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if len(str(seqs[sid].seq)) > 3000])       
  
    fa_fn = path_to_seq_dir + "/" + subject_metagenome_name + ".fa"
    print("Getting length from " + fa_fn + " for " + subject_metagenome_name)
    seqs = SeqIO.index(fa_fn, "fasta")
    #bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) * subject_metagenome_abund[sid] for sid in seqs if sid in subject_metagenome_abund.keys()])    
    s_discard_n = sum([1 for sid in seqs if sid not in subject_metagenome_abund.keys()])
    print("s_discard_n=" + str(s_discard_n))
    if relative_abund:
        bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) * subject_metagenome_abund[sid] for sid in seqs if sid in subject_metagenome_abund.keys()])    
    else:
        bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if sid in subject_metagenome_abund.keys()])    
        #bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if len(str(seqs[sid].seq)) > 3000])  
   
    s_bla_sum = 0.0
    q_bla_sum = 0.0
    for b in bla:
        q_id = b[0]
        s_id = b[1]
        if s_id in subject_metagenome_abund.keys() and q_id in query_metagenome_abund.keys():
        
            s = (int(b[3]) - ( int(b[5]) + int(b[4]) ))
            if relative_abund:
                s_weight = subject_metagenome_abund[s_id]
                q_weight = query_metagenome_abund[q_id]
                s_bla_sum = s_bla_sum + (s_weight * s)
                q_bla_sum = q_bla_sum + (q_weight * s)
            else:
                s_bla_sum = s_bla_sum + s
                q_bla_sum = q_bla_sum + s
                
        #s_bla_sum = s_bla_sum + int(b.split("\t")[3]) - ( int(b.split("\t")[5]) + int(b.split("\t")[4]) )
        #q_bla_sum = q_bla_sum + int(b.split("\t")[3]) - ( int(b.split("\t")[5]) + int(b.split("\t")[4]) )   
    score[subject_metagenome_name] = s_bla_sum / bin_lens[subject_metagenome_name]
    score[query_metagenome_name] = q_bla_sum / bin_lens[query_metagenome_name]
        
    #mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / ((bin_lens[qid] + bin_lens[sid]) / 2)
            
    return score


    
"""   
    /home/siukinng/MG/scaffolds/CrossValidation/db/combined
    for f in *.fa;do ~/tools/blast/bin/makeblastdb -dbtype nucl -in $f;done
    
    library(ggplot2)
    pcoa <- read.table("PCoA_ident70_minlen30_eval10-8.txt", sep="\t", header=T, stringsAsFactors=F)
    pcoa$sample_id <- sapply(1 : nrow(pcoa), function(i) { strsplit(pcoa[i,1],"_")[[1]][1]} )
    pcoa$substrate <- sapply(1 : nrow(pcoa), function(i) { strsplit(pcoa[i,1],"_")[[1]][2]} )
    ggplot(pcoa, aes(x=PCO1, y=PCO2, color=substrate, shape=sample_id)) + geom_point()
    
"""


"""
import calculate_mtx
import os

sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2"]

weighted_score_mtx = [[0 for x in range(len(sample_ids))] for x in range(len(sample_ids))] 
#fn_suffix = ".bla.b1000.0p80.0.bla"
fn_suffix = ".bla.b1000.0.bla"
for i, sample_id_s in enumerate(sample_ids):
    for j, sample_id_q in enumerate(sample_ids):
        print("Processing " + sample_id_s + " and " + sample_id_q)
        bla_fn = sample_id_q + "." + sample_id_s + fn_suffix
        if os.path.isfile(bla_fn):
            #scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=True)
            scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=True, path_to_abund_dir="/home/siukinng/MG/scaffolds_2000", path_to_seq_dir="/home/siukinng/MG/scaffolds_2000/CrossValidation/db/combined", abund_dir_suffix="_2000");
            weighted_score_mtx[i][j] = scores[sample_id_s]
            weighted_score_mtx[j][i] = scores[sample_id_q]
            print(sample_id_s + " vs " + sample_id_q + " = " + str(weighted_score_mtx[i][j]) + ", " + str(weighted_score_mtx[j][i]))


with open("weighted_score_mtx", "w") as OUT:
    OUT.write("Label\t" + "\t".join(sample_ids) + "\n")
    for i, sample_id in enumerate(sample_ids):
        OUT.write(sample_id + "\t" + "\t".join([str(s) for s in weighted_score_mtx[i]]) + "\n")


score_mtx = [[0 for x in range(len(sample_ids))] for x in range(len(sample_ids))] 


for i, sample_id_s in enumerate(sample_ids):
    for j, sample_id_q in enumerate(sample_ids):
        print("Processing " + sample_id_s + " and " + sample_id_q)
        bla_fn = sample_id_q + "." + sample_id_s + fn_suffix
        if os.path.isfile(bla_fn):
            #scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=False)
            scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=False, path_to_abund_dir="/home/siukinng/MG/scaffolds_2000", path_to_seq_dir="/home/siukinng/MG/scaffolds_2000/CrossValidation/db/combined", abund_dir_suffix="_2000");

            score_mtx[i][j] = scores[sample_id_s]
            score_mtx[j][i] = scores[sample_id_q]
            print(sample_id_s + " vs " + sample_id_q + " = " + str(score_mtx[i][j]) + ", " + str(score_mtx[j][i]))

with open("score_mtx", "w") as OUT:
    OUT.write("Label\t" + "\t".join(sample_ids) + "\n")
    for i, sample_id in enumerate(sample_ids):
        OUT.write(sample_id + "\t" + "\t".join([str(s) for s in score_mtx[i]]) + "\n")


"""


"""
import glob

abund_fns = glob.glob("*_5000/*.abund")
for abund_fn in abund_fns:
    print("Processing " + abund_fn + "\n")
    with open(abund_fn) as IN:
        abund = IN.read().splitlines()
    total_abund = sum([float(a.split("\t")[1]) for a in abund])
    abund = {a.split("\t")[0]:float(a.split("\t")[1])/total_abund for a in abund}

    with open(abund_fn.replace(".abund",".relative_abund"), "w") as OUT:
        for ak in abund.keys():
            OUT.write(ak + "\t" + str(abund[ak]) + "\n")

"""


"""
from __future__ import division
import sys
sys.modules[__name__].__dict__.clear()


bla_fn = "SWH-Xyl_Y1.SWH-Cell_Y2.bla.b1000.0.bla"
sample_id = "SWH-Cell_Y2"

with open(bla_fn) as IN:
    bla = IN.read().splitlines()
bla = [b.split("\t") for b in bla]


abund_fn = "/disk/rdisk08/siukinng/MG/scaffolds_5000/" + sample_id + "_5000/" + sample_id + ".relative_abund"
with open(abund_fn) as IN:
    abund = IN.read().splitlines()
abund = {a.split("\t")[0]:float(a.split("\t")[1]) for a in abund}   

len_abund = 0.0
for b in bla:
    id = b[1]
    l = int(b[3]) - (int(4) + int(b[5]))
    print("id=" + id + ", len=" + str(l) + ", abund=" + str(abund[id]) + ", a=" + str((abund[id] * l))) 
    len_abund = len_abund + (abund[id] * l)
    

from Bio import SeqIO
fa_fn = "/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined/" + sample_id + ".fa"
seqs = SeqIO.index(fa_fn, "fasta")
total_len_abund = 0.0
total_len = 0
for s_id in seqs:
    l = len(str(seqs[s_id].seq))
    a = abund[s_id]
    print("s_id=" + s_id + ", len=" + str(l) + ", abund=" + str(a) + ", a=" + str((a * l))) 
    total_len = total_len + l
    total_len_abund = total_len_abund + (a * l)
print(len_abund / total_len_abund)
    
"""
