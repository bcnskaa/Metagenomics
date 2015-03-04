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

cache = {}
def calculate_mtx_with_coverage(bla_fn, metagenome_fa_fn_ext=".fa", subject_metagenome_fa_fn=None, query_metagenome_fa_fn=None, subject_metagenome_abund_fn=None, query_metagenome_abund_fn=None, raw_ofn=None, relative_abund=True, path_to_abund_dir="/home/siukinng/MG/scaffolds_5000", path_to_seq_dir="/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined", abund_dir_suffix="_5000", id_separator='.'):
    import os
    
    id = (os.path.basename(bla_fn)[::-1].split(".",1)[1])[::-1]
    
    query_metagenome_name = id.split(id_separator)[0]
    subject_metagenome_name = id.split(id_separator)[1]
    #query_metagenome_name = bla_fn.split(id_separator)[0]
    #subject_metagenome_name = bla_fn.split(id_separator)[1]
    
#    if query_metagenome_name not in cache.keys():
        
    
    if query_metagenome_abund_fn is None:
        query_metagenome_abund_fn = path_to_abund_dir + "/" + query_metagenome_name + abund_dir_suffix + "/" + query_metagenome_name + ".abund"
    if subject_metagenome_abund_fn is None:
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
    if query_metagenome_fa_fn is None:
        fa_fn = path_to_seq_dir + "/" + query_metagenome_name + metagenome_fa_fn_ext
    else:
        fa_fn = query_metagenome_fa_fn
        
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

    if subject_metagenome_fa_fn is None:
        fa_fn = path_to_seq_dir + "/" + subject_metagenome_name + metagenome_fa_fn_ext
    else:
        fa_fn = subject_metagenome_fa_fn
        
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
   
    discard_bla_n = 0
    raw_data = []
    s_bla_sum = 0.0
    q_bla_sum = 0.0
    for b in bla:
        q_id = b[0]
        s_id = b[1]
        if s_id in subject_metagenome_abund.keys() and q_id in query_metagenome_abund.keys():
            s = (int(b[3]) - ( int(b[5]) + int(b[4]) ))
            raw = [q_id, s_id]
            if relative_abund:
                s_weight = subject_metagenome_abund[s_id]
                q_weight = query_metagenome_abund[q_id]
                s_bla_sum = s_bla_sum + (s_weight * s)
                q_bla_sum = q_bla_sum + (q_weight * s)
                raw.append(str((q_weight * s))) 
                raw.append(str((s_weight * s)))
            else:
                s_bla_sum = s_bla_sum + s
                q_bla_sum = q_bla_sum + s
                raw.append(str(s))
                raw.append(str(s))
            raw_data.append(raw)
        else:
            discard_bla_n = discard_bla_n + 1 
            
        #s_bla_sum = s_bla_sum + int(b.split("\t")[3]) - ( int(b.split("\t")[5]) + int(b.split("\t")[4]) )
        #q_bla_sum = q_bla_sum + int(b.split("\t")[3]) - ( int(b.split("\t")[5]) + int(b.split("\t")[4]) )   
    print("Dicard_bla_n=" + str(discard_bla_n))
    
    if raw_ofn is not None:
        with open(raw_ofn, "w") as OUT:
            for raw in raw_data:
                OUT.write("\t".join(raw) + "\n")

    score[subject_metagenome_name] = s_bla_sum / bin_lens[subject_metagenome_name]
    score[query_metagenome_name] = q_bla_sum / bin_lens[query_metagenome_name]
        
    #mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / ((bin_lens[qid] + bin_lens[sid]) / 2)
            
    return score
# def calculate_mtx_with_coverage(bla_fn, metagenome_fa_fn_ext=".fa", subject_metagenome_fa_fn=None, query_metagenome_fa_fn=None, subject_metagenome_abund_fn=None, query_metagenome_abund_fn=None, raw_ofn=None, relative_abund=True, path_to_abund_dir="/home/siukinng/MG/scaffolds_5000", path_to_seq_dir="/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined", abund_dir_suffix="_5000", id_separator='.'):
#     import os
#     
#     id = (os.path.basename(bla_fn)[::-1].split(".",1)[1])[::-1]
#     
#     query_metagenome_name = id.split(id_separator)[0]
#     subject_metagenome_name = id.split(id_separator)[1]
#     #query_metagenome_name = bla_fn.split(id_separator)[0]
#     #subject_metagenome_name = bla_fn.split(id_separator)[1]
#     
#     if query_metagenome_abund_fn is None:
#         query_metagenome_abund_fn = path_to_abund_dir + "/" + query_metagenome_name + abund_dir_suffix + "/" + query_metagenome_name + ".abund"
#     if subject_metagenome_abund_fn is None:
#         subject_metagenome_abund_fn = path_to_abund_dir + "/" + subject_metagenome_name + abund_dir_suffix + "/" + subject_metagenome_name + ".abund"
#  
#  
#     score = {query_metagenome_name:0.0, subject_metagenome_name:0.0}
# 
#     with open(bla_fn) as IN:
#         bla = IN.read().splitlines()
#     bla = [b.split("\t") for b in bla]
#     
#     if len(bla) == 0:
#         return score
#     
#     
#     with open(query_metagenome_abund_fn) as IN:
#         query_metagenome_abund = IN.read().splitlines()
#     if len(query_metagenome_abund) > 0:
#         query_metagenome_abund = {a.split("\t")[0]:float(a.split("\t")[1]) for a in query_metagenome_abund}
#     else:
#         print("query metagenome " + query_metagenome_name + " does not have abundance information.")
#         return score
#     
#     
#     with open(subject_metagenome_abund_fn) as IN:
#         subject_metagenome_abund = IN.read().splitlines()
#     if len(subject_metagenome_abund) > 0:
#         subject_metagenome_abund = {a.split("\t")[0]:float(a.split("\t")[1]) for a in subject_metagenome_abund}
#     else:
#         print("subject metagenome " + subject_metagenome_name + " does not have abundance information.")
#         return score
# 
# 
#     if relative_abund:
#         query_metagenome_abund_total = sum([float(query_metagenome_abund[q]) for q in query_metagenome_abund.keys()])
#         query_metagenome_abund = {q : float(query_metagenome_abund[q]) / query_metagenome_abund_total for q in query_metagenome_abund.keys()}
#         subject_metagenome_abund_total = sum([float(subject_metagenome_abund[q]) for q in subject_metagenome_abund.keys()])
#         subject_metagenome_abund = {q : float(subject_metagenome_abund[q]) / subject_metagenome_abund_total for q in subject_metagenome_abund.keys()}
# 
# 
#     bin_lens = {}
#     if query_metagenome_fa_fn is None:
#         fa_fn = path_to_seq_dir + "/" + query_metagenome_name + metagenome_fa_fn_ext
#     else:
#         fa_fn = query_metagenome_fa_fn
#         
#     print("Getting length from " + fa_fn + " for " + query_metagenome_name)
#     seqs = SeqIO.index(fa_fn, "fasta")
#     #bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) * query_metagenome_abund[sid] for sid in seqs if sid in query_metagenome_abund.keys()]) 
#     
#     q_discard_n = sum([1 for sid in seqs if sid not in query_metagenome_abund.keys()])
#     print("q_discard_n=" + str(q_discard_n))
#     if relative_abund:
#         bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) * query_metagenome_abund[sid] for sid in seqs if sid in query_metagenome_abund.keys()]) 
#     else:
#         bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if sid in query_metagenome_abund.keys()]) 
#         #bin_lens[query_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if len(str(seqs[sid].seq)) > 3000])       
# 
#     if subject_metagenome_fa_fn is None:
#         fa_fn = path_to_seq_dir + "/" + subject_metagenome_name + metagenome_fa_fn_ext
#     else:
#         fa_fn = subject_metagenome_fa_fn
#         
#     print("Getting length from " + fa_fn + " for " + subject_metagenome_name)
#     seqs = SeqIO.index(fa_fn, "fasta")
#     #bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) * subject_metagenome_abund[sid] for sid in seqs if sid in subject_metagenome_abund.keys()])    
#     s_discard_n = sum([1 for sid in seqs if sid not in subject_metagenome_abund.keys()])
#     print("s_discard_n=" + str(s_discard_n))
#     if relative_abund:
#         bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) * subject_metagenome_abund[sid] for sid in seqs if sid in subject_metagenome_abund.keys()])    
#     else:
#         bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if sid in subject_metagenome_abund.keys()])    
#         #bin_lens[subject_metagenome_name] = sum([len(str(seqs[sid].seq)) for sid in seqs if len(str(seqs[sid].seq)) > 3000])  
#    
#     discard_bla_n = 0
#     raw_data = []
#     s_bla_sum = 0.0
#     q_bla_sum = 0.0
#     for b in bla:
#         q_id = b[0]
#         s_id = b[1]
#         if s_id in subject_metagenome_abund.keys() and q_id in query_metagenome_abund.keys():
#             s = (int(b[3]) - ( int(b[5]) + int(b[4]) ))
#             raw = [q_id, s_id]
#             if relative_abund:
#                 s_weight = subject_metagenome_abund[s_id]
#                 q_weight = query_metagenome_abund[q_id]
#                 s_bla_sum = s_bla_sum + (s_weight * s)
#                 q_bla_sum = q_bla_sum + (q_weight * s)
#                 raw.append(str((q_weight * s))) 
#                 raw.append(str((s_weight * s)))
#             else:
#                 s_bla_sum = s_bla_sum + s
#                 q_bla_sum = q_bla_sum + s
#                 raw.append(str(s))
#                 raw.append(str(s))
#             raw_data.append(raw)
#         else:
#             discard_bla_n = discard_bla_n + 1 
#             
#         #s_bla_sum = s_bla_sum + int(b.split("\t")[3]) - ( int(b.split("\t")[5]) + int(b.split("\t")[4]) )
#         #q_bla_sum = q_bla_sum + int(b.split("\t")[3]) - ( int(b.split("\t")[5]) + int(b.split("\t")[4]) )   
#     print("Dicard_bla_n=" + str(discard_bla_n))
#     
#     if raw_ofn is not None:
#         with open(raw_ofn, "w") as OUT:
#             for raw in raw_data:
#                 OUT.write("\t".join(raw) + "\n")
# 
#     score[subject_metagenome_name] = s_bla_sum / bin_lens[subject_metagenome_name]
#     score[query_metagenome_name] = q_bla_sum / bin_lens[query_metagenome_name]
#         
#     #mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / ((bin_lens[qid] + bin_lens[sid]) / 2)
#             
#     return score


"""

import calculate_mtx
calculate_mtx.blast_bins("/home/siukinng/MG/scaffolds_5000/GZ-Cell_Y2_5000", "/home/siukinng/MG/scaffolds_5000/GZ-Cell_Y1_5000")

import calculate_mtx
calculate_mtx.process_all_samples()


import calculate_mtx
sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0"]
calculate_mtx.process_all_samples_to_reference(sample_ids, "/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta")

import calculate_mtx
sample_ids = ["GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2"]
calculate_mtx.process_all_samples_to_reference(sample_ids, "/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta")

import calculate_mtx
sample_ids = ["SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1"]
calculate_mtx.process_all_samples_to_reference(sample_ids, "/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta")

import calculate_mtx
sample_ids = ["SWH-Cell_Y2", "SWH-Cell55_Y2"]
calculate_mtx.process_all_samples_to_reference(sample_ids, "/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta")



import calculate_mtx
sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0"]
calculate_mtx.process_all_samples_to_reference(sample_ids, db_fn="/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta", fa_fn_ext=".fa")

import calculate_mtx
sample_ids = ["GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2"]
calculate_mtx.process_all_samples_to_reference(sample_ids, db_fn="/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta", fa_fn_ext=".fa")

import calculate_mtx
sample_ids = ["SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2"]
calculate_mtx.process_all_samples_to_reference(sample_ids, db_fn="/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta", fa_fn_ext=".fa")


import calculate_mtx
sample_ids = ["SWH-Xyl_Y1", "SWH-Cell55_Y2"]
calculate_mtx.process_all_samples_to_reference(sample_ids, db_fn="/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes.fasta", fa_fn_ext=".fa")


import calculate_mtx
sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
calculate_mtx.process_all_samples_to_reference(sample_ids, db_fn="/home/siukinng/db/BioProject_Prokaryotes/all_prokaryotes+BacterialDB.fasta", fa_fn_ext=".fasta")



"""
def process_all_samples_to_reference(sample_ids, db_fn="/home/siukinng/db/BacteriaDB/all_fna.fna", binning_len=5000, fa_fn_ext=".fasta"):
    import glob 
    #sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    map_to_ref_dir = "./blast_to_ref"
    if not os.path.isdir(map_to_ref_dir):
        os.makedirs(map_to_ref_dir)
    
    for sample_id in sample_ids:
        fa_fns = glob.glob("/home/siukinng/MG/scaffolds_" + str(binning_len) + "/"+sample_id+"_" + str(binning_len) + "/*"+fa_fn_ext)
        print("Number of fasta files to be processed for " + sample_id + ": " + str(len(fa_fns)))
        for fa_fn in fa_fns:
            print("Processing " + fa_fn)
            bin_id = (os.path.basename(fa_fn)).replace(fa_fn_ext, "")
            #out_fn = bin_id + "+all_fna.bla"
            out_fn = bin_id + "+" + os.path.basename(db_fn) + ".bla"
            mg_pipeline.blastn(fa_fn, db_fn, outdir=map_to_ref_dir, outfn=out_fn, num_threads=16, perc_identity=75, max_target_seqs=1, evalue=1e-50)


"""
# Bin to bin
"""
def process_all_samples():
    sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    for s_sample_id in sample_ids:
        for q_sample_id in sample_ids:
            blast_bins("/home/siukinng/MG/scaffolds_5000/"+q_sample_id+"_5000", "/home/siukinng/MG/scaffolds_5000/"+s_sample_id+"_5000", out_dir=".", seq_fn_ext=".fasta")

 

import mg_pipeline
import glob
import os
def blast_bins(q_dir, s_dir, out_dir=".", seq_fn_ext=".fasta"):    
    s_fa_fns = glob.glob(s_dir + "/*" + seq_fn_ext)
    q_fa_fns = glob.glob(q_dir + "/*" + seq_fn_ext)
    
    s_sample_id = os.path.basename(s_dir)
    s_sample_outdir = out_dir + "/" + s_sample_id
    if not os.path.isdir(s_sample_outdir):
        os.makedirs(out_dir + "/" + s_sample_id)
    for s_fa_fn in s_fa_fns:
        cmd = "ln -s " + s_fa_fn + " " + s_sample_outdir + "/" + os.path.basename(s_fa_fn)
        os.system(cmd)
    
    q_sample_id = os.path.basename(q_dir)
    q_sample_outdir = out_dir + "/" + q_sample_id
    if not os.path.isdir(q_sample_outdir):
        os.makedirs(out_dir + "/" + q_sample_id)
    for q_fa_fn in q_fa_fns:
        cmd = "ln -s " + q_fa_fn + " " + q_sample_outdir + "/" + os.path.basename(q_fa_fn)
        os.system(cmd)
            
    s_fa_fns = glob.glob(s_sample_outdir + "/*" + seq_fn_ext)
    q_fa_fns = glob.glob(q_sample_outdir + "/*" + seq_fn_ext)
    for s_fa_fn in s_fa_fns:
        cmd = "~/tools/blast/bin/makeblastdb -dbtype nucl -in " + s_fa_fn
        os.system(cmd)
        
    for q_fa_fn in q_fa_fns:
        cmd = "~/tools/blast/bin/makeblastdb -dbtype nucl -in " + q_fa_fn
        os.system(cmd)
         
    for s_fa_fn in s_fa_fns:
        s_id = os.path.basename(s_fa_fn).replace(seq_fn_ext, "")
        for q_fa_fn in q_fa_fns:
            q_id = os.path.basename(q_fa_fn).replace(seq_fn_ext, "")
            out_fn = q_id + "+" + s_id + ".bla"
            mg_pipeline.blastn(q_fa_fn, s_fa_fn, outdir=out_dir, outfn=out_fn, num_threads=16, perc_identity=75, max_target_seqs=1, evalue=1e-50)
            


"""

import calculate_mtx
q_id = "SWH-Cell_Y2"
s_id = "GZ-Cell_Y2"
raw_fn = q_id + "." + s_id + ".bla.b1000.0.bla.raw"
q_scaffold2tax_fn = "/home/siukinng/MG/scaffolds_5000/" + q_id + "_5000/" + q_id + ".scaffold2tax_f"
s_scaffold2tax_fn = "/home/siukinng/MG/scaffolds_5000/" + s_id + "_5000/" + s_id + ".scaffold2tax_f"
ss = calculate_mtx.map_raw_with_taxon(raw_fn, q_scaffold2tax_fn, s_scaffold2tax_fn)



for s in :


"""
def map_raw_with_taxon(raw_fn, q_scaffold2tax_fn, s_scaffold2tax_fn, out_fn=None):
    # Import raw data
    with open(raw_fn) as IN:
        raw = IN.read().splitlines()
    raw = [r.split("\t") for r in raw]
    
    # import scaffold map for query ids
    with open(q_scaffold2tax_fn) as IN:
        q_scaffold2tax = IN.read().splitlines()
    q_scaffold2tax = {r.split("\t")[0]:r.split("\t")[1] for r in q_scaffold2tax}
    
    # import scaffold map for subject ids
    with open(s_scaffold2tax_fn) as IN:
        s_scaffold2tax = IN.read().splitlines()
    s_scaffold2tax = {r.split("\t")[0]:r.split("\t")[1] for r in s_scaffold2tax}

    # Check out_fn    
    if out_fn is None:
        out_fn = raw_fn + ".mapped"
    
    q_total = {}
    s_total = {}
    
    discard_n = 0
    # Perform mapping
    with open(out_fn, "w") as OUT:
        for r in raw:
            q_id = r[0]
            s_id = r[0]
            if q_id not in q_scaffold2tax.keys():
                discard_n = discard_n + 1
                continue
            q_tax = q_scaffold2tax[q_id]
            if s_id not in s_scaffold2tax.keys():
                discard_n = discard_n + 1
                continue
            s_tax = s_scaffold2tax[s_id]
            q_val = float(r[2])
            s_val = float(r[3])
            
            if q_tax not in q_total.keys():
                q_total[q_tax] = 0.0
            q_total[q_tax] = q_total[q_tax] + q_val
            
            if s_tax not in s_total.keys(): 
                s_total[s_tax] = 0.0
            s_total[s_tax] = s_total[s_tax] + s_val  
             
            OUT.write("\t".join([q_tax, s_tax, r[2], r[3]]) + "\n")
            
    #for q in q_total.keys():
    #    q_total[q] = q_total[q]
    q_total_len = sum([q_total[q] for q in q_total.keys()])
    q_total = {q:q_total[q]/q_total_len for q in q_total.keys()}
    
    s_total_len = sum([s_total[s] for s in s_total.keys()])
    s_total = {s:s_total[s]/s_total_len for s in s_total.keys()}
       
    return [q_total, s_total] 
    


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

sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]

weighted_score_mtx = [[0 for x in range(len(sample_ids))] for x in range(len(sample_ids))] 
#fn_suffix = ".bla.b1000.0p80.0.bla"
fn_suffix = ".bla.b1000.0.bla"
for i, sample_id_s in enumerate(sample_ids):
    for j, sample_id_q in enumerate(sample_ids):
        print("Processing " + sample_id_s + " and " + sample_id_q)
        bla_fn = sample_id_q + "." + sample_id_s + fn_suffix
        if os.path.isfile(bla_fn):
            #scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=True)
            scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, raw_ofn=bla_fn + ".raw", relative_abund=True, path_to_abund_dir="/home/siukinng/MG/scaffolds_5000", path_to_seq_dir="/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined", abund_dir_suffix="_5000");
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
            scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=False, path_to_abund_dir="/home/siukinng/MG/scaffolds_5000", path_to_seq_dir="/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined", abund_dir_suffix="_5000");

            score_mtx[i][j] = scores[sample_id_s]
            score_mtx[j][i] = scores[sample_id_q]
            print(sample_id_s + " vs " + sample_id_q + " = " + str(score_mtx[i][j]) + ", " + str(score_mtx[j][i]))

with open("score_mtx", "w") as OUT:
    OUT.write("Label\t" + "\t".join(sample_ids) + "\n")
    for i, sample_id in enumerate(sample_ids):
        OUT.write(sample_id + "\t" + "\t".join([str(s) for s in score_mtx[i]]) + "\n")


"""


def calculate_weighted_score_mtx(weighted_score_mtx_fn="weighted_score_mtx", fn_suffix = ".bla.b1000.0.bla",  id_separator="."):
    import calculate_mtx
    import os
    
    sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    
    weighted_score_mtx = [[0 for x in range(len(sample_ids))] for x in range(len(sample_ids))] 
    #fn_suffix = ".bla.b1000.0p80.0.bla"
    fn_suffix = ".bla.b1000.0.bla"
    for i, sample_id_s in enumerate(sample_ids):
        for j, sample_id_q in enumerate(sample_ids):
            print("Processing " + sample_id_s + " and " + sample_id_q)
            bla_fn = sample_id_q + id_separator + sample_id_s + fn_suffix
            if os.path.isfile(bla_fn):
                #scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=True)
                scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, raw_ofn=bla_fn + ".raw", relative_abund=True, path_to_abund_dir="/home/siukinng/MG/scaffolds_5000", path_to_seq_dir="/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined", abund_dir_suffix="_5000");
                weighted_score_mtx[i][j] = scores[sample_id_s]
                weighted_score_mtx[j][i] = scores[sample_id_q]
                print(sample_id_s + " vs " + sample_id_q + " = " + str(weighted_score_mtx[i][j]) + ", " + str(weighted_score_mtx[j][i])) 
    
    with open(weighted_score_mtx_fn, "w") as OUT:
        OUT.write("Label\t" + "\t".join(sample_ids) + "\n")
        for i, sample_id in enumerate(sample_ids):
            OUT.write(sample_id + "\t" + "\t".join([str(s) for s in weighted_score_mtx[i]]) + "\n")






   
def calculate_score_mtx(score_mtx_fn="score_mtx", fn_suffix = ".bla.b1000.0.bla", id_separator="."):
    import calculate_mtx
    import os
    
    sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    score_mtx = [[0 for x in range(len(sample_ids))] for x in range(len(sample_ids))] 

    for i, sample_id_s in enumerate(sample_ids):
        for j, sample_id_q in enumerate(sample_ids):
            print("Processing " + sample_id_s + " and " + sample_id_q)
            bla_fn = sample_id_q + id_separator + sample_id_s + fn_suffix
            if os.path.isfile(bla_fn):
                #scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=False)
                scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, relative_abund=False, path_to_abund_dir="/home/siukinng/MG/scaffolds_5000", path_to_seq_dir="/home/siukinng/MG/scaffolds_5000/CrossValidation/db/combined", abund_dir_suffix="_5000");
    
                score_mtx[i][j] = scores[sample_id_s]
                score_mtx[j][i] = scores[sample_id_q]
                print(sample_id_s + " vs " + sample_id_q + " = " + str(score_mtx[i][j]) + ", " + str(score_mtx[j][i]))
    
    if score_mtx_fn is None:
        score_mtx_fn = "score_mtx"
        
    with open(score_mtx_fn, "w") as OUT:
        OUT.write("Label\t" + "\t".join(sample_ids) + "\n")
        for i, sample_id in enumerate(sample_ids):
            OUT.write(sample_id + "\t" + "\t".join([str(s) for s in score_mtx[i]]) + "\n")  
    
    
"""
import os 
cmd = "ln -s ~/tools/scripts/calculate_mtx.py"
os.system(cmd)
cmd = "ln -s ~/tools/scripts/mg_pipeline.py"
os.system(cmd)

import calculate_mtx

calculate_mtx.calculate_score_mtx_for_all_bins(bla_dir="bins", binning_size="2000")

"""
def calculate_score_mtx_for_all_bins(binning_size="5000", score_mtx_fn="score_mtx_bins", weighted_score_mtx_fn="weighted_score_mtx_bins", bla_dir=".", fn_suffix=".bla", id_separator="+"):
    import calculate_mtx
    import os
    import glob
    
    path_to_abund_dir="/home/siukinng/MG/scaffolds_" + binning_size
    sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    bin_ids = {}
    melt_bin_ids = []
    bin_id2_sample_id = {}
    
    bla_fns = glob.glob(bla_dir + "/*.bla")
    
    # estimate the number of bin groups for each sample
    for bla_fn in bla_fns:
        bin_id = (os.path.basename(bla_fn)).split(id_separator)[0]
        sample_id = bin_id.split(".")[0]
        bin_id2_sample_id[bin_id] = sample_id
        if sample_id not in bin_ids.keys():
            bin_ids[sample_id] = []
        if bin_id not in bin_ids[sample_id]:
            bin_ids[sample_id].append(bin_id)
    
    # Construct the melt bin_ids
    for sample_id in sample_ids:
        print(sample_id + ": " + str(len(bin_ids[sample_id])))
        melt_bin_ids = melt_bin_ids + bin_ids[sample_id]
    
    melt_bin_ids = sorted(melt_bin_ids)
    
#     for s_sample_id in sample_ids:
#         for q_sample_id in sample_ids:
#             for s_bin_id in bin_ids[s_sample_id]:
#                 for q_bin_id in bin_ids[q_sample_id]:
#                     bla_fn = q_bin_id + id_separator + s_bin_id + ".bla"

    score_mtx = [[0.0 for x in range(len(melt_bin_ids))] for x in range(len(melt_bin_ids))] 
    weighted_score_mtx = [[0.0 for x in range(len(melt_bin_ids))] for x in range(len(melt_bin_ids))] 
    for i, s_bin_id in enumerate(melt_bin_ids):
        for j, q_bin_id in enumerate(melt_bin_ids):
            bla_fn = bla_dir + "/" + q_bin_id + id_separator + s_bin_id + ".bla"
            if os.path.isfile(bla_fn): 
                scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, metagenome_fa_fn_ext=".fasta", subject_metagenome_fa_fn=path_to_abund_dir+"/"+bin_id2_sample_id[s_bin_id] + "_" + binning_size + "/" + s_bin_id+".fasta", query_metagenome_fa_fn=path_to_abund_dir+"/"+bin_id2_sample_id[q_bin_id] + "_" + binning_size + "/" + q_bin_id+".fasta", subject_metagenome_abund_fn=path_to_abund_dir+"/"+bin_id2_sample_id[s_bin_id] + "_"+binning_size+"/" + bin_id2_sample_id[s_bin_id]+".abund", query_metagenome_abund_fn=path_to_abund_dir+"/"+bin_id2_sample_id[q_bin_id] + "_"+binning_size+"/" + bin_id2_sample_id[q_bin_id]+".abund", relative_abund=False, path_to_abund_dir="/home/siukinng/MG/scaffolds_"+binning_size, path_to_seq_dir="/home/siukinng/MG/scaffolds_"+binning_size, abund_dir_suffix="_"+binning_size, id_separator=id_separator);
                score_mtx[i][j] = scores[s_bin_id]
                score_mtx[j][i] = scores[q_bin_id]
                 
                weighted_scores = calculate_mtx.calculate_mtx_with_coverage(bla_fn, metagenome_fa_fn_ext=".fasta", subject_metagenome_fa_fn=path_to_abund_dir+"/"+bin_id2_sample_id[s_bin_id] + "_" + binning_size +"/" + s_bin_id+".fasta", query_metagenome_fa_fn=path_to_abund_dir+"/"+bin_id2_sample_id[q_bin_id] + "_" + binning_size +"/" + q_bin_id+".fasta", subject_metagenome_abund_fn=path_to_abund_dir+"/"+bin_id2_sample_id[s_bin_id] + "_"+binning_size+"/" + bin_id2_sample_id[s_bin_id]+".abund", query_metagenome_abund_fn=path_to_abund_dir+"/"+bin_id2_sample_id[q_bin_id] + "_"+binning_size+"/" + bin_id2_sample_id[q_bin_id]+".abund", relative_abund=True, path_to_abund_dir="/home/siukinng/MG/scaffolds_"+binning_size, path_to_seq_dir="/home/siukinng/MG/scaffolds_"+binning_size, abund_dir_suffix="_"+binning_size, id_separator=id_separator);
                weighted_score_mtx[i][j] = weighted_scores[s_bin_id]
                weighted_score_mtx[j][i] = weighted_scores[q_bin_id]                 
                print(s_bin_id + " vs " + q_bin_id + " = " + str(score_mtx[i][j]) + ", " + str(score_mtx[j][i]) + " weighted: " + str(weighted_score_mtx[i][j]) + ", " + str(weighted_score_mtx[j][i]))
                print(s_bin_id + " vs " + q_bin_id + " weighted: " + str(weighted_score_mtx[i][j]) + ", " + str(weighted_score_mtx[j][i]))
    
    if score_mtx_fn is None:
         score_mtx_fn = "score_mtx_bins"
    with open(score_mtx_fn, "w") as OUT:
         OUT.write("Label\t" + "\t".join(melt_bin_ids) + "\n")
         for i, bin_id in enumerate(melt_bin_ids):
             OUT.write(bin_id + "\t" + "\t".join([str(s) for s in score_mtx[i]]) + "\n")  
     
    if weighted_score_mtx_fn is None:
        weighted_score_mtx_fn = "weighted_score_mtx_bins"
        
    with open(weighted_score_mtx_fn, "w") as OUT:
        OUT.write("Label\t" + "\t".join(melt_bin_ids) + "\n")
        for i, bin_id in enumerate(melt_bin_ids):
            OUT.write(bin_id + "\t" + "\t".join([str(s) for s in weighted_score_mtx[i]]) + "\n")  
 
  

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



"""
import os
import glob
from Bio import SeqIO

fa_fns = glob.glob("*.fasta")
ids = {(f.split("_")[2].split("+")[0][0])+(f.split("_")[2].split("+")[1][0:2]):f for f in fa_fns}

out_fn = "testtest"
list_map = {}
combined_seqs = {}

for id in ids:
    seqs = SeqIO.index(ids[id], "fasta")

    for seq_id in seqs.keys():
        new_seq_id = id + "_" + seq_id.split("|")[1]
        list_map[new_seq_id] = seqs[seq_id].description
        combined_seqs[new_seq_id] = str(seqs[seq_id].seq)


with open(out_fn, "w") as OUT:
    for seq_id in combined_seqs.keys():
        OUT.write(">" + seq_id + "\n" + combined_seqs[seq_id] + "\n\n")
        
    
with open(out_fn + ".map", "w") as OUT: 
    for seq_id in list_map.keys():
        OUT.write(seq_id + "\t" + list_map[seq_id] + "\n")
    


import os
import glob
from Bio import SeqIO

sample_id = "SC2"
#old_sample_id = sample_id
old_sample_id = "SWH-Cell_Y2"
faa_fn = sample_id + "_meta.pfam.faa.raw" 
faa_ofn = sample_id + "_meta.pfam.faa" 

#seqs = SeqIO.index(faa_fn, "fasta")
seqs = remove_duplicate(faa_fn)
with open(faa_ofn, "w") as OUT:
    for seq_id in seqs.keys():
        #OUT.write(">" + seq_id.replace("scaffold", old_sample_id) + "\n" + str(seqs[seq_id].seq).replace("*", "") + "\n")
        OUT.write(">" + seq_id.replace(old_sample_id, sample_id) + "\n" + (seqs[seq_id]).replace("*", "") + "\n")


"""

def remove_duplicate(faa_fn):
    from Bio import SeqIO
    import os
    seq_list = {}
    for record in SeqIO.parse(open(faa_fn),'fasta'):
        print(record.id)
        if record.id not in seq_list.keys():
            seq_list[record.id] = str(record.seq)
    return seq_list



def rename_bla_fns(dir="."):
    import glob
    import sys
    
    bla_fns = glob.glob("*.bla")
    for bla_fn in bla_fns:
        new_bla_fn = "+".join([".".join(bla_fn.split(".")[0:2]), ".".join(bla_fn.split(".")[2:4])])
        cmd = "mv " + bla_fn + " " + new_bla_fn +".bla"
        os.system(cmd)

        