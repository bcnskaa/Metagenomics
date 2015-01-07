import mg_pipeline
from Bio import SeqIO
import os
import glob


"""

"""
def extract_combined_hmm_seq(hmm_id="specI"):
    faa_ext = "faa"
    #hmm_id = "specI"

    dir_suffix = "_5000"    
    sample_ids = glob.glob("*" + dir_suffix)
    sample_ids = [id.replace(dir_suffix, "") for id in sample_ids]
    group_id = "combined"
    hmm_score_threshold = 50.0
    
    hmm_tc_fn = "/home/siukinng/db/Markers/" + hmm_id + "/" + hmm_id + ".tc"
    if not os.path.isfile(hmm_tc_fn):
        hmm_tc_fn = None
        
    for sample_id in sample_ids:
        sample_dir = "./" + sample_id + dir_suffix
        hmm_dir = sample_dir + "/Markers/" + group_id + "/" + hmm_id
        #faa_fn = sample_dir + "/Prodigal/" + sample_id + ".combined.prodigal.faa"
        faa_fn = sample_dir + "/Prodigal/" + sample_id + "." + faa_ext
        seq_outdir = sample_dir + "/DistMtx/" + group_id + "/" + hmm_id 
        
        if not os.path.exists(seq_outdir):
            os.makedirs(seq_outdir)
            
        hmm_orf_dict = mg_pipeline.postprocess_HMMER_search(hmm_dir, dom_overlapping_threshold=40, hmm_score_threshold=hmm_score_threshold, hmm_tc_fn=hmm_tc_fn)
        mg_pipeline.extract_hmm_seq(hmm_orf_dict, faa_fn, output_dir=seq_outdir, replace_id="scaffold", sample_id=sample_id)



"""

"""
def extract_bin_hmm_seq(hmm_id = "specI", dir_suffix = "_5000", hmm_score_threshold=50.0):
    faa_ext = "faa"
    #hmm_id = "Pfam"
    #hmm_score_threshold = 50.0
    hmm_tc_fn = "/home/siukinng/db/Markers/" + hmm_id + "/" + hmm_id + ".tc"
    
    if not os.path.isfile(hmm_tc_fn):
        hmm_tc_fn = None
    
    dom_tbl_suffix = ".dom.tbl"
    #dir_suffix = "_5000"   
    sample_ids = glob.glob("*" + dir_suffix)
    sample_ids = [(os.path.basename(id)).replace(dir_suffix, "") for id in sample_ids]
    
    for sample_id in sample_ids:
        print("Sample ID=" + sample_id)
        sample_dir = "./" + sample_id + dir_suffix
        hmm_dir = sample_dir + "/Markers/bins/" + hmm_id
        bin_hmm_fns = glob.glob(hmm_dir + "/*" + dom_tbl_suffix)
        print("bin_hmm_fns=" + str(len(bin_hmm_fns)))
        for bin_hmm_fn in bin_hmm_fns:
            bin_id = (os.path.basename(bin_hmm_fn)).replace(dom_tbl_suffix, "")
            bin_id = bin_id.replace("." + hmm_id, "")
            seq_outdir = sample_dir + "/DistMtx/bins/" + bin_id + "/" + hmm_id
            if not os.path.exists(seq_outdir):
                os.makedirs(seq_outdir)
            faa_fn = sample_dir + "/Prodigal/bins/" + bin_id + ".faa"
            print("Processing " + bin_id)
            hmm_orf_dict = mg_pipeline.postprocess_HMMER_search_by_fn(bin_hmm_fn, dom_overlapping_threshold=40, hmm_score_threshold=hmm_score_threshold, hmm_tc_fn=hmm_tc_fn)
            mg_pipeline.extract_hmm_seq(hmm_orf_dict, faa_fn, output_dir=seq_outdir, replace_id="scaffold", sample_id=bin_id)



# s = mg_pipeline.summarize_hmm_orf_dict(hmm_orf_dict)
# 
# trim_last_character = True
# 
# import glob
# 
# faa_fn = glob.glob("./*." + faa_ext)
# if len(faa_fn) == 1:
#     faa_fn = faa_fn[0]
# 
# sample_id = faa_fn.replace("./", "")
# sample_id = sample_id.replace("."+faa_ext, "")
# 
# for hmm_id in s.keys():
#     seq_ids = [v[0] for v in s[hmm_id]]
#     seqs = mg_pipeline.pick_seqs(faa_fn, seq_ids)
#     for seq in seqs:
#         #seq.description = seq.id
#         seq.description = ""
#         seq.id = (seq.id).replace("contig-80", sample_id)
#         if trim_last_character:
#             seq.seq = seq.seq[0 : len(seq.seq) - 1]
#         
#     out_fn = hmm_id + "." + faa_ext
#     with open(out_fn, "w") as OUT:
#         SeqIO.write(seqs, OUT, "fasta")  
        
         
        
"""
hmm_orf_dict = mg_pipeline.postprocess_HMMER_search(".", dom_overlapping_threshold=40, hmm_score_threshold=150.0)
s = mg_pipeline.summarize_hmm_orf_dict(hmm_orf_dict)

trim_last_character = True

import glob

faa_fn = glob.glob("./*." + faa_ext)
if len(faa_fn) == 1:
    faa_fn = faa_fn[0]

sample_id = faa_fn.replace("./", "")
sample_id = sample_id.replace("."+faa_ext, "")

for hmm_id in s.keys():
    seq_ids = [v[0] for v in s[hmm_id]]
    seqs = mg_pipeline.pick_seqs(faa_fn, seq_ids)
    for seq in seqs:
        #seq.description = seq.id
        seq.description = ""
        seq.id = (seq.id).replace("contig-80", sample_id)
        if trim_last_character:
            seq.seq = seq.seq[0 : len(seq.seq) - 1]
        
    out_fn = hmm_id + "." + faa_ext
    with open(out_fn, "w") as OUT:
        SeqIO.write(seqs, OUT, "fasta")
"""



"""

for f in *.phy;do echo "Processing $f";echo -e "$f\nP\nP\nY\n" > current.cmd; ~/tools/phylip/exe/protdist < current.cmd > screenout; mv outfile $f.mtx;done

"""


