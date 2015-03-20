from Bio import SeqIO
import glob


import sys
import os
sys.path.append(os.path.abspath("/home/siukinng/tools/scripts"))

import mg_pipeline


"""

"""
def extract_combined_hmm_seq(hmm_id="specI", dir_suffix="_5000"):
    faa_ext = "faa"
    #hmm_id = "specI"

    #dir_suffix = "_5000"    
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
import extract_hmm_seq
extract_hmm_seq.extract_bin_hmm_seq(dir_suffix = "_2000")

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
def export_hmm_profile(profile_acc_id, out_fn=None):
    hmm = extract_hmm_profile_from_db(profile_acc_id)
    if hmm is not None:
        if out_fn is None:
            out_fn = profile_acc_id + ".hmm"
        
        with open(out_fn, "w") as OUT:
            OUT.write(hmm + "//\n")
        
    
        
        

def extract_hmm_profile_from_db(profile_acc_id, hmm_db=None, hmm_db_fn="/home/siukinng/db/Markers/Pfam/Pfam.hmm"):
    import re
    
    if hmm_db is None:
        hmm_db = import_hmm_db(hmm_db_fn)
    
    if profile_acc_id in hmm_db.keys():
        return hmm_db[profile_acc_id]
    else:
        return None
    


def import_hmm_db(hmm_db_fn="/home/siukinng/db/Markers/Pfam/Pfam.hmm", excluding_acc_version=True):
    import re
    
    print("Importing HMM Database")
    
    hmm_db = {}
    with open(hmm_db_fn) as IN:
        hmms = IN.read().split("//\n")
    del hmms[len(hmms) - 1]
    hmm_db = {re.findall("\nACC\s+(.+)\n", hmm)[0].split(".")[0]:hmm for hmm in hmms}
    
    return hmm_db

  

"""
for f in *.phy;do echo "Processing $f";echo -e "$f\nP\nP\nY\n" > current.cmd; ~/tools/phylip/exe/protdist < current.cmd > screenout; mv outfile $f.mtx;done
    
"""
def split_hmm_profiles(outdir=".", hmm_db_fn="/home/siukinng/db/Markers/Pfam/Pfam.hmm"):
    import re
    
    with open(hmm_fn) as IN:
        hmm = IN.read().split("//\n")
    
    for h in hmm:
        hmm_id = re.findall("\nACC\s+(.+?).\w+\n", h)
        if len(hmm_id) > 0:
            hmm_id = hmm_id[0]
            with open(outdir + "/" + hmm_id + ".hmm", "w") as OUT:
                OUT.write(h + "//")
        



def prepare_pfam_db_info(out_fn="/home/siukinng/db/Markers/Pfam/Pfam.hmm.info", pfam_db_fn="/home/siukinng/db/Markers/Pfam/Pfam.hmm"):
    import re
    with open(pfam_db_fn) as IN:
        hmms = IN.read()
    hmms = hmms.split("//\n")
    del hmms[len(hmms) - 1]
    
    #hmm_map = {}
    OUT = open(out_fn, "w")
    for hmm in hmms:
        hmm_name = re.findall("\nNAME\s+(.+)\n", hmm)[0]
        hmm_acc = re.findall("\nACC\s+(.+)\n", hmm)[0]
        hmm_desc = re.findall("\nDESC\s+(.+)\n", hmm)[0]
        OUT.write(hmm_name + "\t" + hmm_acc + "|" + hmm_desc + "\n")
    OUT.close()
        
    
"""

hmm_id=""
hmm_fn="$hmm_id.hmm"
tc_val=$(cat cohesin/PF00963.hmm | grep -e "TC" | tr -s " " | cut -d" " -f2)
for f in ../../../*_2000/Prodigal/bins/*.faa;do bin_id=${f##*/};bin_id=${bin_id/.faa/};echo "~/tools/hmmer/bin/hmmsearch --max -T $tc_val -o $bin_id-$hmm_id.out -A $bin_id-$hmm_id.aln --tblout $bin_id-$hmm_id.tbl --domtbl $bin_id-$hmm_id.dom.tbl --pfamtblout $bin_id-$hmm_id.pfam.tbl $hmm_fn $f";done

"""