import mg_pipeline
from Bio import SeqIO
import os
import glob

faa_ext = "faa"
hmm_id = "specI"

dir_suffix = "_5000"    
sample_ids = glob.glob("*" + dir_suffix)
sample_ids = [id.replace(dir_suffix, "") for id in sample_ids]

for sample_id in sample_ids:
    sample_dir = "./" + sample_id + dir_suffix
    hmm_dir = sample_dir + "/Markers/" + hmm_id
    faa_fn = sample_dir + "/Prodigal/" + sample_id + ".combined.prodigal.faa"
    seq_outdir = sample_dir + "/DistMtx/" + hmm_id 
    if not os.path.exists(seq_outdir):
        os.makedirs(seq_outdir)
    hmm_orf_dict = mg_pipeline.postprocess_HMMER_search(hmm_dir, dom_overlapping_threshold=40, hmm_score_threshold=150.0)
    mg_pipeline.extract_hmm_seq(hmm_orf_dict, faa_fn, output_dir=seq_outdir, replace_id="scaffold", sample_id=sample_id)


"""
 ls *.faa | 
"""

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