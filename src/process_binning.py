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






def generate_ref_genomes(nr_sids, out_fn, db_dir=mg_pipeline.NCBI_BIOPROJECT_DB_HOME):
    print_msg("Preparing " + str(len(nr_sids)) +" reference sequences from " + db_dir)
     
    nr_sids = wash_id(nr_sids)
     
    OUT = open(out_fn, "w")
    status = {sid:0 for sid in nr_sids}
    for sid in nr_sids:
        sid_fn = db_dir + "/" + sid + ".fna"
        s = 0
        if os.path.exists(sid_fn):
            seqs = SeqIO.index(sid_fn, "fasta")
            if len(seqs) == 1:
                for seq in SeqIO.parse(sid_fn, "fasta"):
                    SeqIO.write(seq, OUT, "fasta")
            else:
                [exported_n, seq] = merge_fa.merge(sid_fn, sid)
                merge_fa.export_fasta(seq[0], seq[1], OUT)
            s = s + 1
        status[sid] = s
    OUT.close()
    
    print_msg(str(len([s for s in status.keys() if status[s] > 0])) + " out of " + str(len(status)) + " were found.")
    
    return status
    


"""

BLAST the contigs against NCBI BioProjects and determine the closet reference contig sequences can map to.

 
"""
def mapping_contigs_to_references(contig_fn, out_fn=None, out_dir=None, db_fn=mg_pipeline.NCBI_BIOPROJECT_DB, len_cutoff=3000):
    print_msg("Mapping contigs from " + contig_fn + " to " + db_fn)
    
    subject_id = ((os.path.basename(db_fn)[::-1].split(".", 1))[1])[::-1]
    query_id = ((os.path.basename(contig_fn)[::-1].split(".", 1))[1])[::-1]
    
    if out_fn is None:
        if out_dir is None:
            #out_fn = os.path.dirname(contig_fn) + "/" + query_id
            out_fn = "./" + query_id
        else:
            out_fn = out_dir + "/" + query_id
        out_fn = out_fn + "+" + subject_id + ".bla" 
    else:
        out_dir = os.path.dirname(out_fn)
        if len(out_dir) == 0:
            if out_dir is None:
                out_fn = "./" + out_fn
            else:
                out_fn = out_dir + "/" + out_fn
        

    print_msg("Result will be exported to " + out_fn)
    
    out_fn = mg_pipeline.blastn(contig_fn, db_fn, outdir=out_dir, outfn=out_fn, perc_identity=85, max_target_seqs=1, evalue=1e-30)
    
    blast_res = filter_blast_res.import_blast_res(out_fn, thres_aln_len = len_cutoff)
    
    nr_sids = filter_blast_res.get_nr_ids(blast_res, is_query_id_selected=False)
    
    return nr_sids
    


def wash_id(ids):
    washed_ids = [id.split(".")[0] for id in ids]
    return washed_ids
    
    
    
def print_msg(msg):
    mg_pipeline.print_status(msg)