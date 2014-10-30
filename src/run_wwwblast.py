# import os.path
# import sys
# import getopt
# import csv
# import operator
# import time
# from collections import defaultdict
# from subprocess import call
# 
# from Bio import SeqIO
# from Bio import SeqUtils
# from Bio import Seq
# from Bio.Blast import NCBIXML
# from Bio.Blast import NCBIWWW
# 
# 
# 
# def main(argv):
#     # Import blast
#     bitscore_threshold = 1000
#     
#     blast_fn = "SWH-GZ.bla"
#     res = import_blast(blast_fn)
#     res.sort(key=lambda k : k[3], reverse=True)
#     res2 = [r for r in res if r[3] > bitscore_threshold]
#     
#     res_ids = [r[0] for r in res2]
#     res_ids = list(set(res_ids))
#     
#     fasta_fn = "SWH_scaffold.fa"
#         
#     seq_db = SeqIO.index(fasta_fn, "fasta")
#     
#     ids = import_ids(id_fn)
#     
#     if len(ids) > 0:
#         return 0
#         
#     selected_seqs = [seq_db[id] for id in seq_db.keys() if id in res_ids]
#     
#     for s in selected_seqs:
#         sid = s.id
#         seq = str(s.seq)
#         records = run_ncbi_wwwblast(seq)     
#         candidate_species = []
#         if records is not None:
#             if len(records[0].alignments) > 0:    
#                 for alignment in records[0].alignments:
#                     if alignment.hsps[0].score > blast_bitscore_threshold:
#                         species = str((alignment.title).split("|")[4])
#                         species = species.strip()
#                         candidate_species.append(species)
#                         print(sid + ":" + species)
# 
# 
# 
# 
# def import_blast(blast_fn):
#     with open(blast_fn) as IN:
#         lines = IN.read().splitlines()
#     blast_res = [[res.split()[0], res.split()[1], int(res.split()[3]), float(res.split()[11])] for res in lines]
#     return blast_res
# 
# 
#     
# def run_ncbi_wwwblast(seq, database="nr", blast_program="blastn", expect=10, hitlist_size=50):
#     print_status("Connecting to NCBI Web BLAST: " + blast_program + " e-value=" + str(expect) + " hitlist_size=" + str(hitlist_size))
#     result_handler = NCBIWWW.qblast(program=blast_program, database=database, sequence=seq, expect=expect, hitlist_size=hitlist_size)
#     # Store the wwwblast results to local memory
#     records = list(NCBIXML.parse(result_handler))
#     
#     return records
# 
# 
#     
# def import_ids(id_fn):
#     with open(id_fn) as IN:
#         ids = IN.read().splitlines()    
#     return ids
#     
# 
# # Invoke the main function
# if __name__ == "__main__":
#     main(sys.argv[1:])



import os.path
import sys
import getopt
import csv
import operator
import time
from collections import defaultdict
from subprocess import call

from Bio import SeqIO
from Bio import SeqUtils
from Bio import Seq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


def main(argv):
    if len(argv) < 2 or len(argv) > 3:
        print("Usage: python run_wwwblast.py fasta_fn output_fn [blast_fn] ")
        return
    
    fasta_fn = argv[0]
    output_fn = argv[1]
    blast_fn = None
    
    if len(argv) == 3:
        blast_fn = argv[2]
     
    print("Reading from " + fasta_fn)
    
    if blast_fn is not None:
        candidates = process_selected_by_blast(blast_fn, fasta_fn)
    else:
        candidates = process_all(fasta_fn)
    
    if len(candidates) > 0:
        with open(output_fn, "w") as OUT:
            for k in candidates.keys():
                s = "\n".join(candidates[k])
                if len(s) > 0:
                    OUT.write(s + "\n")
        print("Number of items processed: " + str(len(candidates)))
    
    

def process_all(fasta_fn, selected_ids=None, blast_bitscore_threshold=500, seq_len_threshold=1000):
    seq_db = SeqIO.index(fasta_fn, "fasta")
    
    if selected_ids is not None:
        selected_seqs = [seq_db[id] for id in seq_db.keys() if id in res_ids]
    else:
        selected_seqs = [seq_db[id] for id in seq_db.keys()]
    
    print("Number of sequences to be processed: " + str(len(selected_seqs)))
    
    idx = 0
    candidates = {}
    for s in selected_seqs:
        sid = s.id
        seq = str(s.seq)
        
        if len(seq) < seq_len_threshold:
            continue
        
        idx = idx + 1
        print("Blasting [" + str(idx) + " out of " + str(len(selected_seqs)) + "] " + sid + " (" + str(len(seq)) + " bp).")
        
#         if idx == 3:
#             break
        
        records = run_ncbi_wwwblast(seq)     
        candidate_species = []
        if records is not None:
            if len(records[0].alignments) > 0:    
                for alignment in records[0].alignments:
                    if alignment.hsps[0].score > blast_bitscore_threshold:
                        species = str((alignment.title).split("|")[4])
                        species = species.strip()
                        #candidate_species.append(species)
                        candidate_species.append(sid + ":" + str(len(seq)) + ":" + species + ":" + str(alignment.length) + ":" + str(alignment.hsps[0].score))
                        print(sid + ":" + str(len(seq)) + ":" + species + ":" + str(alignment.length) + ":" + str(alignment.hsps[0].score))
        
        candidates[sid] = candidate_species

    return candidates
 
    

def process_selected_by_blast(blast_fn, fasta_fn, bitscore_threshold=1000, blast_bitscore_threshold=600):
    # Import blast
    #bitscore_threshold = 1000
    #blast_bitscore_threshold = 600
    
    
    res = import_blast(blast_fn)
    res.sort(key=lambda k : k[3], reverse=True)
    res2 = [r for r in res if r[3] > bitscore_threshold]
    
    res_ids = [r[0] for r in res2]
    res_ids = list(set(res_ids))
    
    #fasta_fn = "SWH_scaffold.fa"
    
    
    candidates = process_all_seqs(fasta_fn, res_ids)
#     seq_db = SeqIO.index(fasta_fn, "fasta")
#     
#         
#     selected_seqs = [seq_db[id] for id in seq_db.keys() if id in res_ids]
#     candidates = {}
#     for s in selected_seqs:
#         sid = s.id
#         seq = str(s.seq)
#         records = run_ncbi_wwwblast(seq)     
#         candidate_species = []
#         if records is not None:
#             if len(records[0].alignments) > 0:    
#                 for alignment in records[0].alignments:
#                     if alignment.hsps[0].score > blast_bitscore_threshold:
#                         species = str((alignment.title).split("|")[4])
#                         species = species.strip()
#                         candidate_species.append(species)
#                         print(sid + ":" + species)
#     candidates[sid] = candidate_species

    return candidates



def import_blast(blast_fn):
    with open(blast_fn) as IN:
        lines = IN.read().splitlines()
    blast_res = [[res.split()[0], res.split()[1], int(res.split()[3]), float(res.split()[11])] for res in lines]
    return blast_res


    
def run_ncbi_wwwblast(seq, database="nr", blast_program="blastn", expect=1e-5, hitlist_size=5):
    result_handler = NCBIWWW.qblast(program=blast_program, database=database, sequence=seq, expect=expect, hitlist_size=hitlist_size)
    # Store the wwwblast results to local memory
    records = list(NCBIXML.parse(result_handler))
    
    return records


    
def import_ids(id_fn):
    with open(id_fn) as IN:
        ids = IN.read().splitlines()    
    return ids
    

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])

