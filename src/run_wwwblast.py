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


MAX_SEQ_LEN = 5000
blast_program = "blastn"
blast_db = "nr"

def print_usage():
    print("A simple python script for doing NCBI wwwblast")
    print(" ")
    print("Usage:")
    print("  python " + sys.argv[0] + " -i FASTA_INFILE [-o WWWBLAST_RESULT_OUTFILE] [-b BLAST-INFILE [-1|-2]] [-p] [-q QID-LIST-INFIL] [-l BLACKLIST-INFILE] [-h]")
    print("      -i STRING  Input file containing sequences in FASTA format")
    print("      -o STRING  Output file [optional]")
    print("                 If no output file is provided, result will be exported to a file with the file name of input file plus a suffix .bla.")
    print("      -q STRING  Input file containing fasta header [optional]")
    print("      -b STRING  Input file of blast result (in tabular m6 format) [optional]")
    print("      -1 | -2    blast ids to be used -1=query id or -2=subject id [optional]")
    print("      -l STRING  File name containing ids to be excluded from wwwblast [optional]")
    print("      -p         Input seqeunce(s) is(are) protein sequence and ask wwwblast to do a blastp search [optional, default: nucleotide sequence and blastn]")
    print("      -d STRING  Blast database to search against. [optional, default: nr]")
    print("                 Available databases: nr, ref_seq")
    print("      -h         Print usage")
    print(" ")
    print(" ") 
    print("Ver 0.1")



def main(argv):
    try:
        opts, args = getopt.getopt(argv,"h12pi:o:b:q:l:")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)

#     if len(argv) < 2 or len(argv) > 4:
#         print("Usage: " + sys.argv[0] + " run_wwwblast.py fasta_fn output_fn [blast_fn] [blacklist]")
#         return
    
    fasta_fn = None
    output_fn = None
    blast_fn = None
    list_fn = None
    blacklist_fn = None
    blacklist = None
    search_query_id = True
   
    global blast_program
    global blast_db
    
    blast_program = "blastn"
    blast_db = "nr"

    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            fasta_fn = arg
        elif opt in ("-o", "--outfile"):
            output_fn = arg
        elif opt == "-q":
            list_fn = True
            print "IDs will be imported from " + list_fn
        elif opt == "-b":
            blast_fn = arg       
            print "Search will be based on BLAST results from " + blast_fn + "."
        elif opt == "-1":
            search_query_id = True
        elif opt == "-2": 
            search_query_id = False
        elif opt == "-p":
            print "Search protein sequences"
            blast_program = "blastp"
        elif opt == "-l":          
            blacklist_fn = arg
            with open(blacklist_fn) as IN:
                blacklist = IN.read().splitlines()
            blacklist = list(set(blacklist))
            print "Blacklist imported: " + str(len(blacklist)) + " IDs"
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)

    if fasta_fn is None:
        print "Missing input file"
        print_usage()
        sys.exit(0)
    
    if output_fn is None:
        if list_fn is not None:
            output_fn = list_fn + ".bla"
        elif blast_fn is not None:
            output_fn = blast_fn + ".bla"
        else:
            output_fn = fasta_fn + ".bla"
        print "Result will be exported to " + output_fn
    
    
#     if len(argv) == 3:
#         blast_fn = argv[2]
     
     
    if list_fn is not None:
        with open(lst_fn) as IN:
            ids = IN.read().splitlines()
            
        ids = list(set(ids))
        
        excluded_id_count = len(ids)
        if blacklist is not None:
            ids = [id for id in ids if id not in blacklist]
        print str(excluded_id_count - len(ids)) + " IDs are excluded."
        
        candidates = process_all(fasta_fn, ids)
    elif blast_fn is not None:
        opt = "qid"
        if search_query_id:
            print "Search will be based on Query IDs."
        else:
            print "Search will be based on Subject IDs."
            opt = "sid"
            
        candidates = process_selected_by_blast(blast_fn, fasta_fn, selected_id=opt, excluded_ids=blacklist)
    else:
        candidates = process_all(fasta_fn, excluded_ids=blacklist)
    

    if len(candidates) > 0:
        with open(output_fn, "w") as OUT:
            for k in candidates.keys():
                s = "\n".join(candidates[k])
                if len(s) > 0:
                    OUT.write(s + "\n")
        print("Number of items processed: " + str(len(candidates)))
    
    

def process_all(fasta_fn, selected_ids=None, excluded_ids=None, blast_bitscore_threshold=500, seq_len_threshold=50):
    seq_db = SeqIO.index(fasta_fn, "fasta")
    
    global blast_program
    
    if selected_ids is not None:
        selected_seqs = [seq_db[id] for id in seq_db.keys() if id in selected_ids]
    elif excluded_ids is not None:
        selected_seqs = [seq_db[id] for id in seq_db.keys() if id not in excluded_ids]   
    else:
        selected_seqs = [seq_db[id] for id in seq_db.keys()]
    
    print("Number of sequences to be processed: " + str(len(selected_seqs)))
    
    idx = 0
    candidates = {}
    for s in selected_seqs:
        sid = s.id
        seq = str(s.seq)
        
        if len(seq) < seq_len_threshold:
            print(sid + " is shorter than " + seq_len_threshold + " characters, skipped.")
            continue
        
        print("Blasting [" + str(idx) + " out of " + str(len(selected_seqs)) + "] " + sid + " (" + str(len(seq)) + " bp).")
        idx = idx + 1
                
        if len(seq) > MAX_SEQ_LEN:
            seq = seq[0:MAX_SEQ_LEN]
            
        records = run_ncbi_wwwblast(seq, blast_program=blast_program)     
        candidate_species = []
        if records is not None:
            if len(records[0].alignments) > 0:    
                for alignment in records[0].alignments:
                    if alignment.hsps[0].score > blast_bitscore_threshold:
                        species = str((alignment.title).split("|")[4])
                        species = species.strip()
                        #candidate_species.append(species)
                        candidate_species.append(sid + ":" + str(len(seq)) + ":" + species + ":" + str(alignment.length) + ":" + str(alignment.hsps[0].score) + ":" + str(alignment.hsps[0].identities))
                        print(sid + ":" + str(len(seq)) + ":" + species + ":" + str(alignment.length) + ":" + str(alignment.hsps[0].score) + ":" + str(alignment.hsps[0].identities))
        
        candidates[sid] = candidate_species
        
#         if idx == 4:
#             break

    return candidates
 


def process_selected_by_blast(blast_fn, fasta_fn, selected_id="qid", excluded_ids=None, bitscore_threshold=1000, blast_bitscore_threshold=600):
    # Import blast
    #bitscore_threshold = 1000
    #blast_bitscore_threshold = 600
    
    
    res = import_blast(blast_fn)
    res.sort(key=lambda k : k[3], reverse=True)
    res2 = [r for r in res if r[3] > bitscore_threshold]
    
    if selected_id == "qid":
        res_ids = [r[0] for r in res2]
    else:
        res_ids = [r[1] for r in res2]

    # Create a unique list
    res_ids = list(set(res_ids))
    
    excluded_id_count = len(res_ids)
    if excluded_ids is not None:
        res_ids = [id for id in res_ids if id not in excluded_ids]
        print str(excluded_id_count - len(res_ids)) + " IDs are excluded."
    
    candidates = process_all(fasta_fn, res_ids)
    
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


    
def run_ncbi_wwwblast(seq, seq_id=None, database="nr", blast_program="blastn", expect=1e-5, hitlist_size=5, timeout=50):
    records = None
    
    try:
         result_handler = NCBIWWW.qblast(program=blast_program, database=database, sequence=seq, expect=expect, hitlist_size=hitlist_size)
         # Store the wwwblast results to local memory
         records = list(NCBIXML.parse(result_handler))
    except:
        if seq_id is not None:
            print "Unable to search " + seq_id
        else:
            print "Unable to search the sequence"
    return records


    
def import_ids(id_fn):
    with open(id_fn) as IN:
        ids = IN.read().splitlines()    
    return ids
    

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])


"""


# For processing MW's dataset

with open("GZ_contig-SWH_contig.bla.l3000p95.0.bla.bla") as IN:
    qid = IN.read().splitlines()
qid = {q.split(":",1)[0]:q.split(":")[2] for q in qid}
    
    
with open("qid.nr") as IN:
    gp6 = IN.read().splitlines()
gp6 = {g:"Methanosaeta concilii GP-6, complete genome" for g in gp6}

combined_qids = dict(gp6.items() + qid.items())

with open("GZ_contig-SWH_contig.bla.l3000p95.0.bla") as IN:
    bla = IN.read().splitlines()
bla = {b.split("\t")[0]:b for b in bla}


with open("../GZ/Binning/MaxBin/GZ.abund") as IN:
    gz_abund = IN.read().splitlines()
gz_abund = {g.split("\t")[0]:g.split("\t")[1] for g in gz_abund}

with open("../SWH/Binning/MaxBin/SWH.abund") as IN:
    swh_abund = IN.read().splitlines()
swh_abund = {s.split("\t")[0]:s.split("\t")[1] for s in swh_abund}


for b in bla.keys():
    qid = bla[b].split("\t")[0]
    sid = bla[b].split("\t")[1]
    if b in combined_qids.keys():
        bla[b] = bla[b] + "\t" + combined_qids[b]
    else:
        bla[b] = bla[b] + "\t" + "unclassified"
    bla[b] = bla[b] + "\t" + gz_abund[qid]
    bla[b] = bla[b] + "\t" + swh_abund[sid]
    


with open("GZ_contig-SWH_contig.bla.l3000p95.0.bla.species", "w") as OUT:
    for b in bla:
        OUT.write(bla[b] + "\n")
        

"""


"""
import os
import glob

fns = glob.glob("blast/*.bla")

for fn in fns:
    print("Processing " + fn)
    with open(fn) as IN:
        bla = IN.read().splitlines()
        
    map_fn="/home/siukinng/db/Markers/specI/data/RefOrganismDB.v9.taxids2names"
    with open(map_fn) as IN:
        map = IN.read().splitlines()
    map = {m.split("\t")[0]: m.split("\t")[1] for m in map}
    
    out_fn = fn + ".mapped"
    with open(out_fn, "w") as OUT:
        for b in bla:
            b = b.split("\t")
            id = b[1].split(".")[0]
            if id in map.keys():
                b[1] = b[1] + ":" + map[id]
            else:
                b[1] = b[1] + ":" + "Unknown"
            OUT.write("\t".join(b) + "\n")
        
"""

"""
for f in ../../../*_5000/*.fasta;do
    id=${f##*/};
    id=${id/.fasta/};
    echo $id;
    for bla_fn in *.bla.mapped;do
        cat $bla_fn | grep -e "^$id" | cut -f2 | cut -d':' -f2 >> $id.summary
    done
done

"""

"""
import os
import operator
import glob

summary_ofn = "all.tax"
summary_fns = glob.glob("*.summary")

summary_list = {}
for summary_fn in summary_fns:
    bin_id = summary_fn.replace(".summary", "")
    print("Processing " + bin_id)
    
    with open(summary_fn) as IN:
        summary = IN.read().splitlines()
        
    genus = [s.split(" ")[0] for s in summary]
    
    genus_ids = list(set(genus))
    genus_summary = {id: genus.count(id) for id in genus_ids}
    
    sorted_genus = sorted(genus_summary.items(), key=operator.itemgetter(1))
    sorted_genus.reverse()
    summary_list[bin_id] = sorted_genus[0:10]

bin_ids = summary_list.keys()
bin_ids = sorted(bin_ids)

with open(summary_ofn, "w") as OUT:
    for bin_id in bin_ids:
        bin_summary = summary_list[bin_id]
        summary = [g[0] + ":" +  str(g[1]) for g in bin_summary]
        OUT.write(bin_id + "\t" + "\t".join(summary) + "\n")
        

"""

"""
with open("TIGRFAM_ROLES", "w") as OUT:
    for tigrfam_id in role_links.keys():
        role_id = role_links[tigrfam_id]
        if role_id in role_names.keys():
            OUT.write(tigrfam_id + "\t" + role_names[role_id] + "\n")
        
"""
"""
from __future__ import division
import os
import glob 
import numpy
from Bio import SeqIO

# Calculate bin sequence len
bin_fns = glob.glob("bin_seqs/*.fasta")
bin_lens = {}
for bin_fn in bin_fns:
    print("Getting length from " + bin_fn)
    bin_id = (os.path.basename(bin_fn)).replace(".fasta", "")
    seqs = SeqIO.index(bin_fn, "fasta")
    bin_lens[bin_id] = sum([len(str(seqs[sid].seq)) for sid in seqs])


mtx = {}
bla_fns = glob.glob("*.bla")
for bla_fn in bla_fns:

    id = (os.path.basename(bla_fn)).replace(".bla", "")
    
    qid = ".".join(id.split(".")[0:2])
    sid = ".".join(id.split(".")[2:4])
    
    print("Processing " + id)
    with open(bla_fn) as IN:
        bla_res = IN.read().splitlines()
    if bla_res is None:
        mtx[id] = 0.0
    elif len(bla_res) > 0:
        #mtx[id] = numpy.mean([(int(b.split("\t")[3]) - int(int(b.split("\t")[5]) + int(b.split("\t")[4]))) / int(b.split("\t")[3]) for b in bla_res])
        #mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / ((bin_lens[qid] + bin_lens[sid]) / 2)
        mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / bin_lens[qid]
    else:
        mtx[id] = 0.0
    
with open("bla.mtx", "w") as OUT:
    for k in mtx.keys():
        OUT.write(k + "\t" + k.replace(".", "\t") + "\t" + str(mtx[k]) + "\n")
        
    
"""

"""

import glob

bla_fn = "K00399-mgm4560350.p90.bla"


bla_fns = glob.glob("./*.bla")
for bla_fn in bla_fns:
    process_blast_mapping(bla_fn)


def process_blast_mapping(bla_fn, qid_fn=None, out_fn=None):
    print("Processing " + bla_fn)
    
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
    bla = [b.split("\t") for b in bla]
    sids = [b[1] for b in bla]
    
    
    nr_sids = list(set(sids))
    nr_sid_bla = {}
    for sid in nr_sids:
        cur_idx = -1
        selected_bla = [b for b in bla if b[1] == sid]
        for b in selected_bla:
            if cur_idx == -1:
                nr_sid_bla[sid] = b
                cur_idx = 0
            else:
                if float(b[11]) > float(nr_sid_bla[sid][11]):
                    #print(float(b[11]), " > ", float(nr_sid_bla[sid][11]))
                    nr_sid_bla[sid] = b 
    
    # Import qids
    if qid_fn is None:
        qid_fn = bla_fn.split("-")[0] + ".fasta"
    
    from Bio import SeqIO
    qid_seqs = SeqIO.index(qid_fn, "fasta")
    #qid_seq_ids = list([k for k in qid_seqs.keys()])
    qid_seq_descriptions = list([qid_seqs[k].description for k in qid_seqs.keys()])
    qid_seq_ids = {(qid_seqs[k].description).split(" ", 1)[0]:(((qid_seqs[k].description).split(" ", 1)[1]).split(" GN=")[0]).replace(" OS=", "\t")  for k in qid_seqs.keys()}
    
    
    import numpy
    #nr_qids = [nr_sid_bla[k][0] for k in nr_sid_bla.keys()]
    nr_qids = [k for k in qid_seq_ids.keys()]
    nr_qids = list(set(nr_qids))
    nr_qid_count = {}
    nr_qid_mean_iden = {}
    nr_qid_mean_len = {}
    for qid in nr_qids:
        selected_sid_bla = [nr_sid_bla[k] for k in nr_sid_bla.keys() if nr_sid_bla[k][0] == qid]
        nr_qid_count[qid] = 0
        nr_qid_mean_iden[qid] = 0.0
        nr_qid_mean_len[qid] = 0.0
        
        if len(selected_sid_bla) > 0:
            nr_qid_count[qid] = len(selected_sid_bla)
            nr_qid_mean_iden[qid] = numpy.mean([float(nr_sid_bla[k][2]) for k in nr_sid_bla.keys() if nr_sid_bla[k][0] == qid])
            nr_qid_mean_len[qid] = numpy.mean([float(nr_sid_bla[k][3]) for k in nr_sid_bla.keys() if nr_sid_bla[k][0] == qid])
            #nr_qid_mean_iden[qid] = numpy.mean([float(nr_sid_bla[k][2]) for k in nr_sid_bla.keys() if nr_sid_bla[k][0] == qid])
            #nr_qid_mean_len[qid] = numpy.mean([float(nr_sid_bla[k][3]) for k in nr_sid_bla.keys() if nr_sid_bla[k][0] == qid])
    
    if out_fn is None:
        out_fn = bla_fn + ".summary"
    
    print("Exporting " + str(len(nr_qid_count)) + " results to " + out_fn)
    header = ["Protein_ID", "Protein_Name", "Species", "Hit_Count_by_Read", "Mean_Hit_Identity", "Mean_Hit_Length"]
    with open(out_fn, "w") as OUT:
        OUT.write("\t".join(header) + "\n")
        for qid in nr_qids:
            summary = [qid, qid_seq_ids[qid], str(nr_qid_count[qid]), str(nr_qid_mean_iden[qid]), str(nr_qid_mean_len[qid])]
            OUT.write("\t".join(summary) + "\n")

"""
