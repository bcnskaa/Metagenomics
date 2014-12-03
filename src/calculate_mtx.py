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
        mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / ((bin_lens[qid] + bin_lens[sid]) / 2)
        #mtx[id] = sum([int(b.split("\t")[3]) - (int(b.split("\t")[5]) + int(b.split("\t")[4])) for b in bla_res]) / bin_lens[qid]
    else:
        mtx[id] = 0.0
    
with open("bla.mtx", "w") as OUT:
    for k in mtx.keys():
        OUT.write(k + "\t" + k.replace(".", "\t") + "\t" + str(mtx[k]) + "\n")
