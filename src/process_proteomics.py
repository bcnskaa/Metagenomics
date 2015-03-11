import glob
import os


def map_sid_to_gi(bla_fn, out_fn=None, gi_fn="/home/siukinng/db/BacteriaDB/all_faa.gi.lst"): 
    print("Reading from " + gi_fn)
    with open(gi_fn) as IN:
        gi = IN.read().splitlines()
    gi_db = {g.split("\t")[0]: g.split("\t")[3] + ":" + g.split("\t")[2] for g in gi if len(g.split("\t")) == 4}
    
    
    print("Reading from " + bla_fn)
    bla_lst = {}
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
    for b in bla:
        qid = b.split("\t")[0]
        sid = b.split("\t")[1]
        if qid not in bla_lst.keys():
            bla_lst[qid] = []
        bla_lst[qid].append(sid)
    
    if out_fn is None:
        out_fn = bla_fn + ".summary"
       
    print("Mapping...")
    OUT = open(out_fn,"w")
    
    for qid in bla_lst.keys():
        print("Processing " + qid)
        gi_ctx = []
        for sid in bla_lst[qid]:
            gi_ctx.append(gi_db[sid])
        OUT.write(qid + "\t" + "\t".join(gi_ctx) + "\n")
    
    OUT.close()
    
