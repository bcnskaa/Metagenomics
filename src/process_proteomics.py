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


"""
# extract protein ids from binned faa
rm protein_id2bin_id.map;for f in *.faa;do bin_id=${f};bin_id=${bin_id/.faa/};bin_id=${bin_id/GZ-Cell_Y2./GC2_};ids=$(cat $f | grep -e "^>" | sed "s/>//" |  cut -d" " -f1 | sed "s/scaffold_/GC2_/");for id in ${ids[*]};do echo -e "$id\t$bin_id" >> protein_id2bin_id.map;done;done

"""

def map_protein_ids2bin_id(id_fn, map_fn="protein_id2bin_id.map"):
    with open(map_fn) as IN:
        protein_id2bin_id_map = IN.read().splitlines()
    protein_id2bin_id_map = {p.split("\t")[0]:p.split("\t")[1] for p in protein_id2bin_id_map}
    
    with open(id_fn) as IN:
        ids = IN.read().splitlines() 
        
    OUT = open(id_fn + ".mapped", "w")
    for id in ids:
        if id in protein_id2bin_id_map.keys():
            OUT.write(id + "\t" + protein_id2bin_id_map[id]+"\n")
        else:
            OUT.write(id + "\t" + "Not binned\n")        
    OUT.close()
        
    
    
    