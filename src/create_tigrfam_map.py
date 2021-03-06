# create_tigrfam_map.py
import os


with open("TIGRFAMs_ROLE_NAMES" ) as IN:
    lines = IN.read().splitlines()

mainroles = {l.split("\t")[1] : l.split("\t")[3] for l in lines if l.split("\t")[2] == "mainrole:"}
subroles = {l.split("\t")[1] : l.split("\t")[3] for l in lines if l.split("\t")[2] == "sub1role:"}


with open("TIGRFAMs_ROLE_LINK") as IN:
    lines = IN.read().splitlines()

tigrfam_roles = {l.split("\t")[0] : l.split("\t")[1] for l in lines}


with open("TIGRFAMs_GO_LINK") as IN:
    tigrfam_go_lines = IN.read().splitlines()

#tigrfam_go = {l.split("\t")[0] : l.split("\t")[1] for l in tigrfam_go}
tigrfam_go = {}
for l in tigrfam_go_lines:
    tid = l.split("\t")[0]
    go_id = l.split("\t")[1]
    if tid in tigrfam_go.keys():
        tigrfam_go[tid].append(go_id)
    else:
        tigrfam_go[tid] = [go_id]


with open("kegg2go") as IN:
    kegg2go_lines = IN.read().splitlines()

del kegg2go_lines[0]
del kegg2go_lines[0]
#kegg2go = {l.split(" ; ")[1] + ":" + (l.split(" > ")[0]).strip() : (l.split(" > ")[0]).strip() for l in kegg2go_lines}
#kegg2go = {l.split(" ; ")[1] : [(l.split(" > ")[0]).strip()] for l in kegg2go_lines}
kegg2go = {}
for l in kegg2go_lines:
    go_id = l.split(" ; ")[1]
    kegg_id = (l.split(" > ")[0]).strip()
    if go_id in kegg2go.keys():
        kegg2go[go_id].append(kegg_id)
    else:
        kegg2go[go_id] = [kegg_id]

tigrfam_info = {}
# Populate the tigrfam_info dict
for tid in tigrfam_roles.keys():
    tigrfam_info[tid] = []

for tid in tigrfam_info.keys():
    role_id = tigrfam_roles[tid]
    
    if role_id in mainroles.keys():
        tigrfam_info[tid].append(mainroles[role_id])
    else:
        tigrfam_info[tid].append("Unknown")
        
    if role_id in subroles.keys():
        tigrfam_info[tid].append(subroles[role_id])
    else:
        tigrfam_info[tid].append("Unknown")
    
    if tid in tigrfam_go.keys():
        go_ids = tigrfam_go[tid]
        kegg_terms = ""
        for go_id in go_ids:
            #kegg2go_id = tid + ":" + go_id
            if go_id in kegg2go.keys():
                if len(kegg_terms) == 0:
                    kegg_terms = ",".join(kegg2go[go_id])
                else:
                    kegg_terms = kegg_terms + "," + ",".join(kegg2go[go_id])
        if len(kegg_terms) == 0:
            kegg_terms = "NA"
        tigrfam_info[tid].append(kegg_terms)
    else:
        tigrfam_info[tid].append("NA")
    #tigrfam_info[tid].append([mainroles[role_id], subroles[role_id]], kegg2go[go_id])


with open("TIGRFAM.info", "w") as OUT:
    for tid in tigrfam_info.keys():
        OUT.write(tid + "\t" + "\t".join(tigrfam_info[tid]) + "\n")


"""
TIGR01981 -> 76 mainrole: Biosynthesis of cofactors, prosthetic groups, and carriers, sub1role: other
GO:0005198
GO:0016226

"""

"""

with open("kegg2go") as IN:
    kegg2go_lines = IN.read().splitlines()

kegg2go_lines = [l for l in kegg2go_lines if not l.startswith("!")]
kegg2go = {}
for l in kegg2go_lines:
    go_id = l.split(" ; ")[1]
    kegg_id = ((l.split(" > ")[0]).strip()).replace("KEGG:", "")
    if go_id in kegg2go.keys():
        kegg2go[go_id].append(kegg_id)
    else:
        kegg2go[go_id] = [kegg_id]

        

with open("tigrfams2go") as IN:
    tigrfams2go_lines = IN.read().splitlines()

kegg_match = []

tigrfams2go_map = {}
tigrfams_names = {}
go_map = {}
for l in tigrfams2go_lines:
    if not l.startswith("!"):
        terms = l.split(" > ")
        (tid, name) = terms[0].split(" ", 1)
        tid = tid.replace("JCVI_TIGRFAMS:", "")
        (go_term, go_id) = terms[1].split(" ; ")
        go_term = go_term.replace("GO:", "")
        go_map[go_id] = go_term
        kegg_id = ""
        if go_id in kegg2go.keys():
            kegg_id = "|".join(kegg2go[go_id])
            kegg_match.append(tid)
            
        if tid not in tigrfams2go_map.keys():
            if len(kegg_id) > 0:
                tigrfams2go_map[tid] = [go_id + ":" + go_term + ":" + kegg_id]
            else:
                tigrfams2go_map[tid] = [go_id + ":" + go_term]
        else:
            tigrfams_names[tid] = name
            if len(kegg_id) > 0:
                tigrfams2go_map[tid].append(go_id + ":"+ go_term + ":" + kegg_id)
            else:
                tigrfams2go_map[tid].append(go_id + ":"+ go_term)

     
"""
 
 
"""
import re
from Bio import SeqIO


fn = "TIGR00016"

seqs = SeqIO.index(fn + ".faa", "fasta")
seed_seqs = []
for k in seqs.keys():
    if re.search("Seed", k):
        seed_seqs.append(seqs[k])

with open(fn + ".seed.faa", "w") as OUT:
    for s in seed_seqs:
        SeqIO.write(s, OUT, "fasta")
        
        
"""
