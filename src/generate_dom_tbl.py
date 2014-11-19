import glob
import process_HMM
import os

hmm_id = "dbCAN"
fns = glob.glob("*_5000/Markers/" + hmm_id + "/*.dom.tbl")
#fns = glob.glob("*_5000/Markers/bins/" + hmm_id + "/*.dom.tbl")
ids = []
hmm_score_threshold=200

dom_tbls = {}
for fn in fns:
   hmm_orf_dict = process_HMM.postprocess_HMMER_search(fn, hmm_score_threshold=hmm_score_threshold)
   id = (os.path.basename(fn)).replace(".dom.tbl", "")
   id = id.replace("."+hmm_id, "")
   ids.append(id)
   dom_tbl = process_HMM.generate_dom_tbl(hmm_orf_dict)
   dom_tbls[id] = dom_tbl
   

dom_ids = []
for k in dom_tbls.keys():
    dom_ids.extend(dom_tbls[k].keys())
dom_ids = list(set(dom_ids))


dom_list = {}
dom_list["header"] = [0 for i in range(0, len(fns))]
for dom_id in dom_ids:
    dom_list[dom_id] = [0 for i in range(0, len(fns))]
    
for idx, id in enumerate(ids):
    dom_list["header"][idx] = id
    print("Processing " + id)
    for dom_id in dom_ids:
        if dom_id in dom_tbls[id].keys():
            #print("Dom=" + dom_id)
            seq_ids = [d[0] for d in dom_tbls[id][dom_id]]
            seq_ids = list(set(seq_ids))
            #dom_list[dom_id][idx] = len(dom_tbls[id][dom_id])
            dom_list[dom_id][idx] = len(seq_ids)

with open("/home/siukinng/db/Markers/" + hmm_id + "/" + hmm_id + ".hmm.info") as IN:
    hmm_info = IN.read().splitlines()
hmm_info = {v.split("\t")[0]:v.split("\t")[1] for v in hmm_info}
    
# if hmm_id == "specI":
#     with open("/home/siukinng/db/Markers/specI/specI.hmm.info") as IN:
#         specI_info = IN.read().splitlines()
#     specI_info = {v.split("\t")[0]:v.split("\t")[1] for v in specI_info}

with open(hmm_id + ".tbl", "w") as OUT:
    OUT.write("HEADER\tDESC\t" + "\t".join(dom_list["header"]) + "\n")
    for dom_id in dom_ids:
        if dom_id in hmm_info.keys():
            info = hmm_info[dom_id]
        else:
            info = dom_id
        OUT.write(dom_id + "\t" + info + "\t" + "\t".join(str(v) for v in dom_list[dom_id]) + "\n") 

"""


"""