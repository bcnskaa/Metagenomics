##
# KEGG REST API:
# http://www.kegg.jp/kegg/docs/keggapi.html
#
##
import glob
from itertools import izip


def process_seed():
    
    largest_gi = 815659664 + 1
    gi2taxid = [0] * (largest_gi + 1)
    gi2seed = [0] * (largest_gi + 1) 
    gi2taxid_map_fn = "gi_taxid_prot.dmp"
    gi2seed_map_fn="/home/siukinng/db/Markers/ncbi_nr/mapping/gi2seed.map"
    
    
    ## Import gi2seed_map
    with open(gi2seed_map_fn) as IN:
        s_vals = IN.read().splitlines()
    
    print("Number of rows imported: " + str(len(s_vals)))
    for v in s_vals:
        [gi, seed_id] = map(int, v.split("\t"))
        gi2seed[gi] = seed_id
    
    
    ## Import gi2taxid_map
    with open(gi2taxid_map_fn) as IN:
        vals = IN.read().splitlines()
    
    
    print("Number of rows imported: " + str(len(vals)))
    for v in vals:
        [gi, taxid] = map(int, v.split("\t"))
        gi2taxid[gi] = taxid

    
    m8_fns = glob.glob("/home/siukinng/MG/m8/scaffolds/*.m8")
  
    
    for m8_fn in m8_fns:
        print("Processing " + m8_fn)
        with open(m8_fn) as IN:
            m8 = IN.read().splitlines()
        m8_map = [0] * len(m8)
    
    
        #
        skipped_n = 0
        processed_n = 0
        for i, m in enumerate(m8):
            m = m.split("\t")
            # Query ID
            qid = m[0]
            # GI
            gi = int(m[1].split("|")[1])
            tax_id = ""
            seed_id = ""
            try:
                # tax_id
                tax_id = gi2taxid[gi]
                # seed_id
                seed_id = gi2seed[gi]
                processed_n += 1
            except:
                skipped_n += 1
            m8_map[i] = [qid, gi, tax_id, seed_id]
    
        out_fn = m8_fn + ".map"
        print("Exporting results to " + out_fn )
        with open(out_fn , "w") as OUT:
            OUT.write("#read_id\tgi\ttax_id\tseed_id\n")    
            for m in m8_map:
                OUT.write("\t".join(map(str, m)) + "\n")
        print("Number of processed: " + str(processed_n) + " (" + str(skipped_n) + " skipped)")





"""
import process_m8
import glob


m8_fns = glob.glob("*.m8")

gi2go_map = process_m8.import_gi2go()
for m8_fn in m8_fns:
    mm = process_m8.map_m8_gi2go(m8_fn, m8_fn+".gi2go", gi2go_map)


"""
def import_gi2go(gi2go_map_fn="/home/siukinng/db/Markers/Uniprot/idmapping_selected.tab.gi2go"):
    gi2go_map = []
    with open(gi2go_map_fn) as IN:
        gi2go_map = IN.read().splitlines()
    gi2go_map = [g.split("\t") for g in gi2go_map]
    gi2go_map = {int(g[0]):g[1] for g in gi2go_map if len(g) == 2 and len(g[0]) > 0}
    
    return gi2go_map
    


"""
Given a gi2go_map
"""
def map_m8_gi2go(m8_fn, gi2go_ofn, gi2go_map=None, ):
    if gi2go_map is None:
        print("Importing GI2GO map...")
        gi2go_map = import_gi2go()
    
    print("Processing " + m8_fn + " and exporting to " + map_ofn)
    IN = open(m8_fn, "r")
    OUT = open(gi2go_ofn, "w")
    m8_gi2go = []
    skipped_n = 0
    processed_n = 0
    OUT.write("#read_id\tsubject_id\tgi\tGO\n")
    for l in IN:
        l = l.split("\t")
        gi = int(l[1].split("|")[1])
        go=""
        
        try:
            go = gi2go_map[gi]
        except:
            skipped_n += 1
        processed_n += 1
        m8_gi2go.append([l[0], l[1], gi, go])
        OUT.write("\t".join(map(str, [l[0], l[1], gi, go])) + "\n")

    OUT.close()
    IN.close()
        
    print("Processed=" + str(processed_n) + " (" + str(skipped_n) + " skipped)")
    
    return m8_gi2go



## Assume gi2go_fn does not have header line
def merge_gi2go_map(map_fn, gi2go_fn, out_fn):
    IN_map = open(map_fn)
    IN_gi2go = open(gi2go_fn)
    OUT = open(out_fn, "w")
    
    print("Processing " + map_fn + " and " + gi2go_fn)
    
    # Wash the comment line from map_fn
    flush = IN_map.readline()
    OUT.write("\t".join(["#", "read_id", "subject_id", "gi", "go", "tax_id", "seed_id"]) + "\n")
    for map_l, gi_l in izip(IN_map, IN_gi2go):
         map_l = map_l.strip().split("\t")
         gi_l = gi_l.strip().split("\t")
         
         if map_l[1] != gi_l[2]:
             print("Map and Gi2GO files are not synchronized, abort.")
             break
         
         OUT.write("\t".join(gi_l + map_l[2:4]) + "\n")
    
    IN_map.close()
    IN_gi2go.close()
    OUT.close()
    
         
         
         
    
    


"""
hmm_orf_dict = mg_pipeline.postprocess_HMMER_search_by_fn("all_samples.renamed+Pfam.dom.tbl", hmm_score_threshold=10, hmm_tc_fn="/home/siukinng/db/Markers/Pfam/Pfam.tc")

hmm_dom_tbl = mg_pipeline.generate_dom_tbl(hmm_orf_dict, hmm_accession_as_key=False)


"""