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


"""
import glob
import process_m8

m8_fns = glob.glob("*.m8")
for m8_fn in m8_fns:
    process_m8.merge_gi2go_map(m8_fn+".map", m8_fn+".gi2go", m8_fn + ".gi2go+map", dump_read_id=True)

"""
## Assume gi2go_fn does not have header line
def merge_gi2go_map(map_fn, gi2go_fn, out_fn, dump_read_id=False):

    IN_map = open(map_fn)
    IN_gi2go = open(gi2go_fn)
    OUT = open(out_fn, "w")
    if dump_read_id:
        dump_out_fn = out_fn + ".id.dump"
        OUT_dump = open(dump_read_id, "w")
    
    # Wash the comment line from map_fn
    flush = IN_map.readline()
   
    print("Processing " + map_fn + " and " + gi2go_fn)
    
    if dump_read_id:
        OUT.write("\t".join(["#line_id", "gi", "go", "tax_id", "seed_id"]) + "\n")
        OUT_dump.write("\t".join(["#line_id", "read_id", "subject_id"]) + "\n")
        line_id = 0
        for map_l, gi_l in izip(IN_map, IN_gi2go):
             map_l = map_l.rstrip("\n").split("\t")
             gi_l = gi_l.rstrip("\n").split("\t")
             
             if map_l[1] != gi_l[2]:
                 print("Map and Gi2GO files are not synchronized, abort.")
                 break
             line_id += 1
             OUT_dump.write(str(line_id) + "\t" + gi_l[1:2] + "\n")
             OUT.write(str(line_id)+ "\t" + "\t".join(gi_l[2:] + map_l[2:4]) + "\n")
    else:   
        OUT.write("\t".join(["#read_id", "subject_id", "gi", "go", "tax_id", "seed_id"]) + "\n")
        for map_l, gi_l in izip(IN_map, IN_gi2go):
             map_l = map_l.rstrip("\n").split("\t")
             gi_l = gi_l.rstrip("\n").split("\t")
             
             if map_l[1] != gi_l[2]:
                 print("Map and Gi2GO files are not synchronized, abort.")
                 break
             
             OUT.write("\t".join(gi_l + map_l[2:4]) + "\n")
    
    IN_map.close()
    IN_gi2go.close()
    OUT.close()
    OUT_dump.close()
    


"""
Given the gi2go_fn, this function will extract the specified functional trait (GO, SEED, PFAM, COG)

Usage:
import process_m8
import glob
import os

gi2go_fns = glob.glob("*.gi2go+map")
sample_ids = [(os.path.basename(f)).split(".")[0] for f in gi2go_fns]
[sample_list, trait_list, sample_ids, tax_ids, trait_ids] = process_m8.convert_gi2go_to_FDcompatible(gi2go_fns, tax_col_id="tax_id", trait_col_id="go", sample_ids=sample_ids)

"""
def convert_gi2go_to_FDcompatible(gi2go_fns, tax_col_id, trait_col_id, taxids_clusters=None, excluding_zero=False, tax_val_delim=None, sample_ids=None, ofn_prefix=None, sample_ofn=None, trait_ofn=None, sample_trait_ofn=None, trait_val_delim="; ", missing_value="NA"):
    import os
    
    if sample_ids is None:
        sample_ids = {os.path.splitext((os.path.basename(fn)))[0]:fn for fn in gi2go_fns}
    
    print(str(len(sample_ids)) + " ids found:")
    print("-----------------------------------\n" + "\n".join(sample_ids) + "\n")
    
    if ofn_prefix is None:
        ofn_prefix = "samples"
        if taxids_clusters is not None:
            ofn_prefix = ofn_prefix + "-clustered"
        
    if sample_ofn is None:
        sample_ofn = ofn_prefix + "-" + tax_col_id + "+" + trait_col_id + ".samples"
    if sample_trait_ofn is None:
        sample_trait_ofn = ofn_prefix + "-" + tax_col_id + "+" + trait_col_id + ".samples.trait"
    if trait_ofn is None:
        trait_ofn = ofn_prefix + "-" + tax_col_id + "+" + trait_col_id + ".traits"
            
    # Sample_list format:
    # RIDS:CID       tax_id_1    tax_id_2      ...
    # sample_id_1    1            0            ...
    # sample_id_2    3            12           ...
    # ...            ...          ...          ...
    sample_list = {}
    
    # Tax_list format:
    # RID:CID     trait_id_1    trait_id_2     ...
    # tax_id_1    10            21             ...
    # tax_id_2    1             400            ...
    # ...         ...           ...            ...
    trait_list = {}
    # Trait list for each sample
    sample_trait_list = {}   
    
    
    tax_ids = []
    trait_ids = []

    # Summary of clusters
    cluster_stats = {}
    if taxids_clusters is not None:
        for tid in taxids_clusters.keys():
            cluster_id = taxids_clusters[tid]
            try:
                cluster_stats[cluster_id][tid] = 0
            except:
                cluster_stats[cluster_id] = {}
                cluster_stats[cluster_id][tid] = 0
                
    
    # Iterate all the samples
    for sample_id in sample_ids.keys():
        print("Processing " + sample_id)
        
        gi2go_fn = sample_ids[sample_id]
        
        # Read the infile
        IN = open(gi2go_fn)
        header = IN.readline()
        
        # Remove comment and split column ids
        if not header.startswith("#"):
            print("Either header line is missing or does not start with '#', " + gi2go_fn + ",  skipped.")
            continue
        
        header = header.replace("#", "").split("\t")
        
        # Search if tax_col_id exists
        try:
            tax_col_idx = header.index(tax_col_id)
        except:
            print("Tax column ID not found " + gi2go_fn + ": " + tax_col_id + ", file skipped.")
            continue
        
        # Search if trait_col_id exists
        try:
            trait_col_idx = header.index(trait_col_id)
        except:
            print("Trait column ID not found in " + gi2go_fn + ": " + trait_col_id + ", file skipped.")
            continue
        
        # Initialize lists for the current sample_id
        sample_list[sample_id] = {}
        sample_trait_list[sample_id] = {}
        
        skipped_n = 0
        processed_n = 0
        # Proces the content
        for line in IN:
            #if line.startswith("#"):
            #    continue
            line = line.split("\t")
            
            tax_val = int(line[tax_col_idx])
            trait_vals = line[trait_col_idx].split(trait_val_delim)
            
            # Reassign the taxids according to specified cluster
            if taxids_clusters is not None:
                try:
                    cluster_id = taxids_clusters[tax_val]
                    cluster_stats[cluster_id][tax_val] += 1
                except:
                    skipped_n += 1
                    continue
                # Reassign taxid to cluster_id
                tax_val = cluster_id
            
            processed_n += 1
            
            # Skip if nothing presents in the tax value column or trait column
            #if len(tax_val) == 0 or len(trait_vals) == 0 or (len(trait_vals) == 1 and len(trait_vals[0]) == 0):
            #if len(tax_val) == 0 or len(trait_vals) == 0:
            if len(trait_vals) == 0:
                skipped_n += 1
                continue
            
            # Adding the tax val
            try:
                sample_list[sample_id][tax_val] += 1
            except:
                sample_list[sample_id][tax_val] = 1
            
            # Adding traits to trait tables
            for trait_val in trait_vals:
                # Handling the case that the trait value is empty
                if len(trait_val) == 0:
                    trait_val = missing_value
            
                # Adding the trait to the sample trait table    
                try:
                    sample_trait_list[sample_id][trait_val] += 1
                except:
                    sample_trait_list[sample_id][trait_val] = 1
            
                # Adding the trait to the general trait table
                try:
                    trait_list[tax_val][trait_val] += 1
                except:
                    try:
                        trait_list[tax_val][trait_val] = 0
                    except:
                        trait_list[tax_val] = {}
                    trait_list[tax_val][trait_val] = 1
        
        print("Sample " + sample_id + ": " + str(processed_n) + " (" + str(skipped_n) + " skipped)")
        
        # Append tax and trait ids to lists
        tax_ids.extend(sample_list[sample_id].keys())
        trait_ids.extend(sample_trait_list[sample_id].keys())
        
        IN.close()
        
    # Non-redundant list
    trait_ids = list(set(trait_ids))
    tax_ids = list(set(tax_ids))
    
    print("Total number of traits: " + str(len(trait_ids)))
    print("Total number of taxa: " + str(len(tax_ids)))
    
    
    
    # Filter cluster by their counts
    # Obtain totals of each cluster
    
    
    ###########################  
    # Export the sample table
    ###########################
    print("Exporting sample table to " + sample_ofn)
    except_n = 0
    with open(sample_ofn, "w") as OUT:
        OUT.write("\t".join(["#sample_id:"+tax_col_id] + tax_ids) + "\n")
        for sample_id in sample_ids:
            tax_vals = [0 for i in tax_ids]
            for i, tax_id in enumerate(tax_ids):
                processed_n += 1
                try:
                    #tax_vals[i] = len(sample_list[sample_id][tax_id])
                    tax_vals[i] = sample_list[sample_id][tax_id]
                except:
                    except_n += 1
            OUT.write("\t".join([sample_id] + map(str, tax_vals)) + "\n")
    print("Processed: " + str(processed_n) + ", skipped: " + str(except_n))
    
    
    ###########################
    # Export the sample trait table
    ###########################  
    print("Exporting sample specific trait table to " + sample_trait_ofn)
    except_n = 0
    processed_n = 0
    with open(sample_trait_ofn, "w") as OUT:
        OUT.write("\t".join(["#sample_id", trait_col_id, "count"]) + "\n")
        for sample_id in sample_ids:
            for trait_id in trait_ids:
                processed_n += 1
                trait_n = 0
                try:
                    #trait_n = len(trait_list[sample_id][trait_val])
                    trait_n = sample_trait_list[sample_id][trait_id]
                except:
                    except_n += 1
                 
                if excluding_zero and trait_n == 0:
                    continue
                 
                OUT.write("\t".join([sample_id, trait_id, str(trait_n)]) + "\n")
    print("Processed: " + str(processed_n) + ", skipped: " + str(except_n))
     
     
    ###########################
    # Export the general trait table
    ###########################  
    print("Exporting tax specific trait table to " + trait_ofn)
    except_n = 0
    processed_n = 0
    with open(trait_ofn, "w") as OUT:
        OUT.write("\t".join(["#tax_id:" + trait_col_id] + trait_ids) + "\n")
        for tax_id in tax_ids:
            trait_ns = [0 for t in trait_ids]
            for i, trait_id in enumerate(trait_ids):
                processed_n += 1
                 
                try:
                    trait_ns[i] = trait_list[tax_id][trait_id]
                except:
                    except_n += 1
                 
            OUT.write("\t".join([tax_id] + map(str, trait_ns)) + "\n")
                 
    print("Processed: " + str(processed_n) + ", skipped: " + str(except_n))    
    
#    return [sample_list, sample_trait_list, trait_list, sample_ids, tax_ids, trait_ids]
    return [sample_list, sample_trait_list, sample_ids, tax_ids, trait_ids, cluster_stats]


    

# def add_pfam_to_gi2go(gi2go_fn, go2pfam_map=None, out_fn=None, dump_read_id=True):
#     if out_fn is None:
#         out_fn = gi2go_fn + "+pfam"
#     
#     if go2pfam_map is None:
#         go2pfam_map = import_pfam2go()
#     
#     IN_gi2go = open(gi2go_fn)
#     OUT = open(out_fn, "w")
#     
#     print("Appending Pfam data to " + gi2go_fn + " and exporting to " + out_fn)
#     
#     # Wash the comment line from map_fn
#     nomatch_n = 0
#     processed_n = 0
#     
#     OUT.write("\t".join(["#read_id", "subject_id", "gi", "go", "tax_id", "seed_id", "pfam"]) + "\n")
#     for gi_l in IN_gi2go:
#         gi_l = gi_l.split("\t")
#         gos = gi_l[3].split("; ")
#         
#         pfam = []
#         for go in gos:
#             try:
#                 pfam.extend(list(go2pfam_map[go]["Pfam"].keys()))
#                 processed_n += 1
#             except:
#                 nomatch_n += 1
#                 
#         OUT.write("\t".join(gi_l + [";".join(pfam)]) + "\n")
# 
#     print("Processed: " + str(processed_n) + ", no match="+str(nomatch_n))
#     
#     IN_gi2go.close()
#     OUT.close()
   
           
         
def import_pfam2go(pfam2go_fn="/home/siukinng/db/Markers/GeneOntology/pfam2go"):
    print("Importing from " + pfam2go_fn)
    with open(pfam2go_fn) as IN:
         pfam2go_map = IN.read().splitlines()
    # Wash header lines
    pfam2go_map = [(p.split(" > ")[0]).split(" ") + (p.split(" > ")[1]).split(" ; ") for p in pfam2go_map if not p.startswith("!")]
    
    # Construct the go2pfam_map
    go2pfam_map = {}
    duplicated_n = 0
    processed_n = 0
    for p in pfam2go_map:
        go = p[3]
        pfam = p[0].replace("Pfam:", "")
        pfam_name = p[1]
        go_ctx = p[2].replace("GO:", "")
        try:
            go2pfam_map[go]["Pfam"][pfam] = pfam_name
            #go2pfam_map[go]["Pfam_Output"].append(pfam)
        except:
            go2pfam_map[go] = {}
            go2pfam_map[go]["GO"] = go_ctx
            go2pfam_map[go]["Pfam"] = {}
            go2pfam_map[go]["Pfam"][pfam] = pfam_name
            #go2pfam_map[go]["Pfam_Output"] = [pfam]
            
            
        processed_n += 1
        
    print("Number of GO processed: " + str(processed_n) + " (" + str(duplicated_n) + " duplicated)")
    
    return go2pfam_map
        
   
   
"""
import process_m8
import glob
import os

gi2go_fns = glob.glob("*.gi2go+map")

[sample_list, trait_list, sample_ids, tax_ids, trait_ids] = process_m8.convert_gi2go_to_FDcompatible(gi2go_fns, "tax_id", "go")
[cluster_ids, cluster2taxid, taxid2cluster] = process_m8.cluster_taxids_by_tax_lineage(tax_ids)
[sample_list, trait_list, sample_ids, tax_ids, trait_ids, cluster_stats] = process_m8.convert_gi2go_to_FDcompatible(gi2go_fns, "tax_id", "go", taxids_clusters=taxid2cluster)

"""
def cluster_taxids_by_tax_lineage(taxids, taxid2lineage_map=None, tax_level="g", excluding_missing_value=False, missing_value_label="Unassigned"):
    print("Processing " + str(len(taxids)) + " taxid using tax_level: " + tax_level)
    
    # Convert taxids into integer for efficient processing
    taxids = map(int, taxids)
    
    tax_level_tag = tax_level + "__"
    
    # If taxids sharing a same tax_level are clustered
    taxid2cluster = {}
    cluster2taxid = {}
    cluster_ids = []
    if taxid2lineage_map is None:
        taxid2lineage_map = import_taxid2lineage()
    
    
    print("Clustering")
    skipped_n = 0
    processed_n = 0
    for taxid in taxids:
        processed_n += 1
        try:
            lineage = taxid2lineage_map[taxid]
        except:
            skipped_n += 1
            continue
        
        # Extract the value at specified tax level
        tax_level_val = [l for l in lineage.split(";") if l.startswith(tax_level_tag)]
        
        # tax value goes empty, so marked "Unassigned"
        if len(tax_level_val) == 0:
            tax_level_val = missing_value_label
        
        tax_level_val = tax_level_val[0]
        taxid2cluster[taxid] = tax_level_val
        try:
            cluster2taxid[tax_level_val].append(taxid)
        except:
            cluster2taxid[tax_level_val] = []
            cluster2taxid[tax_level_val].append(taxid)
            cluster_ids.append(tax_level_val)
    
    
    
    print("Processed: " + str(processed_n) + ", skipped: " + str(skipped_n))
    print("Number of cluster: " + str(len(cluster2taxid)))
     
    return [cluster_ids, cluster2taxid, taxid2cluster]    
    
    
 


def import_taxid2lineage(taxid2lineage_fn="/home/siukinng/db/Markers/ncbi_nr/mapping/taxid2tax_lineage.map"):
    print("Importing taxid to lineage map from " + taxid2lineage_fn)
    with open(taxid2lineage_fn) as IN:
        taxid2lineage_map = IN.read().splitlines()
    
    if taxid2lineage_map[0].startswith("#"):
        del taxid2lineage_map[0]
    
    taxid2lineage_map = [t.split("\t") for t in taxid2lineage_map]
    taxid2lineage_map = {int(t[0]):t[1] for t in taxid2lineage_map}
    
    print(str(len(taxid2lineage_map)) + " row imported.")
    
    return taxid2lineage_map
        
"""
hmm_orf_dict = mg_pipeline.postprocess_HMMER_search_by_fn("all_samples.renamed+Pfam.dom.tbl", hmm_score_threshold=10, hmm_tc_fn="/home/siukinng/db/Markers/Pfam/Pfam.tc")

hmm_dom_tbl = mg_pipeline.generate_dom_tbl(hmm_orf_dict, hmm_accession_as_key=False)


"""