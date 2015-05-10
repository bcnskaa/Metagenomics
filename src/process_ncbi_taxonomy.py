from __future__ import print_function
from __future__ import division

import sys
import os
from threading import Thread

from ete2 import Tree
    

"""
nodes.dmp
---------

This file represents taxonomy nodes. The description for each node includes 
the following fields:

     tax_id                            -- node id in GenBank taxonomy database
     parent tax_id                     -- parent node id in GenBank taxonomy database
     rank                              -- rank of this node (superkingdom, kingdom, ...) 
     embl code                         -- locus-name prefix; not unique
     division id                       -- see division.dmp file
     inherited div flag  (1 or 0)      -- 1 if node inherits division from parent
     genetic code id                   -- see gencode.dmp file
     inherited GC  flag  (1 or 0)      -- 1 if node inherits genetic code from parent
     mitochondrial genetic code id     -- see gencode.dmp file
     inherited MGC flag  (1 or 0)      -- 1 if node inherits mitochondrial gencode from parent
     GenBank hidden flag (1 or 0)      -- 1 if name is suppressed in GenBank entry lineage
     hidden subtree root flag (1 or 0)  -- 1 if this subtree has no sequence data yet
     comments                           -- free-text comments and citations

names.dmp
---------
Taxonomy names file has these fields:

    tax_id          -- the id of node associated with this name
    name_txt        -- name itself
    unique name     -- the unique variant of this name if name not unique
    name class      -- (synonym, common name, ...)

division.dmp
------------
Divisions file has these fields:
    division id                -- taxonomy database division id
    division cde                -- GenBank division code (three characters)
    division name                -- e.g. BCT, PLN, VRT, MAM, PRI...
    comments

gencode.dmp
-----------
Genetic codes file:

    genetic code id                -- GenBank genetic code id
    abbreviation                -- genetic code name abbreviation
    name                    -- genetic code name
    cde                    -- translation table for this genetic code
    starts                    -- start codons for this genetic code

delnodes.dmp
------------
Deleted nodes (nodes that existed but were deleted) file field:

    tax_id                    -- deleted node id

merged.dmp
----------
Merged nodes file fields:

    old_tax_id                              -- id of nodes which has been merged
    new_tax_id                              -- id of nodes which is result of merging

citations.dmp
-------------
Citations file fields:

    cit_id                    -- the unique id of citation
    cit_key                    -- citation key
    pubmed_id                -- unique id in PubMed database (0 if not in PubMed)
    medline_id                -- unique id in MedLine database (0 if not in MedLine)
    url                    -- URL associated with citation
    text                    -- any text (usually article name and authors)
                        -- The following characters are escaped in this text by a backslash:
                        -- newline (appear as "\n"),
                        -- tab character ("\t"),
                        -- double quotes ('\"'),
                        -- backslash character ("\\").
    taxid_list                -- list of node ids separated by a single space

"""
map_fn_path="/home/siukinng/db/Markers/ncbi_nr/mapping/ncbi_tax"
prot_gi2taxid_map = {}



"""
"""
def assign_tax_family_to_m8_results(m8_fn, lineage_ofn=None, prot_gi2taxid_map=None):
    #global prot_gi2taxid_map
    
    with open(m8_fn) as IN:
        m8 = IN.read().splitlines()
    m8 = {m.split("\t")[0]: [int(m.split("\t")[1].split("|")[1]), m] for m in m8}
    
    skipped_n = 0
    processed_n = 0
    for id in m8.keys():
        gi = m8[id][0]
        taxid = ""
        try:
            taxid = prot_gi2taxid_map[gi]
            processed_n += 1
        except:
            skipped_n += 1
        m8[id].append(taxid)

    nr_taxids = list(set([m[2] for k,m in m8.items()]))

    [tax_nodes, tax_names] = generate_gi_tax_lineage_map()
    tax_lineages = get_tax_lineage(nr_taxids, tax_nodes, prot_gi2taxid_map, tax_names)
    
    print("Assigned taxid: Processed=" + str(processed_n) + ", skipped=" + str(skipped_n))
  
    
    skipped_n = 0
    processed_n = 0
   
    # Assign tax_lineage
    for id in m8.keys():
        taxid = m8[id][2] 
        tax_lineage = ""
        try:
            tax_lineage = tax_lineages[taxid]
            processed_n += 1
        except:
            skipped_n += 1
        m8[id].append(tax_lineage)
        
    print("Assigned Tax Lineage: Processed=" + str(processed_n) + ", skipped=" + str(skipped_n))

    m8_tax_lineage = m8
    
    [finalized_scafolds, scaffolds] = infer_scaffold_tax_lineage(m8_tax_lineage)
    
    if lineage_ofn is None:
        lineage_ofn = m8_fn + ".lineage"
        
    print("Exporting lineage results to " + lineage_ofn)
    with open(lineage_ofn, "w") as OUT:
        OUT.write("scaffold_id\tlineage\n")
        for scaffold_id in finalized_scafolds.keys():
            OUT.write(scaffold_id + "\t" + finalized_scafolds[scaffold_id] + "\n")
    print(str(len(finalized_scafolds)) + " lines wrote.") 
    return [finalized_scafolds, scaffolds, m8_tax_lineage]



"""
m8_tax_lineage: keys are derived from ORFs predicted from scaffolds by Prodigal 
key format: scaffold_XX_X
"""
#def infer_scaffold_tax_lineage(m8_tax_lineage, max_consistent_lineage=1, min_consistent_count=2, cutoff_perc=0.8):
def infer_scaffold_tax_lineage(m8_tax_lineage, max_consistent_lineage=5, cutoff_perc=0.8):
    import operator
    
    scaffolds = {}
    for id in m8_tax_lineage.keys():
        scaffold_id = id[::-1].split("_", 1)[1][::-1]
        if scaffold_id not in scaffolds.keys():
            #print("Adding " + scaffold_id)
            scaffolds[scaffold_id] = []
        lineage = m8_tax_lineage[id][3]
        lineage_family = ""
        lineage_order = ""
        if len(lineage) > 0:
            lineage_order = lineage.split(";")[3].replace("o__", "")
            lineage_family = lineage.split(";")[4].replace("f__", "")
        #scaffolds[scaffold_id].append(lineage_family)
        scaffolds[scaffold_id].append(lineage_order)
    
    finalized_scafolds = {}
    for id in scaffolds.keys():
        # Count the number of unique lineage
        total_lineages = len([lineage for lineage in scaffolds[id] if len(lineage) > 0])
        nr_lineages = [lineage for lineage in list(set(scaffolds[id])) if len(lineage) > 0]
        finalized_scafolds[id] = ""
        
        if len(nr_lineages) > max_consistent_lineage or len(nr_lineages) == 0:
            continue
    
        lineage_counts = {lineage: scaffolds[id].count(lineage) for lineage in nr_lineages}
        
        #if sorted(lineage_counts.items(), key=operator.itemgetter(1), reverse=True)[0][1] < min_consistent_count:
        highest_lineage_count = sorted(lineage_counts.items(), key=operator.itemgetter(1), reverse=True)[0][1]
        highest_lineage = sorted(lineage_counts.items(), key=operator.itemgetter(1), reverse=True)[0][0]
        highest_lineage_perc = highest_lineage_count / total_lineages
        print("Highest_lineage_count=" + str(highest_lineage_count) + " or " + str(highest_lineage_perc) +" (" + highest_lineage + ")")
        #if sorted(lineage_counts.items(), key=operator.itemgetter(1), reverse=True)[0][1] / total_lineages < cutoff_perc:
        
        if highest_lineage_perc < cutoff_perc:
            continue
        
        finalized_scafolds[id] = highest_lineage
    
    return [finalized_scafolds, scaffolds]
            


"""
Map the protein GI to taxid
"""
def get_taxid_from_protgi(prot_gi, gi2taxid_map=prot_gi2taxid_map):
    if len(gi2taxid_map) == 0:
        gi2taxid_map = import_prot_gi2taxid_map()
    
    print("Searching taxid for " + str(prot_gi))
    
    taxid = 0
    try:
        taxid = gi2taxid_map[int(prot_gi)]
    except Exception as e:
        print(taxid + " not existed.")
    #if prot_gi in prot_gi2taxid_map.keys():
    #    return prot_gi2taxid_map[prot_gi]
    #else:
    #    return None
    return taxid


"""
gi_taxid_prot_dmp_fn contains prot_gi and their mapping to taxids
"""
def import_prot_gi2taxid_map(gi_taxid_prot_dmp_fn=map_fn_path+"/gi_taxid_prot.dmp", force_reload=False):
    global prot_gi2taxid_map
    
    if len(prot_gi2taxid_map) > 0 or force_reload:
        print("gi2taxid_map is ready")
        return prot_gi2taxid_map
    
    print("Loading gi to taxid map table.") 
    print("Extracting data from " + gi_taxid_prot_dmp_fn)
    
    with open(gi_taxid_prot_dmp_fn) as IN:
        prot_gi2taxid_map = IN.read().splitlines()
    
    print(str(len(prot_gi2taxid_map)) + " lines read.")
    #prot_gi2taxid_map = [n.split("\t") for n in prot_gi2taxid_map]
    #prot_gi2taxid_map = {n[0]:n for n in prot_gi2taxid_map}
    
    #prot_gi2taxid_map = {n.split("\t")[0]:n.split("\t")[1] for n in prot_gi2taxid_map}
    #prot_gi2taxid_map = {int(n.split("\t")[0]):int(n.split("\t")[1]) for n in prot_gi2taxid_map}
    prot_gi2taxid_map2 = {}
    for n in prot_gi2taxid_map:
        n2 = n.split("\t")
        prot_gi2taxid_map2[int(n2[0])] = int(n2[1])

    prot_gi2taxid_map = prot_gi2taxid_map2
    
    print(str(len(prot_gi2taxid_map)) + " rows are extracted.")
    
    return prot_gi2taxid_map



"""
names.dmp contains taxids and their attributes
"""
def import_taxids(names_fn=map_fn_path+"/names.dmp"):
    print("Reading from " + names_fn)
    with open(names_fn) as IN:
        taxids = IN.read().splitlines()
    taxids = {int(n.split("\t|\t")[0]):n.split("\t|\t") for n in taxids}
    
    print(str(len(taxids)) + " names are read.")
    
    return taxids




def find_tax_lineage_for_taxid(taxid, nodes_fn=map_fn_path+"/nodes.dmp"):
    dist = 1 # Default distance
    taxid = str(taxid)
    
    print("Parsing " + nodes_fn + "...")
    with open(nodes_fn) as IN:
        nodes = IN.read().splitlines()
    nodes = {n.split("\t|\t", 1)[0]: n.replace("\t|\t", "\t").replace("|\t","").split("\t") for n in nodes}
    print("Searching id: " + taxid)
    if taxid in nodes.keys():
        print()
        
    else:
        print(taxid + " is not found.")
        return "Unknown"



"""

"""
def construct_tax_tree_from_nodes(nodes_fn=map_fn_path+"/nodes.dmp"):
    dist = 1 # Default distance
    
    print("Parsing " + nodes_fn + "...")
     
    with open(nodes_fn) as IN:
        nodes = IN.read().splitlines()
    nodes = {int(n.split("\t|\t", 1)[0]): n.replace("\t|\t", "\t").replace("|\t","").split("\t") for n in nodes}
    
    print(str(len(nodes)) + " nodes are imported.")
    
    
    # Estimate the number of nodes of each hierarchical level
    print("Estimating the number of nodes...")
    
    # We only select the following level tags
    # Superkingdom;phylum;class;order;family;genus;species
    levels = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    
    #levels = list(set([nodes[taxid][1] for taxid in nodes.keys()]))
    level_ns = {l:0 for l in levels}
    for level in level_ns.keys():
        level_ns[level] = len([1 for taxid in nodes.keys() if nodes[taxid][2] == level])

    print(", ".join(["=".join(map(str,m)) for m in level_ns.items()]))
   
    # Extract the species nodes, and tree construction starts from species level nodes
    species_nodes = [nodes[taxid] for taxid in nodes.keys() if nodes[taxid][2] == "species"]
    print(str(len(species_nodes)) + " species nodes.")   
    
    
    # Construct the tree
    print("Constructing tax tree...")
    tax_tree = Tree()
    tax_tree.name = "root"
    
    node_map = {}
    processed_n = 0
    singleton_n = 0
    for level in levels:
        print("Exploring level:" + level)
        level_processed_n = 0
        level_singleton_n = 0
        
        ns = [nodes[taxid] for taxid in nodes.keys() if nodes[taxid][2] == level]
        
        if level == "superkingdom":
            for n in ns:
                n[1] = "root"
        
        for n in ns:
            taxid = n[0]
            parent_taxid = n[1]
            
            ###
            if len(tax_tree.search_nodes(name=parent_taxid)) == 1:
                (tax_tree&parent_taxid).add_child(name=taxid, dist=dist)
            else:
                 singleton_n += 1
                 level_singleton_n += 1
            level_processed_n += 1
            processed_n += 1
        print("Number of nodes processed ("+level+"): " + str(level_processed_n))      
        print("Number of singleton ("+level +"): " + str(level_singleton_n))  
        
    print("Number of nodes processed (All levels): " + str(processed_n))      
    print("Number of singleton (All levels): " + str(singleton_n))
    
    # Create 
    #tax_tree = Tree()
    
    return tax_tree



"""

"""

#def generate_gi_tax_lineage_map(map_ofn, gi_taxid_prot_dmp_fn=map_fn_path+"/gi_taxid_prot.dmp", nodes_fn=map_fn_path+"/nodes.dmp", names_fn=map_fn_path+"/names.dmp"):
def generate_gi_tax_lineage_map(gi_taxid_prot_dmp_fn=map_fn_path+"/gi_taxid_prot.dmp", nodes_fn=map_fn_path+"/nodes.dmp", names_fn=map_fn_path+"/names.dmp"):
    print("Reading tax names from " + names_fn)
    with open(names_fn) as IN:
        tax_names = IN.read().splitlines()
    tax_names = {int(n.split("\t|\t")[0]):n.split("\t|\t")[0:2] for n in tax_names if n.endswith("scientific name\t|")}
    print(str(len(tax_names)) + " tax names read.")
    
    
    print("Reading tax nodes from " + nodes_fn)
    with open(nodes_fn) as IN:
        tax_nodes = IN.read().splitlines()
    tax_nodes = {int(n.split("\t|\t")[0]): [int(n.split("\t|\t")[1]), n.split("\t|\t")[2]] for n in tax_nodes}
    print(str(len(tax_nodes)) + " tax nodes imported.")
    
    
#     print("Reading protein GIs and their taxids from " + gi_taxid_prot_dmp_fn)
#     with open(gi_taxid_prot_dmp_fn) as IN:
#         gi2taxids = IN.read().splitlines()
#     gi2taxids = {int(n.split("\t")[0]):int(n.split("\t")[1]) for n in gi2taxids}
#     print(str(len(gi2taxids)) + " gi-taxid pairs extracted.")
    
#     print("Estimating number of non-redundant taxids")
#     nr_taxids = list(set([gi2taxids[gi] for gi in gi2taxids]))
#     print(str(len(nr_taxids)) + " nr taxids found.")
#     
#     print("Processing taxid")
#     tax_lineages = get_tax_lineage(nr_taxids, tax_nodes, gi2taxids, tax_names)
#     print("Number of lineages: " + str(len(tax_lineages)))
#     
#     print("Exporting lineage map to " + map_ofn)
#     OUT = open(map_ofn, "w")
#     for gi, taxid in gi2taxids.items():
#         OUT.write("\t".join([str(gi), str(taxid), tax_lineages[taxid]]) + "\n")
#     OUT.close()
    
    
#     print("Processing taxid")
#     for taxid in nr_taxids[0:10]:
#         tax_lineage = []
#         
#         if taxid in tax_nodes.keys():
#             tax_name = ""
#             if taxid in tax_names.key():
#                 tax_name = tax_names[taxid][1]
#     
#             parent_taxid = tax_nodes[taxid][0]        
#             parent_tax_name = ""
#             if parent_taxid in tax_names.key():
#                 parent_tax_name = tax_names[parent_taxid][1]
#                 
#             node_type = tax_nodes[taxid][1]
#             
#             print(str(taxid) + ": tax_name=" + tax_name + " parent=" + parent_tax_name + "(taxid=" + str(parent_taxid) + ") node_type=" + node_type)
#         else:
#             print(taxid + " is not found in tax nodes.")
#     
    #return [tax_nodes, gi2taxids, tax_names, tax_lineages]
    return [tax_nodes, tax_names]




def get_tax_lineage_mt(tid, tax_lineages, tax_nodes, gi2taxids, tax_names):
    # Lineage levels
    tax_ranks = ["k__", "p__", "c__", "o__","f__", "g__", "s__"]
    
    process_n = 0
    skipped_n = 0
    #tax_lineages = {}
    
    taxids = list(tax_lineages.keys())
    print("(" + str(tid) + ") Number of taxids to be processed: " + str(len(taxids)))
#    print("Processing taxids")
    for taxid in taxids:
        process_n += 1
        
        if process_n % 10000 == 0:
            print(str(tid) + ": " + str(process_n) + " processed")
        
        # k__;p__;c__;o__;f__;g__;s__
        lineage = ["","","","","","",""]
        try:
            cur_taxid = taxid
            node_type = tax_nodes[cur_taxid][1]
            
            while(node_type != "no rank"):
                node_type = tax_nodes[cur_taxid][1]
                tax_name = tax_names[cur_taxid][1]
            
                if node_type == "superkingdom":
                    lineage[0] = tax_name    
                elif node_type == "phylum":
                    lineage[1] = tax_name    
                elif node_type == "class":
                    lineage[2] = tax_name    
                elif node_type == "order":
                    lineage[3] = tax_name    
                elif node_type == "family":
                    lineage[4] = tax_name    
                elif node_type == "genus":
                    lineage[5] = tax_name    
                elif node_type == "species":
                    lineage[6] = tax_name
                else:
                    0
                
                cur_taxid = tax_nodes[cur_taxid][0]
            tax = ";".join(z[0] + z[1] for z in zip(tax_ranks, lineage))
            tax_lineages[taxid] =  tax_ranks[0] + lineage[0] + ";" + tax_ranks[1] + lineage[1] + ";" + tax_ranks[2] + lineage[2] + ";" + tax_ranks[3] + lineage[3] + ";" + tax_ranks[4] + lineage[4] + ";" + tax_ranks[5] + lineage[5] + ";" + tax_ranks[6] + lineage[6]
                
        except:
            skipped_n += 1

    print("Number of processed tax ids: " + str(process_n) + " ("+ str(skipped_n) +")")      
        
    return tax_lineages



"""

Construct the lineage

"""
def get_tax_lineage(taxids, tax_nodes, gi2taxids, tax_names, num_threads=16):
    taxid_n = len(taxids)
    taxid_n_per_t = taxid_n // num_threads
    print("Number of taxids: " + str(taxid_n))
    print("Number of threads: " + str(num_threads))
    print("Number of taxids per thread: " + str(taxid_n_per_t))
    
    
    threads = [None] * num_threads
    #thread_results = [None] * num_threads
    thread_tax_lineages = dict(zip(taxids, ["" for x in range(taxid_n)]))
    thread_results = [ dict(thread_tax_lineages.items()[x:x+taxid_n_per_t]) for x in xrange(0, taxid_n, taxid_n_per_t)]
    
    for i in range(num_threads):
        # get_tax_lineage_mt(tid, tax_lineage, tax_nodes, gi2taxids, tax_names): 
        threads[i] = Thread(target=get_tax_lineage_mt, args=(i, thread_results[i], tax_nodes, gi2taxids, tax_names))
        threads[i].start()

    # do some other stuff
    for i in range(num_threads):
        threads[i].join()

    tax_lineages = {}
    for thread_result in thread_results:
        tax_lineages.update(thread_result)
    
    
    return tax_lineages


    
# def get_tax_lineage(taxids, tax_nodes, gi2taxids, tax_names):
#     # Lineage levels
#     tax_ranks = ["k__", "p__", "c__", "o__","f__", "g__", "s__"]
#     
#     process_n = 0
#     skipped_n = 0
#     tax_lineages = {}
#     
# #    print("Processing taxids")
#     for taxid in taxids:
#         # k__;p__;c__;o__;f__;g__;s__
#         lineage = ["","","","","","",""]
#         
#         if taxid in tax_nodes.keys():
# #             tax_name = ""
# #             if taxid in tax_names.keys():
# #                 tax_name = tax_names[taxid][1]
# #             node_type = tax_nodes[taxid][1]
#             
#             #tax_lineage[node_type] = tax_name
#             #tax_lineage.append(node_type+"="+tax_name)
#             
#             cur_taxid = taxid
#             node_type = tax_nodes[cur_taxid][1]
#             while(node_type != "no rank"):
#                 node_type = tax_nodes[cur_taxid][1]
#                 tax_name = ""
#                 if cur_taxid in tax_names.keys():
#                     tax_name = tax_names[cur_taxid][1]
# 
#                 if node_type == "superkingdom":
#                     lineage[0] = tax_name    
#                 elif node_type == "phylum":
#                     lineage[1] = tax_name    
#                 elif node_type == "class":
#                     lineage[2] = tax_name    
#                 elif node_type == "order":
#                     lineage[3] = tax_name    
#                 elif node_type == "family":
#                     lineage[4] = tax_name    
#                 elif node_type == "genus":
#                     lineage[5] = tax_name    
#                 elif node_type == "species":
#                     lineage[6] = tax_name
#                 else:
#                     0
#                 
#                 cur_taxid = tax_nodes[cur_taxid][0]      
#                 
#                 #tax_lineage[node_type] = tax_name
#                 #tax_lineage.append(node_type+"="+tax_name+"("+str(parent_taxid)+")")
#             #print(", ".join([type + "=" + tax_lineage[type] for type in tax_lineage.keys()]))
#             tax = ";".join(z[0] + z[1] for z in zip(tax_ranks, lineage))
#             tax_lineages[taxid] = tax
#             #print(str(taxid) + ": " + tax)
#             
#             #tax_lineage = tax_lineage[::-1]
#             #print(str(taxid) + ": " + ", ".join(tax_lineage))
#                   
# #             parent_tax_name = ""
# #             if parent_taxid in tax_names.keys():
# #                 parent_tax_name = tax_names[parent_taxid][1]
# #             parent_node_type = tax_nodes[parent_taxid][1]
# #             
# #             print(str(taxid) + ": tax_name=" + tax_name + ", node_type=" + node_type + ", parent=" + parent_tax_name + "(taxid=" + str(parent_taxid) + ", node_type=" + parent_node_type + ")")
#         else:
#             skipped_n += 1
#             #print(str(taxid) + " is not found in tax nodes.")    
#         
#         process_n += 1
#        
#     print("Number of processed tax ids: " + str(len(process_n)) + " ("+ str(len(skipped_n)) +")")      
#         
#     return tax_lineages



