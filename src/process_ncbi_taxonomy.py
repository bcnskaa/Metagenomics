import sys
import os

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


"""
gi_taxid_prot_dmp_fn contains prot_gi and their mapping to taxid
"""
def import_prot_gi2taxid_map(gi_taxid_prot_dmp_fn):
    print("Extracting data from " + gi_taxid_prot_dmp_fn)
    
    with open(names_fn) as IN:
        prot_gi2taxid = IN.read().splitlines()
    prot_gi2taxid = {n.split("\\t")[0]:n.split("\t")[1] for n in prot_gi2taxid}
    
    print(str(len(prot_gi2taxid)) + " rows are extracted.")
    
    return prot_gi2taxid



"""
names.dmp contains taxids and their attributes
"""
def import_taxids(names_fn="names.dmp"):
    print("Reading from " + names_fn)
    with open(names_fn) as IN:
        taxids = IN.read().splitlines()
    taxids = {n.split("\t|\t")[0]:n.split("\t|") for n in taxids}
    
    print(str(len(taxids)) + " names are read.")
    
    return taxids


"""

"""
def construct_tax_tree_from_nodes(nodes_fn="nodes.dmp"):
    dist = 1 # Default distance
    
    print("Parsing " + nodes_fn + "...")
     
    with open(nodes_fn) as IN:
        nodes = IN.read().splitlines()
    nodes = {n.split("\t|\t", 1)[0]: n.replace("\t|\t", "\t").replace("|\t","").split("\t") for n in nodes}
    
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
    
    



