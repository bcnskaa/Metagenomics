import networkx
import numpy
#http://networkx.github.io/documentation/latest/overview.html



class DistMtx:
    # Instance variable
    #self.distmtx
    #self.header
    #self.file_source
    #self.mtx_pairs
    
    #def __init__(self, dist_mtx_fn, header=True, is_float=True, inversing_distance=False):
    def __init__(self, dist_mtx_fn, header=True, inversing_distance=False):
        self.distmtx = [[]]
        self.header = []
        self.file_source = dist_mtx_fn
        self.mtx_pairs = {}
        self.inversing_distance = inversing_distance 
        
        offset = 0.0
        if inversing_distance:
            offset = 1.0
        
        print("Paring distance matrix from " + self.file_source)
        
        with open(dist_mtx_fn) as IN:
            raw_mtx = IN.read().splitlines()
        self.header = raw_mtx[0].split("\t")[1:]
        del raw_mtx[0]
        
        #if is_float:
        self.distmtx = [[0.0 for i in xrange(len(self.header))] for j in xrange(len(self.header))]
        #else:
        #self.distmtx = [[0 for i in xrange(len(self.header))] for j in xrange(len(self.header))]
        
        #raw_mtx = {m.split("\t")[0]:m.split("\t")[1:] for m in raw_mtx}  
        for row in raw_mtx:
            items = row.split("\t")
            id = items[0]
            #if is_float:
            vals = [abs(offset - float(v)) for v in items[1:]]
            #else:
            #   vals = [abs(offset - int(v)) for v in items[1:]]
                
            idx = self.header.index(id)
            self.distmtx[idx] = vals
            
            for i, item in enumerate(items[1:]):
                row_id = id
                col_id = self.header[i]
                
                #if is_float:
                self.mtx_pairs[(row_id, col_id)] = abs(offset - float(item))
                #else:
                #    self.mtx_pairs[(row_id, col_id)] = int(item)



    """
    return the element at the position i, j
    """     
    def get(self, i, j):
        if i > 0 and j > 0 and i < len(self.header) and j < len(self.header):
            return self.distmtx[i][j]
        else:
            print("Warning: " + str(i) + ", " + str(j) + " not found.")
            return None



    def get_by_ids(self, row_id, col_id):
        if row_id not in self.header or col_id not in self.header:
            print("Warning: " + row_id + ", " + col_id + " not found.")
            return None
        
        row_idx = self.header.index(row_id)



    def getIdx(self, id):
        if id in self.header:
            return self.header.index(id)
        else:
            return -1



    def size(self):
        return len(self.header)



    def getIDs(self):
        return self.header
    
    
    def get_pair_values(self, pair_keys):
        pair_values = []
        for pk in pair_keys:
            if pk in self.mtx_pairs.keys():
                pair_values.append(pk + (self.mtx_pairs[pk], ))
        return pair_values


    
    def get_pair_value(self, pair_key):
        if pair_key in self.mtx_pairs.keys():
            return self.mtx_pairs[pair_key]
        else:
            return None
    
    
    def get_matched_pair(self, pair_key):
        vals = []
        if pair_key in self.mtx_pairs.keys():
            vals.append(self.mtx_pairs[pair_key])
        else:
            return None
        
        pair_key = (pair_key[1], pair_key[0])
        if pair_key in self.mtx_pairs.keys():
            vals.append(self.mtx_pairs[pair_key])
        else:
            return None
    
        return vals
       
        
    """
    Get the averaged match pair
    """
    def get_average_matched_pair(self, pair_key):
        val = 0.0
        if pair_key in self.mtx_pairs.keys():
             val = val + self.mtx_pairs[pair_key]
        else:
            return None
        
        pair_key = (pair_key[1], pair_key[0])
        if pair_key in self.mtx_pairs.keys():
             val = val + self.mtx_pairs[pair_key]
        else:
            return None
    
        return val / 2.0
    
    
    
    
    """
    Get the pair keys whose values are greater than the cutoff value
    """
    def extract_pair_keys(self, dist_cutoff=0.5, excluding_diagonal=False, is_greater_than=True):  
        if is_greater_than:
            pair_keys = [k for k in self.mtx_pairs.keys() if self.mtx_pairs[k] >= dist_cutoff]
        else:
            pair_keys = [k for k in self.mtx_pairs.keys() if self.mtx_pairs[k] <= dist_cutoff] 
        
        if excluding_diagonal:
            nr_pair_keys = []
            for pk in pair_keys:
                if pk[0] != pk[1]:
                    nr_pair_keys.append(pk)
                    
            return nr_pair_keys 
        else:
            return pair_keys
    
    
    
    def extract_pair_key_values(self, dist_cutoff=0.5, excluding_diagonal=False, is_greater_than=True):  
        #pair_keys = [k for k in self.mtx_pairs.keys() if self.mtx_pairs[k] >= dist_cutoff]
        if is_greater_than:
            pair_key_values = [k + (self.mtx_pairs[k], ) for k in self.mtx_pairs.keys() if self.mtx_pairs[k] >= dist_cutoff]
        else:
            pair_key_values = [k + (self.mtx_pairs[k], ) for k in self.mtx_pairs.keys() if self.mtx_pairs[k] <= dist_cutoff]
        
            
        if excluding_diagonal:
            nr_pair_key_values = []
            for pk in pair_key_values:
                if pk[0] != pk[1]:
                    nr_pair_key_values.append(pk)
                    
            pair_key_values = nr_pair_key_values
        
        return pair_key_values
    
        #pair_key_values = []
        #for k in pair_keys:
        #    pair_key_values.append(k + (self.mtx_pairs[k],))
        
                    
                
    def get_pair_keys(self):
        return self.mtx_pairs.keys()
 


    def get_nr_ids_from_pair_keys(self, pair_keys):
        nr_ids = [pk[0] for pk in pair_keys] + [pk[1] for pk in pair_keys]
        nr_ids = list(set(nr_ids))
        return nr_ids
    
    

def import_node_abund(abund_fn="all.updated.tax.abund.lst"):
    with open(abund_fn) as IN:
        abund = IN.read().splitlines()
    abund = {a.split("\t")[1] : float(a.split("\t")[3]) for a in abund}
    return abund
        
        
 
def process_graph(cutoff=0.8, excluding_diagonal=True, gml_outfn=None, node_factor=1000):
    import networkx as nx
    import matplotlib.pyplot as plt
    import re
    
    G = nx.DiGraph()
    
    dixt_mtx_fn = "dist_2000.mtx" 
    abund_fn="all.updated.tax.abund.lst"
    
    abund = import_node_abund(abund_fn = abund_fn)
    dm = DistMtx(dixt_mtx_fn, inversing_distance=True)
    
    links = dm.extract_pair_keys(cutoff, excluding_diagonal=excluding_diagonal, is_greater_than=False)
    weighted_links = dm.get_pair_values(links)
    linked_nodes = dm.get_nr_ids_from_pair_keys(links)
    #G.add_nodes_from(linked_nodes)
    for n in linked_nodes:
        G.add_node(n, Weight=int(node_factor*float(abund[n])))
        
    #G.add_weighted_edges_from(weighted_links)
    for l in weighted_links:
        G.add_edge(l[0], l[1], len=float(l[2]))
        #G.add_edge(l[0], l[1], len=float(l[2]), weight=float(l[2]))
    G.add_edges_from(links)
    
    # Ancestral node
    abund["GS0"] = 1.0
    links_to_GS0 = [("GS0", n, 0.01) for n in linked_nodes if re.match("^G.1", n)]
    G.add_node("GS0", Weight=int(node_factor*float(abund["GS0"])))
    #G.add_weighted_edges_from(links_to_GS0)    
    for l in links_to_GS0:
        G.add_edge(l[0], l[1], len=float(l[2]))
        #G.add_edge(l[0], l[1], len=float(l[2]), weight=float(l[2]))
    
    
    abund["SS0"] = 1.0
    links_to_SS0 = [("SS0", n, 0.01) for n in linked_nodes if re.match("^S.1", n)]
    G.add_node("SS0", Weight=int(node_factor*float(abund["SS0"])))
    #G.add_weighted_edges_from(links_to_SS0)
    for l in links_to_SS0:
        G.add_edge(l[0], l[1], len=float(l[2]))
        #G.add_edge(l[0], l[1], len=float(l[2]), weight=float(l[2]))
    
    
    all_links = links + links_to_GS0 + links_to_SS0
    
    
    #pos = nx.spring_layout(G) # positions for all nodes
    pos=nx.spring_layout(G,iterations=1000)

    # nodes
    #nx.draw_networkx_nodes(G,pos,node_size=80)
    nx.draw_networkx_nodes(G, pos, node_size=[float(abund[g])*node_factor for g in G.nodes()], alpha=0.7)

    # edges
    nx.draw_networkx_edges(G,pos,edgelist=all_links)

    # labels
    nx.draw_networkx_labels(G,pos,font_size=10,font_family='sans-serif')

    plt.axis('off')
    
    #nx.draw(G)
    plt.show()
    
    if gml_outfn is None:
        gml_outfn = "nodes_" + str(cutoff) + ".gml"
    print("Exporting to " + gml_outfn)
    nx.write_gml(G, gml_outfn)
    
    return G
    
    
# 
#http://stackoverflow.com/questions/1851296/specified-edge-lengths-on-networkx-igraph-python/1898456#1898456
def import_dist_mtx(dist_mtx_fn="dist_2000.mtx"):
    dist_mtx = DistMtx(dist_mtx_fn)
    
    
def test_code():
    import pylab as p
    import networkx as nx

    G = nx.Graph()
    G.add_edge("H","C")
    G.add_edge("B","C")
    G.add_edge("A","B")
    G.add_edge("A","H")
    G.add_edge("B","D")
    G.add_edge("E","D")
    pos=nx.spring_layout(G,iterations=1000)
    nx.draw_networkx_nodes(G, pos, alpha=0.7)
    nx.draw_networkx_edges(G,pos,width=1)
    nx.draw_networkx_labels(G,pos,font_size=10,font_family='sans-serif')
    
    p.show()
    G.clear()
    