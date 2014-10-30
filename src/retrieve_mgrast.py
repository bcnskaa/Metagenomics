from __future__ import division
import os
import sys
import ijson
import urllib2
import math


data_item_url = 'data.item.url'
data_item_file_name = 'data.item.file_name'


def main(arg):
    meta_ids = import_meta_ids(arg[0])
    if len(meta_ids) > 0:
        for meta_id in meta_ids:
            e = retrieve_from_mgrast(meta_id)
        
            if len(e) == 2:
                cmd = "wget " + e[data_item_url] + " -O " + e[data_item_file_name]
                print cmd
                os.system(cmd)
                
        

def import_meta_ids(fn):
    with open(fn) as IN:
        meta_ids = IN.read().splitlines()
    return meta_ids


def retrieve_from_mgrast(meta_id, stage_name="upload"):
    url = "http://api.metagenomics.anl.gov/1/download/mgm" + meta_id
    print("Retrieving " + meta_id + " from " + url)
    
    res = urllib2.urlopen(url)
    #page = res.read()
    
    #ip = ijson.parse(page)
    ip = ijson.parse(res)
    
    elements = {}
    state = 0
    for p, e, v in ip:
        #print str(p) + ":" + str(e)
        if state == 2:
            #print str(p) + ":" + str(e)
            if str(p) == data_item_file_name:
                elements[data_item_file_name] = str(v)
            if str(p) == data_item_url:
                elements[data_item_url] = str(v)

        if state == 1:
            if str(p) == 'data.item.stage_name' and str(v) == stage_name:
                state = 2
                
        if str(p) == 'data.item' and str(e) == 'start_map':
            #print("start_map")
            state = 1

        if str(p) == 'data.item' and str(e) == 'end_map':
            #print("end_map")
            state = 0
            
    return elements
        


"""

import fetch_mgrast
import glob

ids = fetch_mgrast.import_meta_ids("meta_ids.lst2")
fns = glob.glob("sequences/*.fna")
fns = [f.replace("sequences/mgm", "") for f in fns]
fns = [f.replace(".050.upload.fna", "") for f in fns]

ids = [id for id in ids if id not in fns]

fetch_mgrast.generate_run_cmd(ids)

"""
def generate_run_cmd(ids, outfn_prefix="run_cmd", chunk=10, file="050.1"):
    fns = [outfn_prefix + "." + str(i) + ".sh" for i in range(1, chunk)]
    
    if len(ids) == 0:
        print("No file is found.")
        return None
    
    size = int(math.ceil(len(ids) / chunk))
    steps = [i * size for i in range(0, chunk)]
    steps = steps + [len(ids)]
    
    count = 0
    i = 0
    for j in steps[1:len(steps)]:
        out_fn = outfn_prefix + "." + str(count) + ".sh" 
        print("Processing " + str(i) + ":" + str(j) + " to " + out_fn)
        with open(out_fn, "w") as OUT:
            for idx in range(i, j):
                OUT.write("echo \"Downloding " + ids[idx] + "\"\n")
                OUT.write("wget -q http://api.metagenomics.anl.gov/1/download/mgm" + ids[idx] + "?file=" + file + " -O mgm" + ids[idx] + ".upload.fna\n")
        i = j
        count = count + 1
        

    

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    
# 
# 
# """
#     Upload to MG-Rast
#     
#     webkey=EWbN2FHyJ5M3zrxgb2FDBNRSk
#     fn=P1_contigs.fa
#     curl -H "auth: $webkey"  -X POST -F "upload=@./$fn" "http://api.metagenomics.anl.gov/1/inbox/" > curl_output.txt
#     
# """
# 
# 
# def retrieve(argv): 
#     project_id = argv[0]
#     
#     print("Project ID=", project_id)
#     
#     cwd = os.getcwd()
#     os.mkdir(project_id)
#     os.chdir(project_id)
#     
#     # Import meta ids for the project
#     with open("../" + str(project_id) +".lst") as IN:
#         meta_ids = IN.read().splitlines()
#     meta_ids = [l for l in meta_ids if len(l) > 1]
#     
#     print("Meta_ids=" + str(len(meta_ids)))
#     
#     for meta_id in meta_ids:
#         print("Processing " + meta_id)
#         
#         cmd = "rm .listing"
#         os.system(cmd)
#         
#         cmd = "wget --no-remove-listing ftp://ftp.metagenomics.anl.gov/projects/" + str(project_id) + "/" + meta_id + "/raw/"
#         print(cmd)
#         os.system(cmd)
#        
#         with open(".listing") as IN:
#             lines = IN.read().splitlines()
#         
#         fn = [l.split()[8] for l in lines if ".fna.gz" in l][0]
#         
#         print("Retrieving " + fn)
#         
#         cmd = "wget ftp://ftp.metagenomics.anl.gov/projects/" + str(project_id) + "/" + meta_id + "/raw/" + fn
#         print(cmd)
#         os.system(cmd)
#         
#     os.chdir(cwd)
# 
# 
# """
#  Given an abundance summary produced by maxbin, abundance information
#  of sequence header of contig file (exported from idba_ud) can be included.
# """
# def prepare_seq_cov(maxbin_abund_fn, contig_fn, output_fn):
#     contig_fn = "../../idba_ud/contig.fa"
#     maxbin_abund_fn = "GZ_.abund"
#     
#     with open(maxbin_abund_fn) as IN:
#         lines = IN.read().splitlines()     
#     abunds = {v.split("\t")[0]:v.split("\t")[1] for v in lines}
#     
#     mod_seqs = []
#     for seq in SeqIO.parse(contig_fn, "fasta"):    
#         if seq.id in abund.keys():
#             a = str(int(round(float(abunds[seq.id]))))
#             seq.id = seq.id + "_[cov=" + a + "]"
#             seq.description = ""
#             mod_seqs.append(seq)
#     
#     OUT = open(output_fn, "w")
#     SeqIO.write(mod_seqs, OUT, "fasta")
#     OUT.close()
# 
# 
# 
# def retrieve_taxonomic_hit(mg_rast_id):
#     link = "http://metagenomics.anl.gov/metagenomics.cgi?page=MetagenomeOverview&metagenome=" + mg_rast_id + "&action=chart_export&name=organism_domain_hits&file=download." + mg_rast_id + ".organism_domain_hits"

    