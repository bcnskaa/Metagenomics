# 2014 SKWoolf bcnskaa AT gmail DOT com

import os.path
import sys
import getopt
from collections import defaultdict


# Require biopython
# To run this script, path to biopython libraries has to be included in PYTHONPATH
from Bio import SeqIO
from Bio.Seq import Seq


# Path to circos
CIRCOS_HOME="~/tools/circos/"

"""
 Import all the contig
"""
def import_contig_info_from_maxbin(maxbin_dir, fna_fn_ext="fna"):
    print("Reading contig info from " + maxbin_dir + "...")
    
    fns = glob.glob(maxbin_path + "/*." + fn_ext)
    
    if len(fns) == 0:
        print("No *." + fna_fn_ext + " is found in " + maxbin_dir + ", abort now.")
        return None
    
    bin_groups = {}
    for j, fn in enumerate(fns):
        seq_idx = SeqIO.index(fn, "fasta")
        id = (fn[::-1].split("/",1)[0])[::-1]
        id = id.replace("."+fn_ext, "")
        
        selected_contigs = {k:len(seq_idx[k].seq) for k in seq_idx.keys()}
        
        bin_groups[id] = selected_contigs
        
    return bin_groups


"""

"""
def import_abundance_from_maxbin(maxbin_dir, summary_fn_ext="summary"):
    print("Reading abundances from " + maxbin_dir + "...")   
    




    

"""
 
"""
def others():
    
    # Generate karyotypes
    chroms = []
    #selected_contig_ids = {}
    selected_contig_ids = []
    for j, fn in enumerate(fns):
        seq_idx = SeqIO.index(fn, "fasta")
        id = (fn[::-1].split("/",1)[0])[::-1]
        id = id.replace("."+fn_ext, "")
        
        #color_code = colors[j]
        color_code = "black"
        selected_contigs = [k for k in seq_idx.keys() if len(seq_idx[k].seq) > chrom_len_threshold]
    
        selected_contig_ids.extend(selected_contigs)
        chroms.extend(["chr - " + k.replace("-", "_") + " " + id.replace("-", "_") + " 0 " + str(len(seq_idx[k])) + " " + color_code for i, k in enumerate(selected_contigs)])



     
    
"""
chrom_len_threshold = 55000
run_id = "SWH-xyl35_S3.002"
query_seq_len = 10000


# Prepare color codes from Circos's color palettes
color_fn = "/home/siukinng/tools/circos/etc/colors.conf"
lines = []


with open(color_fn) as IN:
    lines.extend(IN.read().splitlines())


lines = [l for l in lines if not l.startswith("<") and not l.startswith("#") and len(l) > 2]
colors = [(l.split("=",1)[0]).strip() for l in lines]
colors.remove("white")


from random import shuffle
shuffle(colors)


import glob
from Bio import SeqIO

# Maxbin path to subject sequences
maxbin_path = "."
maxbin_path = "/home/siukinng/samples/lab/P1/MaxBin"

fn_ext = "fna"
fns = glob.glob(maxbin_path + "/*." + fn_ext)

if len(fns) == 0:
    fn_ext = "fasta"
    fns = glob.glob(maxbin_path + "/*." + fn_ext)
    

# Generate karyotypes
chroms = []
#selected_contig_ids = {}
selected_contig_ids = []
for j, fn in enumerate(fns):
    seq_idx = SeqIO.index(fn, "fasta")
    id = (fn[::-1].split("/",1)[0])[::-1]
    id = id.replace("."+fn_ext, "")
    #id = id.replace(".", "_")
    color_code = colors[j]
    selected_contigs = [k for k in seq_idx.keys() if len(seq_idx[k].seq) > chrom_len_threshold]
    #selected_contig_ids[id] = selected_contigs
    selected_contig_ids.extend(selected_contigs)
    #chroms.extend(["chr - " + id + "." + str(i) + " " + str(i) + " 0 " + str(len(seq_idx[k])) + " " + color_code for i, k in enumerate(seq_idx) if len(seq_idx[k]) > chrom_len_threshold])
    #chroms.extend(["chr - " + id + "." + str(i) + " " + str(i) + " 0 " + str(len(seq_idx[k])) + " " + color_code for i, k in enumerate(selected_contigs)])
    #chroms.extend(["chr - " + k + " " + str(i) + " 0 " + str(len(seq_idx[k])) + " " + color_code for i, k in enumerate(selected_contigs)])
    chroms.extend(["chr - " + k.replace("-", "_") + " " + id.replace("-", "_") + " 0 " + str(len(seq_idx[k])) + " " + color_code for i, k in enumerate(selected_contigs)])

# Query chromosomome
chroms.extend(["chr - " + run_id.replace("-", "_") + " " + run_id.replace("-", "_") + " 0 " + str(query_seq_len) + " red"])

# Export karyotype to a file
with open(run_id + ".karyotype.txt", "w") as OUT:
    for chrom in chroms:
        print>>OUT, chrom



# Read the blast results, and construct circos links based on selected_contig ids

blast_fn = "/home/siukinng/samples/lab/GZ-cell_S1/900/MaxBin/combined-combined.bla.xml"
blast_fn = "/home/siukinng/samples/lab/SWH-xyl35_S3/Binning/Map_P1/SWH-xyl35_S3.001-P1_contig-MaxBin.All.xml"
blast_fn = "/home/siukinng/samples/lab/SWH-xyl35_S3/Binning/Map_P1/SWH-xyl35_S3.007-P1_contig-MaxBin.All.xml"

cutoff_range=1500


from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SearchIO

query_species = "SWH_xyl35-S3"

# Generate band
blast_res = SearchIO.index(blast_fn, "blast-xml")
links = []
for key in blast_res.keys():
    res = blast_res[key]
    if len(res.hits) > 0:
        for hit in res.hits:
            for hsp in hit.hsps:
                #if hsp.hit_id != hsp.query_id:
                range = hsp.query_range[1] - hsp.query_range[0]
                if range > cutoff_range and hsp.hit_id in selected_contig_ids:
                    print query_species + "\t1\t10\t" + hsp.hit_id + "\t" + str(hsp.hit_range[0]) + "\t" + str(hsp.hit_range[1])
                    #links.append(query_species + "\t1\t10\t" + hsp.hit_id + "\t" + str(hsp.hit_range[0]) + "\t" + str(hsp.hit_range[1]))
                    links.append(query_species.replace("-", "_") + " 1 " + str(query_seq_len) + " " + (hsp.hit_id).replace("-", "_") + " " + str(hsp.hit_range[0] + 1) + " " + str(hsp.hit_range[1]))


with open(run_id + ".links.txt", "w") as OUT:
    for link in links:
        print>>OUT, link



"""

def main(argv):
    print("")  



"""

"""
def generate_conf():
    print("Generating configuration file...")



"""

"""
def generate_karyotype():
    print("Generating configuration file...")



# Print the usage of this script 
def print_usage():
    print("A utility for generating Circos diagram.")
    print(" ")
    print("Usage:")
    print(" python " + sys.argv[0] + " READ_1.fq READ_2.fq [OUTPUT_PREFIX]")
    #print("  python -i BLAST-RESULT-INFILE -o FILTER-OUTFILE [-q] [-b BITSCORE-CUTOFF] [-l ALIGNMENT-LENGTH-CUTOFF] [-p PERCENTAGE-IDENTITY]")
    #print("      -i STRING  Input file generated by BLAST with -m6 option")
    print(" ")
    print(" ") 
    print("Ver 0.2b")
    


# Print status
def print_status(msg):
    caller_name = inspect.stack()[1][3]
    msg = "[ " + caller_name + " ] " + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "  " + msg
    print(msg)



def getOverlap(x, y):
    #print "x[1]=", x[1], ", y[1]=", y[1], ", x[0]=", x[0], ", y[0]=", y[0]
    return max(0, min(x[1], y[1]) - max(x[0], y[0]))



"""
 Extract 
     ~/tools/blast/bin/blastn -db combined.seq -query combined.seq -out combined-combined.bla.xml -outfmt 5 -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 3 

"""
def extract_links(blast_fn, cutoff_range=200, discard_self=None):
    blast_res = SearchIO.index(blast_fn, "blast-xml")
    
    query_species = "species1"
    hit_species = "species2"
    
    for key in blast_res.keys():
        res = blast_res[key]
        if len(res.hits) > 0:
            for hit in res.hits:
                print "hit:"
                for hsp in hit.hsps:
                    range = hsp.query_range[1] - hsp.query_range[0]
                    if range > cutoff_range:
                        #print hsp.query_id + "@" + str(hsp.query_range[0]) + "-" + str(hsp.query_range[1]) + "\t" + hsp.hit_id + "@" + str(hsp.hit_range[0]) + "-" + str(hsp.hit_range[1]) + "\t" + str(range)
                        #print query_species + "\t1\t10\t" + hit_species + "\t" + str(hsp.hit_range[0]) + "\t" + str(hsp.hit_range[1])
                        print query_species + "\t1\t10\t" + hsp.hit_id + "\t" + str(hsp.hit_range[0]) + "\t" + str(hsp.hit_range[1])



def generate_band(blast_fn, mapped_col="blue", unmapped_col="red"):
    with open(blast_fn) as IN:
        lines = IN.read().splitlines()
        
    bands = [[int(l.split()[8]), int(l.split()[9])] for l in lines]
    
    # Check reverse
    for i, b in enumerate(bands):
        if b[0] > b[1]:
            bands[i] = [b[1], b[0]]
            
    # Sort
    bands = sorted(bands, key=lambda v: v[0])
        
    # Combined overlapping regions
    combined_bands = []
    for b in bands:
        i = len(combined_bands)
        print(i)
        if i == 0:
            combined_bands.append([b[0], b[1]])
        else:
            print(str(i) + " "+ str((combined_bands[i - 1])[1]))
            if b[0] < (combined_bands[i - 1])[1]:
                print("overlapped")
                (combined_bands[i - 1])[1] = b[1]
            else:
                combined_bands.append([b[0], b[1]])
    
    karyotypes = []
    for b in combined_bands:
        band_id = "m" + str(len(karyotypes))
        karyotypes.append("band chr1 " + band_id + " " + band_id + " " + str(b[0]) + " " + str(b[1]) + " " + mapped_col)
    
    for i, b in enumerate(combined_bands):
        if i > 0:
            band_id = "u" + str(i-1)
            spos = (combined_bands[i - 1][1] + 1)
            epos = b[0]
            karyotypes.append("band chr1 " + band_id + " " + band_id + " " + str(spos) + " " + str(epos) + " " + unmapped_col)
     
    return karyotypes  



"""
from Bio import SeqUtils

seq = list(SeqIO.read("/disk/rdisk08/siukinng/db/BacteriaDB/all_fna/Propionibacterium_acnes_TypeIA2_P_acn33_uid80745/NC_016516.fna", "fasta"))

total_gc = SeqUtils.GC(seq)
step = 1000
gc_ctx = SeqUtils.GC_skew(seq, step)
output = ["chr1 " + str(i * step) + " " + str((i + 1) * step) + " " + str(gc + total_gc) for i, gc in enumerate(gc_ctx)]
with open("gc.txt", "w") as OUT:
    for gc in output:
        OUT.write(gc + "\n")


with open("NC_016516.ffn") as IN:
    lines = IN.read().splitlines()
    
genes = [((((l.split()[0]).split("|")[4]).replace(":","")).replace("c", "")).split("-") for l in lines if l.startswith(">")]
genes = [[int(g[0]), int(g[1])] for g in genes]

reverse_n = 0
for i, g in enumerate(genes):
    if g[1] < g[0]:
        genes[i] = [g[0], g[1]]
        reverse_n = reverse_n + 1

output = ["chr1 " + str(g[0]) + " " + str(g[1]) for g in genes]

with open("genes.txt", "w") as OUT:
    for g in output:
        OUT.write(g + "\n")

"""
# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    
    
    
"""

~/tools/blast/bin/blastn -query GZ-Cell_Y2_scaffold.fa -db GZ_scaffold.fa -outfmt 6 -out GZ-Cell_Y2-GZ.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2
~/tools/blast/bin/blastn -query GZ-Cell_Y2_scaffold.fa -db SWH_scaffold.fa -outfmt 6 -out GZ-Cell_Y2-SWH.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2
~/tools/blast/bin/blastn -query GZ-Xyl_Y2_scaffold.fa -db SWH_scaffold.fa -outfmt 6 -out GZ-Xyl_Y2-SWH.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2
~/tools/blast/bin/blastn -query GZ-Xyl_Y2_scaffold.fa -db GZ_scaffold.fa -outfmt 6 -out GZ-Xyl_Y2-GZ.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2


~/tools/blast/bin/blastn -query SWH-Cell_Y2_scaffold.fa -db GZ_scaffold.fa -outfmt 6 -out SWH-Cell_Y2-GZ.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2
~/tools/blast/bin/blastn -query SWH-Cell_Y2_scaffold.fa -db SWH_scaffold.fa -outfmt 6 -out SWH-Cell_Y2-SWH.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2
~/tools/blast/bin/blastn -query SWH-Xyl_Y2_scaffold.fa -db SWH_scaffold.fa -outfmt 6 -out SWH-Xyl_Y2-SWH.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2
~/tools/blast/bin/blastn -query SWH-Xyl_Y2_scaffold.fa -db GZ_scaffold.fa -outfmt 6 -out SWH-Xyl_Y2-GZ.bla -num_threads 16 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2


"""