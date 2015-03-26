from __future__ import division
import os
import glob 
import numpy
import sys
from Bio import SeqIO



HOME = "/home/siukinng"
TOOLS_HOME = HOME + "/tools"
SCRIPTS_HOME = TOOLS_HOME + "/scripts"
#EMIRGE_HOME = TOOLS_HOME + "/EMIRGE"

sys.path.append(os.path.abspath(SCRIPTS_HOME))

import mg_pipeline




#def process_EMIRGE(sample_ids_abbrev={"GZ-Cell_Y2":"GC2"}, sample_dir="/home/siukinng/MG/scaffolds_2000/", EMIRGE_HOME=mg_pipeline.EMIRGE_HOME):
def process_EMIRGE(sample_ids_abbrev={"GZ-Xyl_Y2":"GX2", "GZ-Xyl_Y1":"GX1", "GZ-Seed_Y0":"GS0", "GZ-Cell_Y1":"GC1", "GZ-Cell_Y2":"GC2", "SWH-Xyl_Y2":"SX2", "SWH-Xyl_Y1":"SX1", "SWH-Seed_Y0":"SS0", "SWH-Cell_Y1":"SC1", "SWH-Cell_Y2":"SC2", "SWH-Cell55_Y2":"S52"}, sample_dir="/home/siukinng/MG/scaffolds_2000/", EMIRGE_HOME=mg_pipeline.EMIRGE_HOME):
    import glob
    import os

    #sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    #sample_ids_abbrev = {"GZ-Xyl_Y2":"GX2", "GZ-Xyl_Y1":"GX1", "GZ-Seed_Y0":"GS0", "GZ-Cell_Y1":"GC1", "GZ-Cell_Y2":"GC2", "SWH-Xyl_Y2":"SX2", "SWH-Xyl_Y1":"SX1", "SWH-Seed_Y0":"SS0", "SWH-Cell_Y1":"SC1", "SWH-Cell_Y2":"SC2", "SWH-Cell55_Y2":"S52"}    
    sample_ids = list(sample_ids_abbrev.keys())
    
    for sample_id in sample_ids:
        print("Processing " + sample_id)
        emirge_fa_fn = post_process_EMIRGE(dir=sample_id)
        
        if emirge_fa_fn is None:
            print(sample_id + " skipped")
            continue
        
        if not os.path.isfile(emirge_fa_fn):
            continue
        
        [renamed_emirge_fa_fn, id_map] = rename_16S_2(sample_ids_abbrev[sample_id], in_fn=emirge_fa_fn)    
        if not os.path.isfile(renamed_emirge_fa_fn):
            continue 
        
        print("Mapping " + renamed_emirge_fa_fn + " to bins")
        
        #emirge_fa_fn = sample_id + "/" + sample_id + ".16S.renamed.fasta"
        #emirge_fa_fn = sample_id + "/" + sample_id + ".16S.renamed.fasta"
        
        
        binning_dir = sample_dir
        if not binning_dir.endswith("/"): 
            binning_dir = binning_dir + "/"
        binning_dir = binning_dir + sample_id + "/"
        consolidated_bin_fa_fn = binning_dir + sample_id + ".fasta"
        bla_res = map_emirge_to_consolidated_bin(consolidated_bin_fa_fn, renamed_emirge_fa_fn)
        if bla_res is not None:
            print("Number of blast results: " + str(len(bla_res)))
            
        if consolidated_bin_fa_fn is None:
            print("Unable to process " + consolidated_bin_fa_fn)
            continue
        
        sids = [bla_res[k] for k in bla_res.keys()]
        
        
        print("Getting length of scaffolds...")
        scaffold_lens = mg_pipeline.get_seq_lens(consolidated_bin_fa_fn)
        #scaffold_ids = list(scaffold_lens.keys())
        
        #print("scaffold_ids['1']=" + scaffold_ids[1])
        #for sid in sids:
        #    print(sid)

        bin_map_res = mg_pipeline.map_scaffold_ids2_bin_id(sids, binning_dir)
        
        # Rename the bin_group into bin_group_abbrev
        bin_map_res = {k.replace(sample_id, sample_ids_abbrev[sample_id]):bin_map_res[k] for k in bin_map_res.keys()}
        
#         print("bla_res")
#         for k in bla_res.keys():
#             print(k + "=" + bla_res[k])
#         print("bin_map_res")
#         for k in bin_map_res.keys():
#             print(k + "=" + bin_map_res[k]) 
#         


        OUT = open("./" + sample_id + "/" + sample_id + ".16S.renamed.mapped2bin", "w")
        #map_emirge_to_bins(bins_dir, renamed_emirge_fa_fn, outdir="./"+sample_id)
        tax_map = map_emirge_to_greengene(renamed_emirge_fa_fn, outdir="./"+sample_id)
        for k in tax_map.keys():
            scaffold_id = ""
            scaffold_len = -1
            if k in bla_res.keys():
                scaffold_id = bla_res[k]
                if scaffold_id in scaffold_lens.keys():
                    scaffold_len = scaffold_lens[scaffold_id]
            bin_id = "Not_binned"
            if scaffold_id in bin_map_res.keys():
                bin_id = bin_map_res[scaffold_id]
            tax = tax_map[k]  
            print(k + "\t" + bin_id + "\t" + scaffold_id + "\t" + str(scaffold_len) + "\t" + tax)
            OUT.write(k + "\t" + bin_id + "\t" + scaffold_id + "\t" + str(scaffold_len) + "\t" + tax + "\n")
        OUT.close()
        
#         for k in bla_res.keys():
#             scaffold_id = bla_res[k]
#             bin_id = "Unknown"
#             if scaffold_id in bin_map_res.keys():
#                 bin_id = bin_map_res[scaffold_id]
#             tax = tax_map[k]
#             print(k + "\t" + bin_id + "\t" + tax)
                
                     
        
#         acc_ids = [id_map[id][1] for id in id_map.keys()]
#         acc_ids_2_taxid_map = map_acc_2_taxids(acc_ids)
#         if acc_ids_2_taxid_map is not None:
#             for k in acc_ids_2_taxid_map.keys():
#                 print(k + "=" + acc_ids_2_taxid_map[k])
#             acc_ids_2_tax_map = map_taxid_2_tax(acc_ids_2_taxid_map)
#             print("Tax Info of " + sample_id + ":")
#             for k in acc_ids_2_tax_map.keys():
#                 print(k + "=" + acc_ids_2_tax_map[k])
        
    


def run():
    dir="lab"
    insert_size=500
    insert_size_sd=50
    read_len=90  

    working_dir = "/home/siukinng/MG/EMIRGE"
    
    
    
    

"""


# Run EMIRGE
dir="env"
#dir="lab"
insert_size=500
insert_size_sd=50
read_len=90

cd /home/siukinng/MG/EMIRGE
for id in GZ;do
        #cmd="~/tools/EMIRGE/bin/emirge.py ./$id -1 ~/samples/"$dir"/$id/"$id"_1.trimmed.cleaned.fq -2 ~/samples/"$dir"/$id/"$id"_2.trimmed.cleaned.fq --phred33 -b ~/tools/EMIRGE/
db/SSURef_111_candidate_db.fasta -f ~/tools/EMIRGE/db/SSURef_111_candidate_db.fasta -l $read_len -s $insert_size_sd -i $insert_size -a 16"
        cmd="~/tools/EMIRGE/bin/emirge.py ./$id -1 ~/samples/"$dir"/$id/"$id"_1.fq -2 ~/samples/"$dir"/$id/"$id"_2.fq --phred33 -b ~/tools/EMIRGE/db/SSURef_111_candidate_db.fasta 
-f ~/tools/EMIRGE/db/SSURef_111_candidate_db.fasta -l $read_len -s $insert_size_sd -i $insert_size -a 16"
        eval $cmd
done


# Post-process EMIRGE results
for d in *;do ~/tools/EMIRGE/bin/emirge_rename_fasta.py $d/iter.40 > $d/$d.16S.fasta;done

# Merge all 16S sequences
cmd = 'echo -n "" > all_samples.16S.renamed.fasta'
os.system(cmd)
cmd= "for f in */*.16S.renamed.fasta;do cat $f >> all_samples.16S.renamed.fasta;done"
os.system(cmd)

# Generate a summary of EMIRGE results
cmd='for f in *;do id=${f};if [ -d $f ];then count=$(cat $f/$f.16S.fasta | grep -c ">");echo -e "$f\t$count";fi; done'
os.system(cmd)



# Blast the EMIRGE results
cmd='for d in *;do if [ -d $d ];then ~/tools/blast/bin/blastn -query $d/$d.16S.renamed.fasta -db ~/db/Markers/GreenGene/gg_13_5.fasta -outfmt 6 -out $d/$d.16S.renamed+gg_13_5.bla -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 89 -max_target_seqs 1 -num_threads 16;fi;done'
os.system(cmd)


# Link the fasta sequences
for f in ../../*Xyl_*;do id=${f##*/};ln -fs $f/$id.16S.fasta;done


# Make blast db
for f in *.fasta;do ~/tools/blast/bin/makeblastdb -in $f -dbtype nucl;done


sid="GZ"
qid="SWH"
~/tools/blast/bin/blastn -query $qid.16S.fasta -db $sid.16S.fasta -outfmt 6 -out $qid.16S-$sid.16S.bla -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 90 -max_target_seqs 2 -num_threads 8


#-evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 90 -max_target_seqs 1 -num_threads 14


"""

def extract_tax_from_ncbi(acc_id, db="taxonomy"):
    from Bio import Entrez
    
    
    
    
"""

"""
def extract_acc_from_fasta_ids(fasta_fn):
    from Bio import SeqIO
    
    seqs = SeqIO.index(fasta_fn, "fasta")
    seq_ids = list(seqs.keys())
    
    acc_ids = [id.split("|")[1].split(" ")[0] for id in seq_ids]
    
    return acc_ids



"""

"""
def map_acc_2_taxids(acc_ids, tax_id_fn=mg_pipeline.EMIRGE_HOME+"/db/tax_slv_ssu_nr_119.acc_taxid"):
    with open(tax_id_fn) as IN:
        taxids = IN.read().splitlines()
    taxids = {tid.split("\t")[0]: tid.split("\t")[1] for tid in taxids}
    
    acc_id_2_taxid_map = {i:"" for i in acc_ids}
    for i in acc_id_2_taxid_map.keys():
        if i in taxids.keys():
            acc_id_2_taxid_map[i] = taxids[i]
        else:
            print(i + " is not found.")

    return acc_id_2_taxid_map
    
    
    
def import_tax_tbl(tax_tbl_fn):
    with open(tax_tbl_fn) as IN:
        tax_tbl = IN.read().splitlines()
    tax_tbl = {l.split("\t")[1]:l for l in tax_tbl}
    return tax_tbl

  
   
def map_taxid_2_tax(acc_id_2_taxid_map, tax_tbl_fn=mg_pipeline.EMIRGE_HOME+"/db/tax_ncbi_ssu_ref_119.txt"):
    tax_tbl = import_tax_tbl(tax_tbl_fn)

    acc_id_2_tax_map = {}
    for acc_id in acc_id_2_taxid_map.keys():
        tax_id = acc_id_2_taxid_map[acc_id]
        if acc_id in tax_tbl.keys():
            acc_id_2_tax_map[acc_id] = tax_tbl[tax_id]
        #else:
        #    print("No tax info for " + acc_id + "(" + tax_id + ").")
    return acc_id_2_tax_map




def post_process_EMIRGE(dir=".", out_fn=None, EMIRGE_HOME=mg_pipeline.EMIRGE_HOME):
    import os
    import sys
    
    # Post-process EMIRGE results
    if not os.path.isdir(dir):
        return None
    if not os.path.isdir(dir+"/iter.40"):
        return None
    
    if out_fn is None:
        out_fn = dir + "/all.16S.fasta"
        
    cmd = EMIRGE_HOME + "/bin/emirge_rename_fasta.py " + dir + "/iter.40 > " + out_fn
    os.system(cmd)
    
    return out_fn
   
    
    
"""
import test

bla_fn = "GZ-Cell_Y2.16S-gg_13_5.bla"
with open(bla_fn) as IN:
    bla = IN.read().splitlines()
bla2taxids = {b.split("\t")[0]:b.split("\t")[1] for b in bla}

tax_map = test.get_tax_from_greengene(bla2taxids)

"""
def get_tax_from_greengene(bla2taxids, greengene_tax_fn=mg_pipeline.GREENGENE_DB_TAX):
    with open(greengene_tax_fn) as IN:
        greengene_db = IN.read().splitlines()
    greengene_db = {g.split("\t")[0]:g.split("\t")[1] for g in greengene_db}
    tax_map = {}
    for id in bla2taxids.keys():
        tax = "-"
        tax_id = bla2taxids[id]
        if tax_id in greengene_db.keys():
            tax = greengene_db[tax_id]
        tax_map[id] = tax
    return tax_map



"""
import process_EMIRGE
import glob
import os

#sample_ids = glob.glob("*")

#sample_ids = [f for f in sample_ids if os.path.isdir(f)]
sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
#sample_ids = ["GZ-Cell_Y2", "SWH-Cell_Y2", "GZ-Xyl_Y2"]
for sample_id in sample_ids:
    #emirge_fa_fn = sample_id + "/" + sample_id + ".16S.renamed.fasta"
    emirge_fa_fn = sample_id + "/" + sample_id + ".16S.renamed.fasta"
    
    if os.path.isfile(emirge_fa_fn):
        bins_dir = "/home/siukinng/MG/scaffolds_2000/" + sample_id + "_2000"
        process_EMIRGE.map_emirge_to_bins(bins_dir, emirge_fa_fn)



"""

def map_emirge_to_greengene(emirge_fa_fn, out_bla_fn=None, outdir=None, separator="+", green_gene_db=mg_pipeline.GREENGENE_DB):
    import glob
    import os 

    if not os.path.isfile(emirge_fa_fn):
        print("File not existed: " + emirge_fa_fn + ".")
        return None
    
    if outdir is None:
        outdir = os.path.dirname(emirge_fa_fn)
        if len(outdir) == 0:
            outdir = "."
        outdir = outdir + "/"
        
    if out_bla_fn is None:
        qid = os.path.basename(emirge_fa_fn)[::-1].split(".", 1)[1][::-1]
        sid = os.path.basename(green_gene_db)[::-1].split(".", 1)[1][::-1]
        out_bla_fn = qid + separator + sid + ".bla"
        
    bla_fn = mg_pipeline.blastn(query_fn=emirge_fa_fn, db_fn=green_gene_db, outdir=outdir, outfn=out_bla_fn, max_target_seqs=1, perc_identity=89, num_threads=14, evalue=1e-35)
    if bla_fn is None:
        print("Blasting " + emirge_fa_fn + " failed.")
        return None
        
    cmd = "python " + mg_pipeline.SCRIPTS_HOME + "/filter_blast_res.py -q -i " + bla_fn
    os.system(cmd)
    
    bla_fn = bla_fn + ".q.bla" 
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
    bla2taxids = {b.split("\t")[0]:b.split("\t")[1] for b in bla}
    
    tax_map = get_tax_from_greengene(bla2taxids)
    
    return tax_map


    
    
def map_emirge_to_consolidated_bin(consolidated_bin_fa_fn, emirge_fa_fn, bla_out_fn=None, outdir=None, separator="+"):
    import glob
    import os
    
    if not os.path.isfile(emirge_fa_fn):
        print("File not existed: " + emirge_fa_fn + ".")
        return None
    
    cmd = "/home/siukinng/tools/blast/bin/makeblastdb -in " + emirge_fa_fn + " -dbtype nucl"
    os.system(cmd)
    
    if outdir is None:
        outdir = os.path.dirname(emirge_fa_fn)
        if len(outdir) == 0:
            outdir = "."
        outdir = outdir + "/"
    
    if bla_out_fn is None:
        qid = os.path.basename(consolidated_bin_fa_fn)[::-1].split(".", 1)[1][::-1]
        sid = os.path.basename(emirge_fa_fn)[::-1].split(".", 1)[1][::-1]
        #bla_fn = (os.path.basename(consolidated_bin_fa_fn)).replace(bin_fa_fn_ext, "") + separator + (os.path.basename(emirge_fa_fn)).replace(bin_fa_fn_ext, "") + ".bla"
        bla_fn = qid + separator + sid + ".bla"
        
    bla_fn = mg_pipeline.blastn(query_fn=consolidated_bin_fa_fn, db_fn=emirge_fa_fn, outdir=outdir, outfn=bla_fn, max_target_seqs=1, perc_identity=89, num_threads=14, evalue=1e-35)
        
    if bla_fn is None:
        print("Blasting " + emirge_fa_fn + " against " + consolidated_bin_fa_fn + " failed.")
        return None
        
    cmd = "python " + mg_pipeline.SCRIPTS_HOME + "/filter_blast_res.py -q -i " + bla_fn
    os.system(cmd)
    

    #return bla_fn + ".q.bla"
    bla_fn = bla_fn + ".q.bla"

    
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
    bla = {b.split("\t")[1]:b.split("\t")[0] for b in bla}
    
    return bla
    


    

def map_emirge_to_bins(bins_dir, emirge_fa_fn, outdir=None, separator="+", bin_fa_fn_ext=".fasta"):
    import glob
    import os 

    if not os.path.isfile(emirge_fa_fn):
        print("File not existed: " + emirge_fa_fn + ".")
        return None

    cmd = "/home/siukinng/tools/blast/bin/makeblastdb -in " + emirge_fa_fn + " -dbtype nucl"
    os.system(cmd)
    
    if outdir is None:
        outdir = os.path.dirname(emirge_fa_fn)
        if len(outdir) == 0:
            outdir = "."
        outdir = outdir + "/"
    
    bins_fns = glob.glob(bins_dir + "/*" + bin_fa_fn_ext)
    for bins_fn in bins_fns:
        print("Mapping 16S sequences to " + bins_fn)
        #bla_fn = outdir + (os.path.basename(bins_fn)).replace(".fasta", "") + separator + (os.path.basename(emirge_fa_fn)).replace(".fasta", "") + ".bla"
        #cmd = "/home/siukinng/tools/blast/bin/blastn -query " + bins_fn + " -db " + emirge_fa_fn + " -out " + bla_fn + " -outfmt 6 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 89 -max_target_seqs 1 -num_threads 14"
        
        #cmd = BLAST_HOME + "/bin/blastn -query " + bins_fn + " -db " + emirge_fa_fn + " -out " + bla_fn + " -outfmt 6 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 89 -max_target_seqs 1 -num_threads 14"
        #os.system(cmd)
         
        #blastn(query_fn, db_fn, outdir=".",outfn=None, outfmt=6, num_threads=16, best_hit_score_edge=0.05, best_hit_overhang=0.25, perc_identity=80, max_target_seqs=2, evalue=1e-10):
        
        bla_fn = (os.path.basename(bins_fn)).replace(bin_fa_fn_ext, "") + separator + (os.path.basename(emirge_fa_fn)).replace(bin_fa_fn_ext, "") + ".bla"
        #bla_fn = mg_pipeline.blastn(query_fn=bins_fn, db_fn=emirge_fa_fn, outdir=outdir, outfn=bla_fn, max_target_seqs=1, perc_identity=89, num_threads=14)
        
        #if bla_fn is None:
        #    print("Blasting " + bins_fn + " failed.")
        #    continue
        
        #cmd = "python " + mg_pipeline.SCRIPTS_HOME + "/filter_blast_res.py -q -i " + bla_fn
        #os.system(cmd)
        
        map_emirge_to_consolidted_bin(bins_fn, emirge_fa_fn, bla_out_fn=bla_fn)
        
    return outdir
        


"""
import process_EMIRGE

process_EMIRGE.rename_16S()

sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]



"""
def rename_16S(out_dir = ".", sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]):
    import os
    
    out_dir = "."
    #sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    for sample_id in sample_ids:
        print("Processing " + sample_id)
        #fa_fn = out_dir + "/" + sample_id + "/" + sample_id + ".16S.fasta"
        fa_fn = out_dir + "/" + sample_id + "/" + sample_id + ".16S.fasta"
        fa_outfn = out_dir + "/" + sample_id + "/" + sample_id + ".16S.renamed.fasta"
        
        # cat SWH-Seed_Y0.16S.fasta | sed "s/>/>SWH-Seed_Y0|/" > SWH-Seed_Y0.16S.renamed.fasta
        cmd = "cat " + fa_fn + " | sed 's/>/>" + sample_id + "|/' > " + fa_outfn
        os.system(cmd)
        
        

"""
import process_EMIRGE
sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
sample_ids_abbrev = {"GZ-Xyl_Y2":"GX2", "GZ-Xyl_Y1":"GX1", "GZ-Seed_Y0":"GS0", "GZ-Cell_Y1":"GC1", "GZ-Cell_Y2":"GC2", "SWH-Xyl_Y2":"SX2", "SWH-Xyl_Y1":"SX1", "SWH-Seed_Y0":"SS0", "SWH-Cell_Y1":"SC1", "SWH-Cell_Y2":"SC2", "SWH-Cell55_Y2":"S52"}    
for sample_id in sample_ids:
    print("Processing " + sample_id)
    #fa_fn = sample_id + "/" + sample_id +".16S.fasta"
    fa_fn = sample_id + "/all.16S.fa"
    
    process_EMIRGE.rename_16S_2(sample_ids_abbrev[sample_id], in_fn = fa_fn)
    
"""
def rename_16S_2(id_prefix, in_fn = "all_samples.16S.fasta", out_fn = None):
    import re
    from Bio import SeqIO
    #sample_ids_abbrev = {"GZ-Xyl_Y2":"GX2", "GZ-Xyl_Y1":"GX1", "GZ-Seed_Y0":"GS0", "GZ-Cell_Y1":"GC1", "GZ-Cell_Y2":"GC2", "SWH-Xyl_Y2":"SX2", "SWH-Xyl_Y1":"SX1", "SWH-Seed_Y0":"SS0", "SWH-Cell_Y1":"SC1", "SWH-Cell_Y2":"SC2", "SWH-Cell55_Y2":"S52"}
    
    if not os.path.isfile(in_fn):
        print("File not existed: " + in_fn + ".")
        return None
    
    if out_fn is None:
        if in_fn.endswith(".fasta"):
            out_fn = in_fn.replace(".fasta",".renamed.fasta")
        elif in_fn.endswith(".fa"):    
            out_fn = in_fn.replace(".fa",".renamed.fa")
        else:
            out_fn = in_fn + ".renamed.fasta"
        
    seqs = SeqIO.index(in_fn, "fasta")
    #selected_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "SWH-Xyl_Y2", "SWH-Xyl_Y1"]
    #selected_ids = ["GZ-Cell_Y2", "GZ-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell_Y1", "SWH-Cell55_Y2"]
    
    id_map = {}
    with open(out_fn, "w") as OUT:
        for seq_id in seqs.keys():
            #if seq_id.split("|")[0] not in selected_ids:
            #if seq_id.split("|")[0] not in sample_ids_abbrev.keys():
            #    continue
            #id = sample_ids_abbrev[seq_id.split("|")[0]] + "." + seq_id.split("|")[1]
            id = id_prefix + "." + seq_id.split("|")[0]
            
            rid = seq_id.split("|")[0]
            acc = seq_id.split("|")[1].split(" ")[0]
            
            id_map[id] = [rid, acc, seq_id]
             
            #print("exporting " + id)
            #seqs[seq_id].desc = id
            #seqs[seq_id].id = id
            OUT.write(">"+id +"\n" + str(seqs[seq_id].seq) + "\n")
            #SeqIO.write(seqs[seq_id], OUT, "fasta")
            
    return [out_fn, id_map]


def rename_plotfile_name_to_greengene(bla_infn, infn="plotfile", gg_tax_fn="/home/siukinng/db/Markers/GreenGene/gg_13_5_taxonomy.txt"):
    import re
    with open(infn) as IN:
        intree = IN.read().splitlines()
        
    print("Reading from " + gg_tax_fn)
    with open(gg_tax_fn) as IN:
        gg = IN.read().splitlines()
    gg = {g.split("\t")[0]:g.split("\t")[1].replace(" ", "") for g in gg}
    
    with open(bla_infn) as IN:
        gg2 = IN.read().splitlines()
    gg2 = {g.split("\t")[0]:g.split("\t")[0]+" ("+g.split("\t")[1]+")" for g in gg2}
    
    gg.update(gg2)
    
    with open(infn + ".renamed", "w") as OUT:
        for l in intree:
            m = re.findall("\((.+?)\) show", l)
            if len(m) > 0:
                m = m[0]
                if m in gg.keys():
                    OUT.write("(" + m + "=" + gg[m] + ") show\n")
                else:
                    OUT.write(l + "\n")
            else:
                OUT.write(l + "\n")

        

def rename_intree_name_to_greengene(infn="intree", gg_tax_fn="/home/siukinng/db/Markers/GreenGene/gg_13_5_taxonomy.txt"):
    import re
    with open(infn) as IN:
        intree = IN.read()
        
    print("Reading from " + gg_tax_fn)
    with open(gg_tax_fn) as IN:
        gg = IN.read().splitlines()
    gg = {g.split("\t")[0]:g.split("\t")[1].replace(" ", "") for g in gg}
    
    matches = re.findall("\(+(.+?)\:", intree)
    if matches is not None:
        for m in matches:
            if m in gg.keys():
                print("Replacing " + m + " to " + gg[m])
                intree = intree.replace("("+m+":", "("+gg[m]+":")
    
    with open(infn + ".renamed", "w") as OUT:
        OUT.write(intree)



"""

"""    
def extract_16S_from_greengene(gg_seq_fn="/home/siukinng/db/Markers/GreenGene/gg_13_5.fasta"):
    from Bio import SeqIO
    
    sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    
    for sample_id in sample_ids:
        print("Processing " + sample_id)
        bla_fn = sample_id + "/" + sample_id + ".16S.renamed+gg_13_5.bla"
        with open(bla_fn) as IN:
            bla = IN.read().splitlines()
        sids = [b.split("\t")[1] for b in bla]
        seqs = SeqIO.index(gg_seq_fn, "fasta")
        out_fn = bla_fn + ".sid.fasta"
        with open(out_fn, "w") as OUT:
            for sid in sids:
                if sid in seqs.keys():
                    seq = str(seqs[sid].seq)
                    OUT.write(">" + sid + "\n" + seq + "\n")

    
    #cmd = "cat GZ-Cell_Y1.16S.renamed+gg_13_5.bla  | sed 's/\t/ /g' | cut -d' ' -f2 > GZ-Cell_Y1.16S.renamed+gg_13_5.bla.sid"     
    
    
       
    