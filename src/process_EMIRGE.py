
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



# Generate a summary of EMIRGE results
for f in *;do id=${f};count=$(cat $f/$f.16S.fasta | grep -c ">");echo -e "$f\t$count";done



# Blast the EMIRGE results
for d in *;do ~/tools/blast/bin/blastn -query $d/$d.16S.fasta -db ~/db/Markers/GreenGene/gg_13_5.fasta -outfmt 6 -out $d/$d.16S-gg_13_5.bla -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 90 -max_target_seqs 2 -num_threads 14;done


# Link the fasta sequences
for f in ../../*Xyl_*;do id=${f##*/};ln -fs $f/$id.16S.fasta;done


# Make blast db
for f in *.fasta;do ~/tools/blast/bin/makeblastdb -in $f -dbtype nucl;done


sid="GZ"
qid="SWH"
~/tools/blast/bin/blastn -query $qid.16S.fasta -db $sid.16S.fasta -outfmt 6 -out $qid.16S-$sid.16S.bla -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 90 -max_target_seqs 2 -num_threads 8


-evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 90 -max_target_seqs 1 -num_threads 14


"""

def extract_tax_from_ncbi(acc_id, db="taxonomy"):
    from Bio import Entrez
    
    
    

def extract_acc_from_fasta_ids(fasta_fn):
    from Bio import SeqIO
    
    seqs = SeqIO.index(fasta_fn, "fasta")
    seq_ids = list(seqs.keys())
    
    acc_ids = [id.split("|")[1].split(" ")[0] for id in seq_ids]
    
    return acc_ids




def map_acc_2_taxids(acc_ids, tax_id_fn="tax_slv_ssu_nr_119.acc_taxid"):
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
    
  
    
def map_taxid_2_tax(acc_id_2_taxid_map, tax_tbl_fn="tax_ncbi_ssu_ref_119.txt"):
    tax_tbl = import_tax_tbl(tax_tbl_fn)
    
    acc_id_2_tax_map = {}
    for acc_id in acc_id_2_taxid_map.keys():
        tax_id = acc_id_2_taxid_map[acc_id]
        if acc_id in tax_tbl.keys():
            acc_id_2_tax_map[acc_id] = tax_tbl[tax_id]
    return acc_id_2_tax_map




"""
import test

bla_fn = "GZ-Cell_Y2.16S-gg_13_5.bla"
with open(bla_fn) as IN:
    bla = IN.read().splitlines()
bla2taxids = {b.split("\t")[0]:b.split("\t")[1] for b in bla}

tax_map = test.get_tax_from_greengene(bla2taxids)

"""
def get_tax_from_greengene(bla2taxids, greengene_tax_fn="/home/siukinng/db/Markers/GreenGene/gg_13_5_taxonomy.txt"):
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
import test
import glob


sample_ids = glob.glob(".")

sample_ids = [f for f in sample_ids if os.path.isdir(f)]
for sample_id in sample_ids:
    emirge_fa_fn = sample_id + "/" + sample_id + ".16S.fasta"
    if os.path.isfile(emirge_fa_fn):
        bins_dir = "/home/siukinng/MG/scaffolds_2000/" + sample_id + "_2000"
        test.map_emirge_to_bins(bins_dir, emirge_fa_fn)

"""
def map_emirge_to_bins(bins_dir, emirge_fa_fn):
    import glob
    import os 

    cmd = "/home/siukinng/tools/blast/bin/makeblastdb -in " + emirge_fa_fn + " -dbtype nucl"
    os.system(cmd)
    
    outdir = os.path.dirname(emirge_fa_fn) + "/"
    
    bins_fns = glob.glob(bins_dir + "/*.fasta")
    for bins_fn in bins_fns:
        bla_fn = outdir + (os.path.basename(bins_fn)).replace(".fasta", "") + "-"+ (os.path.basename(emirge_fa_fn)).replace(".fasta", "") + ".bla"
        cmd = "/home/siukinng/tools/blast/bin/blastn -query " + bins_fn + " -db " + emirge_fa_fn + " -out " + bla_fn + " -outfmt 6 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 90 -max_target_seqs 2 -num_threads 14"
        os.system(cmd)
        cmd = "python /home/siukinng/tools/scripts/filter_blast_res.py -q -i " + bla_fn
        os.system(cmd)
        
        
        
        

       
    