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


-evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 90 -max_target_seqs 1 -num_threads 14


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
import process_EMIRGE
import glob
import os

#sample_ids = glob.glob("*")

#sample_ids = [f for f in sample_ids if os.path.isdir(f)]
sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
for sample_id in sample_ids:
    emirge_fa_fn = sample_id + "/" + sample_id + ".16S.renamed.fasta"
    if os.path.isfile(emirge_fa_fn):
        bins_dir = "/home/siukinng/MG/scaffolds_2000/" + sample_id + "_2000"
        process_EMIRGE.map_emirge_to_bins(bins_dir, emirge_fa_fn)



"""
def map_emirge_to_bins(bins_dir, emirge_fa_fn, separator="+"):
    import glob
    import os 

    cmd = "/home/siukinng/tools/blast/bin/makeblastdb -in " + emirge_fa_fn + " -dbtype nucl"
    os.system(cmd)
    
    outdir = os.path.dirname(emirge_fa_fn) + "/"
    
    bins_fns = glob.glob(bins_dir + "/*.fasta")
    for bins_fn in bins_fns:
        bla_fn = outdir + (os.path.basename(bins_fn)).replace(".fasta", "") + separator + (os.path.basename(emirge_fa_fn)).replace(".fasta", "") + ".bla"
        cmd = "/home/siukinng/tools/blast/bin/blastn -query " + bins_fn + " -db " + emirge_fa_fn + " -out " + bla_fn + " -outfmt 6 -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 89 -max_target_seqs 1 -num_threads 14"
        os.system(cmd)
        cmd = "python /home/siukinng/tools/scripts/filter_blast_res.py -q -i " + bla_fn
        os.system(cmd)
        


"""
import process_EMIRGE

process_EMIRGE.rename_16S()

sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]



"""
def rename_16S():
    import os
    
    out_dir = "."
    sample_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "GZ-Seed_Y0", "GZ-Cell_Y1", "GZ-Cell_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    for sample_id in sample_ids:
        print("Processing " + sample_id)
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
    fa_fn = sample_id + "/" + sample_id +".16S.fasta"
    process_EMIRGE.rename_16S_2(sample_ids_abbrev[sample_id], in_fn = fa_fn)
    
"""
def rename_16S_2(id_prefix, in_fn = "all_samples.16S.fasta", out_fn = None):
    import re
    from Bio import SeqIO
    sample_ids_abbrev = {"GZ-Xyl_Y2":"GX2", "GZ-Xyl_Y1":"GX1", "GZ-Seed_Y0":"GS0", "GZ-Cell_Y1":"GC1", "GZ-Cell_Y2":"GC2", "SWH-Xyl_Y2":"SX2", "SWH-Xyl_Y1":"SX1", "SWH-Seed_Y0":"SS0", "SWH-Cell_Y1":"SC1", "SWH-Cell_Y2":"SC2", "SWH-Cell55_Y2":"S52"}
    
    if out_fn is None:
        out_fn = in_fn.replace(".fasta",".renamed.fasta")
    
    seqs = SeqIO.index(in_fn, "fasta")
    #selected_ids = ["GZ-Xyl_Y2", "GZ-Xyl_Y1", "SWH-Xyl_Y2", "SWH-Xyl_Y1"]
    #selected_ids = ["GZ-Cell_Y2", "GZ-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell_Y1", "SWH-Cell55_Y2"]
    
    with open(out_fn, "w") as OUT:
        for seq_id in seqs.keys():
            #if seq_id.split("|")[0] not in selected_ids:
            #if seq_id.split("|")[0] not in sample_ids_abbrev.keys():
            #    continue
            #id = sample_ids_abbrev[seq_id.split("|")[0]] + "." + seq_id.split("|")[1]
            id = id_prefix + "." + seq_id.split("|")[0]
            print("exporting " + id)
            #seqs[seq_id].desc = id
            #seqs[seq_id].id = id
            OUT.write(">"+id +"\n" + str(seqs[seq_id].seq) + "\n")
            #SeqIO.write(seqs[seq_id], OUT, "fasta")




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
    
    
       
    