from __future__ import division

import glob
import os
import numpy
import sys
from Bio import SeqIO


sys.path.append(os.path.abspath("/home/siukinng/tools/scripts"))


import mg_pipeline
import pick_seq
import merge_fa


"""


import process_proteomics
process_proteomics.process_proteomics("identified_ids.lst", "../../db/all_Y1+2/SWH-Cell_Y1+Y2.faa")

process_proteomics.process_proteomics("identified_ids.OMSSA.lst", "../../db/all_Y1+2/SWH-Cell_Y1+Y2.faa")

process_proteomics.process_proteomics("identified_ids.COMET.lst", "../../db/all_Y1+2/SWH-Cell_Y1+Y2.faa")
process_proteomics.process_proteomics("identified_ids.COMET.lst", "../../db/all_Y1+2/GZ-Cell_Y1+Y2.faa")


"""
def process_proteomics(lst_fn, fa_fn, lst_ofn=None, map_to_nr=False, summary_follows_lst_format=False):
    [lst_ofn, nr_ids, export_id_n] = tidy_sids(lst_fn, lst_ofn)
    print("Extracting " + str(len(nr_ids)) + " sequences")
    
    # Extract the sequences specified in the list file (lst_fn) from the sequence db (fa_fn)
    [fa_ofn, export_fa_n] = extract_seq(nr_ids, fa_fn, lst_fn + ".fa")
    
    if map_to_nr:
        # Blast to the NCBI NR database 
        bla_ofn = blast_to_nr(fa_ofn)
    else:
        # Blast to the NCBI BioProject
        bla_ofn = blast_to_bacterial(fa_ofn)
        summary_fn = bla_ofn+".summary"
        map_dat = map_sid_to_gi(bla_ofn, summary_ofn=summary_fn)
        if summary_follows_lst_format:
            format_bla_to_lst_fn(map_dat, lst_fn, summary_ofn + ".reformated")
            
            
            
"""
 Provided lst_fn may not be well-formed, we need to reformat them (excluding blank line, split delimited lines).
"""
def tidy_sids(lst_fn, lst_ofn=None, selected_prefix=None):
    with open(lst_fn) as IN:
        lst = IN.read().splitlines()
    
    ids = []
    for l in lst:
        l = l.split(",")
        ids.extend(l)
    
    nr_ids = [l for l in list(set(ids)) if len(l) > 0]
    
    if selected_prefix is not None:
        nr_ids = [l for l in nr_ids if l.startswith(selected_prefix)]
    
    if lst_ofn is None:
        lst_ofn = lst_fn + ".nr"
    
    export_n = 0
    with open(lst_ofn, "w") as OUT:
        for nr_id in nr_ids:
            OUT.write(nr_id + "\n")
            export_n = export_n + 1
    print(str(export_n) + " row exported.")
    
    return [lst_ofn, nr_ids, export_n]



"""
 
"""
def extract_seq(id_list, fa_fn, fa_ofn, is_append=False):
    export_n = pick_seq.pick_seq(id_list, fa_fn, fa_ofn, is_append)
    print(str(export_n) + " exported.")
    return [fa_ofn, export_n]



"""
"""
def extract_seq_from_db_by_tax(species_lst, db_fa_fn, fa_ofn=None):
    seqs = SeqIO.index(db_fa_fn, "fasta")
    
    seq_ids = list(seqs.keys())
    
    if fa_ofn is None:
        fa_ofn = "selected_seq_ids.fa"   
    
    export_n = 0
    OUT = open(fa_ofn, "w")   
    selected_seq_ids = []
    for species in species_lst:
        for seq_id in seq_ids:
            desc = seqs[seq_id].description
            if species in desc:
                export_n = export_n + 1
                merge_fa.export_fasta(desc, str(seqs[seq_id].seq), OUT=OUT)
    OUT.close()
    
    print(str(export_n) + " sequences exported.")
        
    
    

"""
Do BLAST search to the NR db
"""
def blast_to_nr(fa_fn, out_fn=None, nr_db_fn="/home/siukinng/db/Markers/ncbi_nr/nr"):
    if out_fn is None:
        out_fn = fa_fn + ".bla"
    
    
    #mg_pipeline.blastp(fa_fn, nr_db_fn, outdir=".", outfn=out_fn, perc_identity=89, max_target_seqs=1, evalue=1e-30)
    mg_pipeline.blastp(fa_fn, nr_db_fn, outdir=".", outfn=out_fn, max_target_seqs=1, evalue=1e-30)
    
    return out_fn
    
    

"""
Do BLAST search to the NCBI BioProject db
"""
def blast_to_bacterial(fa_fn, out_fn=None, nr_db_fn="/home/siukinng/db/BacteriaDB/all_faa.gi.faa"):
    if out_fn is None:
        out_fn = fa_fn + ".bla"
        
    #mg_pipeline.blastp(fa_fn, nr_db_fn, outdir=".", outfn=out_fn, perc_identity=89, max_target_seqs=1, evalue=1e-30)
    mg_pipeline.blastp(fa_fn, nr_db_fn, outdir=".", outfn=out_fn, max_target_seqs=1, evalue=1e-30)
    
    return out_fn








"""
Import the blast file and extract sids that then map to the NR list file
"""
def map_sid_to_nr(bla_fn, out_fn=None, gi_fn="/home/siukinng/db/BacteriaDB/all_faa.gi.lst"):
    print("Not implemented")



"""
Reformat the map file to the list file
"""
def format_bla_to_lst_fn(map_dat, lst_fn, reformated_summary_ofn, delim=";", list_delim=","):
#     with open(summary_fn) as IN:
#         map_dat = IN.read().splitlines()
#     map_dat = {s.split("\t")[0]:s.split("\t")[1] for s in map_dat}
#     
    with open(lst_fn) as IN:
        lst_dat = IN.read().splitlines()
        
    OUT = open(reformated_summary_ofn, "w")
    for l in lst_dat:
        if len(l) == 0:
            OUT.write("\n")
            continue
        
        vals = l.split(list_delim)
        #print(l + "\t" + str(len(vals)))
        map_vals = ["Unknown_function" for v in vals]
        
        for i, val in enumerate(vals):
            if val in map_dat.keys():
                map_vals[i] = map_dat[val]
        OUT.write(l +"\t" + delim.join(map_vals) + "\n")
    
    OUT.close()

    return reformated_summary_ofn



"""
"""
def format_summary_fn_to_lst_fn(summary_fn, lst_fn, reformated_summary_ofn=None, delim=";", list_delim=","):
    print("Formating " + summary_fn + " to follow " + lst_fn)
    
    with open(summary_fn) as IN:
        map_dat = IN.read().splitlines()
    map_dat = {s.split("\t")[0]:s.split("\t")[1] for s in map_dat}
    
    #print(list(map_dat.keys())[0] + "=" + map_dat[list(map_dat.keys())[0]])
    
    if reformated_summary_ofn is None:
        reformated_summary_ofn = summary_fn + ".reformated"  
    
    return format_bla_to_lst_fn(map_dat, lst_fn, reformated_summary_ofn, delim, list_delim)
    
    
    

"""
Import the blast file and extract sids that then map to the gi list file
"""
def map_sid_to_gi(bla_fn, id_fn=None, summary_ofn=None, gi_fn="/home/siukinng/db/BacteriaDB/all_faa.gi.lst"): 
    print("Reading from " + gi_fn)
    with open(gi_fn) as IN:
        gi = IN.read().splitlines()
    gi_db = {g.split("\t")[0]: g.split("\t")[3] + ":" + g.split("\t")[2] for g in gi if len(g.split("\t")) == 4}
    
    
    print("Reading from " + bla_fn)
    map_lst = {}
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
        
    for b in bla:
        qid = b.split("\t")[0]
        sid = b.split("\t")[1]
        if qid not in bla_lst.keys():
            map_lst[qid] = []
            map_lst[qid].append(sid)
    
    if summary_ofn is None:
        summary_ofn = bla_fn + ".summary"
       
    print("Mapping...")
    OUT = open(summary_ofn,"w")
    
    for qid in map_lst.keys():
        print("Processing " + qid)
        gi_ctx = []
        for sid in map_lst[qid]:
            gi_ctx.append(gi_db[sid])
        OUT.write(qid + "\t" + "\t".join(gi_ctx) + "\n")
    
    OUT.close()
    
    return map_lst



"""

If the protein sequences were extracted from binned metagenomics data, this function
will map the protein sequences to its assigned bin group. If sequences do not belong to 
any bin group, "Not binned" will be marked.


# extract protein ids from binned faa, and generate the protein_id map file
rm protein_id2bin_id.map;for f in *.faa;do bin_id=${f};bin_id=${bin_id/.faa/};bin_id=${bin_id/GZ-Cell_Y2./GC2_};ids=$(cat $f | grep -e "^>" | sed "s/>//" |  cut -d" " -f1 | sed "s/scaffold_/GC2_/");for id in ${ids[*]};do echo -e "$id\t$bin_id" >> protein_id2bin_id.map;done;done

"""
def map_protein_ids2bin_id(summary_fn, map_fn="protein_id2bin_id.map", id_mapped_fn=None, delimiter="\t"):
    with open(map_fn) as IN:
        protein_id2bin_id_map = IN.read().splitlines()
    protein_id2bin_id_map = {p.split("\t")[0]:p.split("\t")[1] for p in protein_id2bin_id_map}
    
    with open(summary_fn) as IN:
        summary = IN.read().splitlines() 
    ids = {s.split(delimiter)[0]:s for s in summary}
        
    if id_mapped_fn is None:
        id_mapped_fn = summary_fn + ".mapped"
    
    OUT = open(id_mapped_fn, "w")
    for id in ids.keys():
        if id in protein_id2bin_id_map.keys():
            OUT.write(ids[id] + "\t" + protein_id2bin_id_map[id]+"\n")
        else:
            OUT.write(ids[id] + "\t" + "Not binned\n")        
    OUT.close()
        


"""

"""
def rename_seq_id_to_bin_id(bin_fa_fns, header_prefix="scaffold", sample_ids_abbrev={"GZ-Xyl_Y2":"GX2", "GZ-Xyl_Y1":"GX1", "GZ-Seed_Y0":"GS0", "GZ-Cell_Y1":"GC1", "GZ-Cell_Y2":"GC2", "SWH-Xyl_Y2":"SX2", "SWH-Xyl_Y1":"SX1", "SWH-Seed_Y0":"SS0", "SWH-Cell_Y1":"SC1", "SWH-Cell_Y2":"SC2", "SWH-Cell55_Y2":"S52"}):
    print("Combining bin files (" + str(len(bin_fa_fns)) + ")")
    
    bin_id_map = {}
    for bin_fa_fn in bin_fa_fns:
        print("Reading from " + bin_fa_fn)
        sample_id = os.path.basename(bin_fa_fn)
        bin_id = sample_id[::-1].split(".", 1)[1][::-1]
        sample_id = bin_id.split(".", 1)[0]
        new_sample_id = sample_ids_abbrev[sample_id]
        print(sample_id + " to " + new_sample_id)
        if sample_id in sample_ids_abbrev.keys():
            bin_id = bin_id.replace(sample_id, new_sample_id)
            #bin_id = bin_id.replace(".", "_")
            
            print("bin_id=" + bin_id)
            
            seqs = SeqIO.index(bin_fa_fn, "fasta")
            seq_ids = list(seqs.keys())
            
            m = {seq_id.replace(header_prefix, new_sample_id):bin_id for seq_id in seq_ids}
            bin_id_map.update(m)
        else:
            print(sample_id + " does not exist in the provided sample id list.")
    
    return bin_id_map
    


"""

If the fasta headers of the input fasta file (fa_fn) follow with the format used by NCBI, 
the species id can be retrieved.
 
"""
def extract_species_info_from_fasta(fa_fn, out_fn=None):
    import re
    from Bio import SeqIO
    
    print("Processing " + fa_fn)
    
    seqs = SeqIO.index(fa_fn, "fasta")
    seq_ids = list(seqs.keys())
    
    unknown_n = 0
    info_n = 0
    
    if out_fn is None:
        out_fn = fa_fn + ".gi"
    info_map = {}
    for seq_id in seq_ids:
        gi = seq_id.split("|")[1]
        desc = seqs[seq_id].description
        species = re.findall("\[(.+)\]", desc)
        if len(species) == 1:
            species = species[0]
            info_n = info_n + 1
        else:
            species = "Unknown"
            unknown_n = unknown_n + 1
        info_map[gi] = species

    sorted_map = sorted(info_map.items(), key=lambda x:x[1])
    OUT = open(out_fn, "w")
    for (gi, species) in sorted_map:   
        OUT.write(gi+"\t"+species+"\n")
    OUT.close()
    
    print("Export: " + str(info_n) + " (Unknown:" + str(unknown_n) + ")")
        
            
    
#     
# #bam_fns=`ls /scratch/bcnskaa/$TUMOR_TYPE/*.bam`
# bam_fns=`ls ./*.bam`
# 
# ref_fn="$HOME/share/db/UCSC/23_chrs.clean.fa"
# 
# sv_types=("INV")
# #sv_types=("DEL" "DUP" "INV")
# delly_path="$HOME/tools/CNV/delly"
# excluded_list_fn="$HOME/tools/CNV/delly-0.5.6/human.hg19.excl.tsv"
# #export OMP_NUM_THREADS=${#bam_fns[*]} 
# export OMP_NUM_THREADS=32
# 
# chr_ids=("10")
# 
# 
# for chr_id in ${chr_ids[*]};do
#     fn_list=""
#     for f in ${bam_fns[*]};do
#         out_fn=$(basename $f .bam)
#         out_fn="$out_fn.$chr_id.bam_tmp"
#         cmd="~/tools/samtools/samtools view $f $chr_id -b > $out_fn"
#         echo $cmd
#         eval $cmd
#         
#         cmd="~/tools/samtools/samtools index $out_fn"
#         echo $cmd
#         eval $cmd
#         fn_list=$fn_list" $out_fn"
#     done
# 
#     #Generate a command
#     for sv_type in ${sv_types[*]};do
#         out_fn="$TUMOR_TYPE"_$sv_type.$chr_id.vcf
#         touch $out_fn
#         cmd="$delly_path -t $sv_type -o $out_fn -x $excluded_list_fn -g $ref_fn $fn_list" 
#         #cmd="$delly_path -t $sv_type -o "$TUMOR_TYPE"_$sv_type.$chr_id.vcf -x $excluded_list_fn -g $ref_fn $fn_list" 
# 
#         echo $cmd
#         eval $cmd
#     done
# 
#     rm ./*.$chr_id.bam_tmp*
# done
#  
#      
#     