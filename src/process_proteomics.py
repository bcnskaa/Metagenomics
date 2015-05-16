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

# Installation of TPP with Python

# Download the most updated version of Python from
# https://mail.python.org/pipermail/tutor/2002-March/012903.html
./configure --enable-shared
make -j16 altinstall prefix=/home/siukinng/tools/Python-2.7.9 exec-prefix=/home/siukinng/tools/Python-2.7.9

# If there is another version of Python installed in system, we can make an alias 
ln -s /home/siukinng/tools/Python-2.7.9 /home/siukinng/tools/Python
ln -s /home/siukinng/tools/Python/bin/python2.7 /home/siukinng/tools/Python/bin/python

# Add the following lines to .bash_profile 
# alias python=$HOME/tools/Python/bin/python2.7
# export PYTHONPATH=$HOME/tools/Python/lib:$HOME/Python/lib/python2.7/site-packages
# export LD_LIBRARY_PATH=$HOME/tools/Python/lib

# Download the libxml2 library
CFLAGS=-L/home/siukinng/tools/Python/lib ./configure --with-python=/home/siukinng/tools/Python/bin/python --prefix=/home/siukinng/tools/lib
make -j16

# export LD_LIBRARY_PATH=$HOME/tools/Python/lib:$HOME/tools/lib/lib

# 



# Download the precompiled binnary at https://groups.google.com/forum/#!topic/spctools-discuss/qOR0pabxsWo
# https://depot.galaxyproject.org/package/linux/x86_64/tpp/tpp-4.8.0-Linux-x86_64.tar.gz

"""



  
"""


import process_proteomics
process_proteomics.process_proteomics("identified_ids.lst", "../../db/all_Y1+2/SWH-Cell_Y1+Y2.faa")

process_proteomics.process_proteomics("identified_ids.OMSSA.lst", "../../db/all_Y1+2/SWH-Cell_Y1+Y2.faa")

process_proteomics.process_proteomics("identified_ids.COMET.lst", "../../db/all_Y1+2/SWH-Cell_Y1+Y2.faa")
process_proteomics.process_proteomics("identified_ids.COMET.lst", "../../db/all_Y1+2/GZ-Cell_Y1+Y2.faa")


# with map file available
process_proteomics.process_proteomics("interact-GZ-Xyl.ipro.prot.xls", "GZ-Xyl+DECOY.faa", "/home/siukinng/MG/m8/proteins/renamed/all_samples+nr.renamed.m8")

"""
def process_proteomics(lst_fn, fa_fn=None, map_fn=None, lst_ofn=None, map_to_nr=False, summary_follows_lst_format=False):

    if map_fn is not None:
        print("Reading the ID to GI map file from " + map_fn)
        with open(map_fn) as IN:
            qid2gi_map = IN.read().splitlines()
        
        # m8 file containing PIDs to GIs
        qid2gi_map = [m.split("\t") for m in qid2gi_map]
            
#         # Check if the sid in gi|...|ref|... format
#         chk_l = qid2gi_map[0]
#         if chk_l[1].startswith("gi|"):
#             print("Extracting GIs...")
#             for s in qid2gi_map:
#                 s[1] = s[1].split("|")[1]
#             print("Number of GI cleaned: " + str(len(qid2gi_map)))
#         
#         # Extract ID to GI pairs
#         print("Extracting PID-GI pairs...")
#         pid_n = 0
#         pid2gi_pair_n = 0
#         
#         # Protein IDs to GIs 
#         pid2gi_map = {}
#         for m in qid2gi_map:
#             pid = m[0]
#             gi = m[1]
#             try:
#                 pid2gi_map[pid].append(gi)
#             except:
#                 pid2gi_map[pid] = []
#                 pid_n += 1
#                 pid2gi_map[pid].append(gi)
#             
#             pid2gi_pair_n += 1
#             
#         print("Number of PIDs = " + str(pid_n) + ", number of pid-gi pairs = " + str(pid2gi_pair_n))
#         
        
        #print("Importing the GI to Description map")
        #gi2ctx_map = import_gi2ctx_map()

        # def map_sid_to_gi_ctx(sid2gi_map, gi2ctx_map=None):
        #print("Mapping sids to GI description")
        #map_dat = map_id_to_gi_ctx(pid2gi_map, qid2gi_map)
        export_summary(qid2gi_map, lst_fn)
        

    else:    
        [lst_ofn, nr_ids, export_id_n] = tidy_sids(lst_fn, lst_ofn)
        print("Extracting " + str(len(nr_ids)) + " sequences")
        # Extract the sequences specified in the list file (lst_fn) from the sequence db (fa_fn)
        [fa_ofn, export_fa_n] = extract_seq(nr_ids, fa_fn, lst_fn + ".fa")
    
        if map_to_nr:
            # Blast to the NCBI NR database 
            bla_ofn = blast_to_nr(fa_ofn)
        else:
            # Blast to the NCBI BioProject
            print("Blasting...")
            bla_ofn = blast_to_bacterial(fa_ofn)

        
            summary_fn = bla_ofn+".summary"
            # Mapping sid to gi
            print("Mapping sid to gi")
            map_dat = map_sid_to_gi(bla_ofn, summary_ofn=summary_fn)
    
        if summary_follows_lst_format:
            format_bla_to_lst_fn(map_dat, lst_fn, lst_fn + ".summary.reformated")
            
            
            
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
    #export_n = 17063
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
    
    mg_pipeline.blastp(fa_fn, nr_db_fn, outdir=".", outfn=out_fn, max_target_seqs=1, evalue=1e-30, perc_identity=49)
    
    return out_fn
    
    

"""
Do BLAST search to the NCBI BioProject db
"""
def blast_to_bacterial(fa_fn, out_fn=None, nr_db_fn="/home/siukinng/db/BacteriaDB/all_faa.gi.faa"):
    if out_fn is None:
        out_fn = fa_fn + ".bla"
        
    #mg_pipeline.blastp(fa_fn, nr_db_fn, outdir=".", outfn=out_fn, perc_identity=49, max_target_seqs=1, evalue=1e-30)
    
    ####
    #mg_pipeline.blastp(fa_fn, nr_db_fn, outdir=".", outfn=out_fn, max_target_seqs=1, evalue=1e-30)
    
    return out_fn





"""
Import the blast file and extract sids that then map to the NR list file
"""
def map_sid_to_nr(bla_fn, out_fn=None, gi_fn="/home/siukinng/db/BacteriaDB/all_faa.gi.lst"):
    print("Not implemented")



"""
Reformat the list file
"""
def format_bla_to_lst_fn(map_dat, lst_fn, reformated_summary_ofn, delim=";", protein_id_delim=",", list_delim="\t"):
#     with open(summary_fn) as IN:
#         map_dat = IN.read().splitlines()
#     map_dat = {s.split("\t")[0]:s.split("\t")[1] for s in map_dat}
#     

    # Read in the list file
    with open(lst_fn) as IN:
        lst_dat = IN.read().splitlines()
        
    # Prepare the outfile
    OUT = open(reformated_summary_ofn, "w")
  
    processed_n = 0
    skipped_n = 0  
    # Go through every line of the list file 
    for i, l in enumerate(lst_dat):
        # Skip the first line
        if i == 0:
            OUT.write(l + "\n")
            continue
        
        # Skip the line if it contains nothing
        if len(l) == 0:
            OUT.write("\n")
            continue
        
        # delimite the line
        lst_vals = l.split(list_delim)
        
        protein_ids = lst_vals[2].split(protein_id_delim)
        #print(l + ":=" + str(len(lst_vals)))
        map_annotations = ["Unknown_function" for pid in protein_ids]
        
        
        
        # Map the annotation to pid
        # pid == qid
        for i, pid in enumerate(protein_ids):
            processed_n += 1
            try:
                map_annotations[i] = map_dat[pid][0]
            #if pid in map_dat.keys():
            #    # Only select the first description
            #    map_annotations[i] = map_dat[pid][0]
            #else:
            except:
                skipped_n += 1
                
        # Replace the descriptions (column 9) in list file with updated GI description
        lst_vals[9] = ",".join(map_annotations)
        
        OUT.write(list_delim.join(lst_vals) + "\n")
        
        #OUT.write(l +"\t" + delim.join(map_vals) + "\n")
    
    print("Number of processed: " + str(processed_n) + ", number of skipped: " + str(skipped_n))
    
    OUT.close()

    return reformated_summary_ofn



"""
"""
def export_summary(map_dat, lst_fn, summary_ofn=None, reformated_summary_ofn=None, delim=";", list_delim=","):
    if summary_ofn is None:
        summary_ofn = lst_fn + ".summary"
    if reformated_summary_ofn is None:
        reformated_summary_ofn = summary_ofn + ".reformated"  
    
    print("Formating " + summary_ofn + " to follow " + lst_fn)

    map_dat_lst = {}
    for m in map_dat:
        try:
            map_dat_lst[m[0]].append(m[1])
        except:
            map_dat_lst[m[0]] = [m[1]]
  
  
    print("Exporting results to " + summary_ofn)
    skipped_n = 0
    processed_n = 0
    p_value_threshold = get_p_value(lst_fn.replace(".prot.xls", ".pep.summary.txt"))
    ids = pick_protein(lst_fn, p_val_threshold=p_value_threshold)

    with open(summary_ofn, "w") as OUT:
        OUT.write("## Summary generated from " + lst_fn + "; P-value threshold=" + str(p_value_threshold) + "\n")
        for id in ids:
            try:
                OUT.write(id + "\t" + map_dat_lst[id][0] + "\n")   
                processed_n += 1    
            except:
                skipped_n += 1
    print("Number of items exported to " + summary_ofn + ": " + str(processed_n - skipped_n) + " (" + str(skipped_n) + " skipped)")
    
    
#     with open(summary_fn) as IN:
#         map_dat = IN.read().splitlines()
#    map_dat = {s.split("\t")[0]:s.split("\t")[1] for s in map_dat}

            
 
    print("Number of map data: " + str(len(map_dat)))
#    print(list(map_dat.keys())[0] + "=" + map_dat[list(map_dat.keys())[0]])
    
    
    return format_bla_to_lst_fn(map_dat_lst, lst_fn, reformated_summary_ofn, delim, list_delim)
    


def import_gi2ctx_map(gi_fn="/home/siukinng/db/Markers/ncbi_nr/gi2ctx.map"):
#def import_gi2ctx_map(gi_fn="/home/siukinng/db/BacteriaDB/all_faa.gi.lst"):
    with open(gi_fn) as IN:
        gi = IN.read().splitlines()
    #gi2ctx_map = {int(g.split("\t")[0]): g.split("\t")[3] + ":" + g.split("\t")[2] for g in gi if len(g.split("\t")) == 4}
    gi2ctx_map = [g.split("\t") for g in gi]
    #gi2ctx_map = {int(g[0]): g[1] for g in gi2ctx_map}
    
    print("Number of GI: " + str(len(gi2ctx_map)))
    
    return gi2ctx_map



"""
sid2gi_map is a dictionary with sid as key and a list of associated as value
"""
def map_id_to_gi_ctx(id2gi_map, gi2ctx_map=None):
    print("Processing " + str(len(id2gi_map)) + " sids")
    
    if gi2ctx_map is None:
        print("Importing the GI to Description map")
        gi2ctx_map = import_gi2ctx_map()
        
    # Sid mapped to GI descriptions
    #sid_to_gi_ctx = {sid for sid in sid2gi_map.keys()}
    id_to_gi_ctx = {id:[int(gi) for gi in id2gi_map[id]] for id in id2gi_map.keys()}
    
    processed_n = 0
    skipped_n = 0
    print("Mapping ID to GI description...")
    for id, gi_ctx in id_to_gi_ctx.items():  
        #gi_ctx = ["Unknown_function" for gi in sid2gi_map[sid]]

        # Iterate all mapped gi
        for i, gi in enumerate(gi_ctx):
            gi_ctx[i] = "Unknown_function"
            try:
                gi_ctx[i] = gi2ctx_map[gi]
            except:
                print("GI not found: " + str(gi))
                skipped_n += 1
            processed_n += 1  
            
    print("Number of records processed: " + str(processed_n) + " (" + str(skipped_n) + " skipped)")    
    
    return id_to_gi_ctx



"""
Import the blast file and extract sids that are then mapped to the gi list file
"""
def map_sid_to_gi(bla_fn, id_fn=None, summary_ofn=None, gi_map=None): 
    print("Reading GIs of protein sequences from " + gi_fn)
    if gi_map is None:
        gi_map = import_gi_map()
        
    
    print("Reading BLAST results from " + bla_fn)
    map_lst = {}
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
        
    # Blast results
    bla_lst = {}
    
    # Sid mapped to GI descriptions
    qid_to_gi_ctx = {}
    
    for b in bla:
        qid = b.split("\t")[0]
        sid = b.split("\t")[1]
        # Maybe more than more sids hit by qid
        if qid not in bla_lst.keys():
            map_lst[qid] = []
            map_lst[qid].append(sid)
            qid_to_gi_ctx[qid] = []
    
    if summary_ofn is None:
        summary_ofn = bla_fn + ".summary"
       
    print("Mapping...")
    OUT = open(summary_ofn,"w")
    
    for qid in map_lst.keys():
        print("Processing " + qid)
        gi_ctx = []
        # Iterate all mapped sids
        for sid in map_lst[qid]:
            gi_ctx.append(gi_map[sid])
            qid_to_gi_ctx[qid].append(gi_map[sid])
            
        OUT.write(qid + "\t" + "\t".join(gi_ctx) + "\n")
    
    OUT.close()
    
    return qid_to_gi_ctx
    #return map_lst



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
    


def import_result_map(map_fn, discard_not_evidence_row = True, min_percent_coverage=29.0):
    with open(map_fn) as IN:
        map_lst = IN.read().splitlines()
    del map_lst[0]
    map_lst = [m.split("\t") for m in map_lst if len(m) > 15]
    
    if discard_not_evidence_row:
        map_lst = [m for m in map_lst if m[13] == 'Y']
        
    map_lst = [m for m in map_lst if m[13] == 'Y']
    
    map_lst = [m for m in map_lst if float(m[4]) >= min_percent_coverage]
    
    if len(map_lst) == 0:
        return {}
    
    nr_map_lst = {}
    entry_lis = []
    for m in map_lst:
        entry_id = m[0]
        if entry_id not in entry_list:
            entry_list.append(entry_id)
            pep_ids = m[2].split(",")
            
            nr_map_lst.update({pid:m for pid in pep_ids})

    return nr_map_lst




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


def get_p_value(iprot_summary_fn):
    p_val_threshold = 0.0
    with open(iprot_summary_fn) as IN:
        summary = IN.read().splitlines()
        p_val_thread = [s for s in summary if s.startswith("P threshold for protein group FDR ")]
    
    if len(p_val_thread) == 0:
        print("Unable to extract recommended P-value threshold from " + iprot_summary_fn)
        print("default threshold=0.95 will be used.")
        p_val_threshold = 0.95
    else:
        p_val_threshold = float(p_val_thread[0].split("=")[1])
        print("Recommended P-value threshold=" + str(p_val_threshold))  
    return p_val_threshold



"""
Usage: 

import process_proteomics

process_proteomics.pick_protein()

"""
def pick_protein(lst_fn, iprot_summary_fn=None, p_val_threshold=0.95):
    if iprot_summary_fn is not None:
        p_val_threshold=get_p_value(iprot_summary_fn)
    
    with open(lst_fn) as IN:
        lst = IN.read().splitlines()
    lst = [l.split("\t") for l in lst if len(l) > 1 and not l.startswith("entry no.")]
    
    lst = [l for l in lst if float(l[4]) >= p_val_threshold]
    
    # Nr_IDs
    ids = []
    for l in lst:
        ids.extend(l[2].split(","))
    
    nr_ids = list(set(ids))
    
    return nr_ids



def select_m8(ids, selected_m8_ofn, m8_fn="/home/siukinng/MG/m8/proteins/renamed/all_samples+nr.renamed.m8"):
    with open(m8_fn) as IN:
        m8 = IN.read().splitlines()
    m8 = {m.split("\t")[0]: m for m in m8}
    
    skipped_n = 0
    selected_m8 = []
    for id in ids:
        try:
            selected_m8.append(m8[id])
        except:
            skipped_n += 1
    print("Number of selected m8: "+ str(len(selected_m8)) + " (" + str(skipped_n) + " skipped)")
    
    with open(selected_m8_ofn, "w") as OUT:
        for m in selected_m8:
            OUT.write(m + "\n")
            

"""

"""
def extract_faa(m8_fn, seq_map, out_fn=None):
    print("Reading from " + m8_fn)
    with open(m8_fn) as IN:
        m8 = IN.read().splitlines()
    m8 = [m.split("\t") for m in m8]
    print("Number of records read: " + str(len(m8)))
    
    m8_sids = [m[1] for m in m8]
    
    chk_l = m8_sids[0]
    
    if chk_l[1].startswith("gi|"):
        print("Extracting GIs...")
        for s in m8_sids:
            s[1] = s[1].split("|")[1]
        print("Number of GI cleaned: " + str(len(m8_sids)))

    m8_sids = [int(gi) for gi in m8_sids]
    
    
    nr_sids = list(set(m8_sids))
    print("Number of non-redundant ids: " + str(len(nr_sids)))
    
    print("Searching sequences...")
    skipped_n = 0
    processed_n = 0
    selected_seqs = {}
    for sid in nr_sids:
        try:
            selected_seqs[sid] = seq_map[sid]
        except:
            skipped_n += 1
    print("Number of processed: " + str(processed_n) + " (" + str(skipped_n) + " skipped)")

    if out_fn is None:
        out_fn = m8_fn + ".faa"
    
    exported_n = 0
    with open(out_fn, "w") as OUT:
        for sid in selected_seqs.keys():
            OUT.write(">" + sid + "\n" + selected_seqs[sid] + "\n")
            exported_n += 1
    
    print(str(exported_n) + " exported")
    
    return out_fn
            


# Main 
def main(argv):
    process_proteomics(argv[0], argv[1], summary_follows_lst_format=True)



# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])


def bench_marking():
    import time
    from random import randint

    len = 40000000
    test_n = 2000000
    a = {"a" + str(i):"a"+str(i) for i in xrange(len)}
    b = {i:"a"+str(i) for i in xrange(len)}

    t = [randint(0,len) for i in xrange(test_n)]

    start_time = time.time()
    for i in t:
        b[i]
    print("--- %s seconds ---" % (time.time() - start_time))


    t2 = ["a" + str(i) for i in t]
    start_time = time.time()
    for i in t2:
        a[i]
    
    print("--- %s seconds ---" % (time.time() - start_time))

