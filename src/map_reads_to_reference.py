from __future__ import print_function
from __future__ import division

import glob
import os
import numpy
import sys
from Bio import SeqIO


sys.path.append(os.path.abspath("/home/siukinng/tools/scripts"))


import mg_pipeline
import pick_seq
import filter_blast_res
import process_bam
import process_binning


LOG = None


"""

Given a directory which contains files with file extension as specified by fa_fn_ext, this function will map 
the fasta sequences of the files with fa_fn_ext extension to the NCBI BioProject DB. The identified reference
genomes will then serve as a template for mapping reads from fq files with a same file name (XXX_1 and XXX_2) from
the directory as given in fq_dir.
 
"""
def map_all_fa(wd_dir=".", fa_fn_ext="fa", fq_dir=".", fq_fn_ext="fq", ref_fa_fn=None, generate_ref_for_all_fa_fn=False):
    # Extract all fasta sequences from the working directory (wk_dir)
    fa_fns = glob.glob(wk_dir + "/*." + fa_fn_ext)
    
    # Exit if no file with the specified file extension (fa_fn_ext) is found
    if len(fa_fns) == 0:
        mg_pipeline.print_status("No fasta file is found at " + wk_dir)
        return 0
    
    # Check if the sequence of ref genomes is provided.
    # If not, we will map the 
    if ref_fa_fn is None:
        # map all fa_fns to the reference genomes
        sids = map_all_fa_fns_to_references(fa_fns)
        
        ref_fa_fn = wd_dir + "/combined_ref.fasta"
        
        # generate the ref genome file
        status = generate_ref_genomes(sids, ref_fa_fn)
    
        # 
        if len([k for k in status.keys() if status[k] != 0]) == 0:
            mg_pipeline.print_status("Problem on preparing reference genomes, aborted.")
            return None
    else:
        if len(os.path.dirname(ref_fa_fn)) == 0:
            ref_fa_fn = wd_dir + "/" + ref_fa_fn
    
    for fa_fn in fa_fns:
        sample_id = (os.path.basename(fa_fn))[::-1].split(".",1)[1][::-1]
        mg_pipeline.print_status("Processing " + sample_id + " (" + fa_fn + ")")
        
        mg_pipeline.print_status("Checking fq files: ")
        read_1_fn = fq_dir + "/" + sample_id + "_1.fq"
        read_2_fn = fq_dir + "/" + sample_id + "_2.fq"
        if not mg_pipeline.assert_proc(read_1_fn) or not mg_pipeline.assert_proc(read_1_fn):
            mg_pipeline.print_status(read_1_fn + " or " + read_2_fn + " is not found, sample skipped.")
            continue
        else:
            mg_pipeline.print_status(read_1_fn + " and " + read_2_fn + " will be used.")
                
        [sorted_bam_fn, fq_fn, cov_fn] = map_reads_to_reference(sample_id, read_1_fn, read_2_fn, fa_fn, ref_fa_fn)
    
    


def map_reads_to_reference(sample_id, read_1_fq, read_2_fq, contig_fa_fn, ref_fa_fn=None, ref_ids=None, subject_id=None, out_dir="."):
    mg_pipeline.print_status("Initializing mapping process...")
    
    if ref_fa_fn is None:
        if subject_id is not None:
            ref_fa_fn = out_dir + "/" + subject_id + ".fa"
            prefix = sample_id + "+" + subject_id
        else:
            ref_fa_fn = out_dir + "/" + sample_id + ".ref.fa"
            prefix = sample_id + "+" + "ref"
        
        ref_nr_sids = ref_ids
        if ref_nr_sids is None:
            ref_nr_sids = mapping_contigs_to_references(contig_fa_fn)
        status = process_binning.generate_ref_genomes(ref_nr_sids, ref_fa_fn)
        
        if len([k for k in status.keys() if status[k] != 0]) == 0:
            mg_pipeline.print_status("Problem on preparing reference genomes.")
    else:
        prefix = sample_id + "+" + os.path.splitext((os.path.basename(ref_fa_fn)))[0]
        
    sorted_bam_fn = mg_pipeline.run_bwa(prefix, read_1_fq, read_2_fq, ref_fa_fn)
    
    if sorted_bam_fn is None:
        mg_pipeline.print_status("Problem on generating a sorted bam file.")

    
    fq_fn = prefix + ".fq"
    if not mg_pipeline.assert_proc(fq_fn):
        if mg_pipeline.generate_fq_from_bam(sorted_bam_fn, ref_fa_fn, fq_fn) is None:
            mg_pipeline.print_status("Problem on preparing fq file.") 
    else:
        mg_pipeline.print_status(fq_fn + " exists.")
        
    cov_fn = mg_pipeline.generate_coverage_from_bam(sorted_bam_fn) 
    if cov_fn is None:
        mg_pipeline.print_status("Problem on generating a coverage file.")
    
    return [sorted_bam_fn, fq_fn, cov_fn]



def map_read_to_NR(query_fa_fn, db_fn="nr_idx", diamond_program="blastx", max_target_seqs=1, num_threads=40, tmp_dir=".", diamond_path=None):
    
    # ~/tools/megan/diamond blastx --max-target-seqs 1 -p 40 -d nr_idx -q ./$id.fa -a $id+nr -t /scratch/bcnskaa/tmp
    print("")



def assign_taxa_to_reads(m8_fn, prot2gi_fn, gi2tax_fn):
    print("")



"""
Return a list of reads that do not map to any reference sequence
Ref: https://www.biostars.org/p/45654/


samtools view -f4 whole.bam | cut -f10 | sort | uniq -c | sort -nr > unmapped_unique.count
"""
def extract_unmapped_read_ids(bam_fn, out_fn=None, samtool_path=mg_pipeline.SAMTOOLS_HOME):
        
    cmd = samtool_path + "/samtools view -f " + bam_fn + " | cut -f10 | sort | uniq > " + out_fn 
    
    


def extract_multi_mapped_read_ids(bam_fn, samtool_path=mg_pipeline.SAMTOOLS_HOME):
    cmd = samtool_path + "/samtools" 
    

"""
Return a non-redundant list of subject ids mapped by query sequences
"""
def get_all_ref_sids_from_bla_fns(bla_fns):
    sids = []
    for bla_fn in bla_fns:
        print("Processing " + bla_fn)
        ids = filter_blast_res.get_nr_ids_from_fn(bla_fn, is_query_id_selected=False)
        if len(ids) > 0:
            sids = sids + ids
    
    return list(set(sids))

    

"""

Mapping contig sequences to reference genomes available in NCBI BioProject DB, and return a list of

"""
def map_all_fa_fns_to_references(*contig_fa_fns):
    sids = []
    for contig_fa_fn in contig_fa_fns:
        ids = mapping_contigs_to_references(contig_fa_fn)
        if len(ids) > 0:
            sids = sids + ids
    
    return list(set(sids))



"""

BLAST the contigs against NCBI BioProjects and determine the closet reference genome to which contig sequences can map.

"""
def mapping_contigs_to_references(contig_fn, out_fn=None, out_dir=None, len_cutoff=3000, db_fn=mg_pipeline.NCBI_BIOPROJECT_DB):
    print_msg("Mapping contigs from " + contig_fn + " to " + db_fn)
    
    subject_id = ((os.path.basename(db_fn)[::-1].split(".", 1))[1])[::-1]
    query_id = ((os.path.basename(contig_fn)[::-1].split(".", 1))[1])[::-1]
    
    if out_fn is None:
        if out_dir is None:
            #out_fn = os.path.dirname(contig_fn) + "/" + query_id
            out_fn = "./" + query_id
        else:
            out_fn = out_dir + "/" + query_id
        out_fn = out_fn + "+" + subject_id + ".bla" 
    else:
        out_dir = os.path.dirname(out_fn)
        if len(out_dir) == 0:
            if out_dir is None:
                out_fn = "./" + out_fn
            else:
                out_fn = out_dir + "/" + out_fn
        

    print_msg("Result will be exported to " + out_fn)
    
    out_fn = mg_pipeline.blastn(contig_fn, db_fn, outdir=out_dir, outfn=out_fn, perc_identity=85, max_target_seqs=1, evalue=1e-30)
    
    blast_res = filter_blast_res.import_blast_res(out_fn, thres_aln_len = len_cutoff)
    
    nr_sids = filter_blast_res.get_nr_ids(blast_res, is_query_id_selected=False)
    
    return nr_sids





"""

sample_id="GZ-Cell_Y2"
ref_id="Oscillibacter_valericigenes_Sjm18-20"
ref_seq_fn="$ref_id.fasta"
read_1_fn="~/samples/lab/GZ-Cell_Y2/GZ-Cell_Y2_1.trimmed.cleaned.fq"
read_2_fn="~/samples/lab/GZ-Cell_Y2/GZ-Cell_Y2_2.trimmed.cleaned.fq"

cmd="~/tools/bwa/bwa index $ref_seq_fn"
echo $cmd
eval $cmd
cmd="~/tools/bwa/bwa mem -t 14 -k 51 $ref_seq_fn $read_1_fn $read_2_fn > $sample_id.sam"
echo $cmd
eval $cmd
cmd="~/tools/samtools/samtools view -Sb $sample_id.sam > $sample_id.bam"
echo $cmd
eval $cmd
cmd="~/tools/samtools/samtools sort $sample_id.bam $sample_id.sorted"
echo $cmd
eval $cmd
cmd="~/tools/samtools/samtools mpileup -uf $ref_seq_fn $sample_id.sorted.bam | ~/tools/samtools/bcftools/bcftools view -cg - | ~/tools/samtools/bcftools/vcfutils.pl vcf2fq  > $sample_id-$ref_id.fq"
echo $cmd
eval $cmd


"""

"""

import map_reads_to_reference

cutoff_len = 800
all_perc = map_reads_to_reference.process_all_dirs(".", cutoff_len)
sample_ids = ["GZ-Cell_Y2", "GZ-Cell_Y1", "GZ-Seed_Y0", "GZ-Xyl_Y1", "GZ-Xyl_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]

map_reads_to_reference.expprt_perc2(all_perc, "all_perc.l" + str(cutoff_len) + "_wolowermask.stat", sample_ids)

"""

def process_all_dirs(dir, cutoff_len=200, mask_lower_case=False, print_log=True, delimiter="+"):
    import glob
    import os
    
    if print_log:
        LOG = open("all_perc.log", "w")
    else:
        LOG = None
    
    dirs = glob.glob(dir + "/*")
    dirs = [d for d in dirs if os.path.isdir(d)]
    fq_perc = {}
    for dir in dirs:
        percs = process_all(dir, cutoff_len, mask_lower_case=mask_lower_case, delimiter=delimiter)
        if len(percs) > 0:
            #fq_perc.update(percs)
            fq_perc[os.path.basename(dir)] = percs
        
    export_perc(all_perc=fq_perc)
        
    if print_log:
        LOG.close()
        
    return fq_perc



def process_all_fq(fq_dir, cutoff_len, mask_lower_case=False, print_log=True, fq_fn_ext="fq", delimiter="+", washed_ids=False):
    import glob
    import os
    global LOG
    
    if print_log:
        LOG = open("all_perc.log", "w")
    else:
        LOG = None
    
    fq_fns = glob.glob(fq_dir + "/*." + fq_fn_ext)
    fq_perc = process_selected_fq(fq_fns, cutoff_len, mask_lower_case=mask_lower_case, print_log=print_log, fq_fn_ext=fq_fn_ext, delimiter=delimiter, washed_ids=washed_ids)
#     fq_perc = {}
#     for fq_fn in fq_fns:
#         print("Processing " + fq_fn)
#         percs = get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=mask_lower_case)
#         #fq_perc[fq_fn.replace("."+fq_fn_ext, "")] = percs
#         fq_perc[fq_fn.replace("."+fq_fn_ext, "")] = sum([percs[k] for k in percs.keys()])
#         
#         
#     if print_log:
#         LOG.close()
        
    return fq_perc



"""
Provide a list of fq files, 
"""
def process_selected_fq(fq_fns, cutoff_len, combined_seqs_from_same_fn=False, mask_lower_case=False, print_log=True, delimiter="+", washed_ids=True):
    import glob
    import os
    global LOG
    
    if delimiter is None or len(delimiter) == 0:
        delimiter = "#####"
    
    if print_log:
        LOG = open("all_perc.log", "w")
    else:
        LOG = None

    if len(fq_fns) == 0:
        return {}

    fq_perc = {}
    for fq_fn in fq_fns:
        #sample_id = ((os.path.basename(fq_fn)[::-1].split(".", 1))[1])[::-1]
        sample_id = os.path.basename(fq_fn).split(delimiter)[0]
        
        print("Processing " + sample_id + " (" + fq_fn+ ")")
        
        percs = get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=mask_lower_case, washed_ids=washed_ids)
        
        
        print("Percs_size="+str(len(percs)))
        
        #fq_perc[fq_fn.replace("."+fq_fn_ext, "")] = percs
        #fq_perc[fq_fn.replace("."+fq_fn_ext, "")] = sum([percs[k] for k in percs.keys()])
        if combined_seqs_from_same_fn:
            fq_perc[sample_id] = sum([percs[k] for k in percs.keys()])
        else:
            #fq_perc[sample_id] = percs
            for species_id in percs.keys():
                if species_id not in fq_perc.keys():
                    fq_perc[species_id] = {}
                fq_perc[species_id][sample_id] = percs[species_id]
        
    if print_log:
        LOG.close()
        
    return fq_perc 


"""

Given the ids (nr_sids), this function will extract the corresponding sequences from the
NCBI BioProject DB and export to the out file (out_fn). If the reference genome consists of 
many contigs, they will merge into a single sequence.

"""
def generate_ref_genomes(nr_sids, out_fn, db_dir=mg_pipeline.NCBI_BIOPROJECT_DB_HOME):
    print_msg("Preparing " + str(len(nr_sids)) +" reference sequences from " + db_dir)
     
    nr_sids = wash_id(nr_sids)
     
    OUT = open(out_fn, "w")
    status = {sid:0 for sid in nr_sids}
    for sid in nr_sids:
        sid_fn = db_dir + "/" + sid + ".fna"
        s = 0
        if os.path.exists(sid_fn):
            seqs = SeqIO.index(sid_fn, "fasta")
            if len(seqs) == 1:
                for seq in SeqIO.parse(sid_fn, "fasta"):
                    SeqIO.write(seq, OUT, "fasta")
            else:
                [exported_n, seq] = merge_fa.merge(sid_fn, sid)
                merge_fa.export_fasta(seq[0], seq[1], OUT)
            s = s + 1
        status[sid] = s
    OUT.close()
    
    print_msg(str(len([s for s in status.keys() if status[s] > 0])) + " out of " + str(len(status)) + " were found.")
    
    return status
    




"""

import map_reads_to_reference
import glob

fq_fns = glob.glob("*.fq")
fq_fns = fq_fns[4:6]

cutoff_len = 800
percs = map_reads_to_reference.process_selected_fq(fq_fns, cutoff_len)
map_reads_to_reference.export_perc2(percs, "all_perc."+str(cutoff_len)+".stat", export_tax_info=True)


"""
def export_perc2(all_perc, out_fn="all_perc.stat", sample_ids=None, export_tax_info=False):
    print("perc_len=" + str(len(all_perc)))
    
    species_ids = list(all_perc.keys())
    if sample_ids is None:
         sample_ids = []
         for species_id in species_ids:
             sample_ids = sample_ids + list(all_perc[species_id].keys())
         sample_ids = list(set(sample_ids))
    
    with open(out_fn, "w") as OUT:
        if export_tax_info:
            OUT.write("Code\tPhylum\tClass\tSpecies\t" + "\t".join(sample_ids) + "\n")
        else:
            OUT.write("Species\t" + "\t".join(sample_ids) + "\n")
        
        for species_id in species_ids:
            line = [species_id]
            
            # Assign tax info
            if export_tax_info:
                tax = get_tax(species_id)
                if tax is None:
                    tax = "\t"
                line.extend(tax) 
                          
            perc = all_perc[species_id]
            for sample_id in sample_ids:
                if sample_id in perc.keys():
                    line.append(str(perc[sample_id]))
                else:
                    line.append("0.0")
            OUT.write("\t".join(line) + "\n")



def export_perc(out_fn="all_perc.stat", all_perc=None, dir=".", mask_lower_case=False, cutoff_len=150, sample_ids = ["GZ-Cell_Y2", "GZ-Cell_Y1", "GZ-Seed_Y0", "GZ-Xyl_Y1", "GZ-Xyl_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]):
    import map_reads_to_reference

    if all_perc is None:
        print("Export_perc will extract info from the current directory.")
        all_perc = map_reads_to_reference.process_all_dirs(".", cutoff_len, mask_lower_case)
        
    #sample_ids = ["GZ-Cell_Y2", "GZ-Cell_Y1", "GZ", "GZ-Xyl_Y1", "GZ-Xyl_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    
    with open(out_fn, "w") as OUT:
        OUT.write("Species\t" + "\t".join(sample_ids) + "\n")
        for species_id in all_perc.keys():
            line = [species_id]
            perc = all_perc[species_id]
            for sample_id in sample_ids:
                if sample_id in perc.keys():
                    line.append(str(perc[sample_id]))
                else:
                    print(sample_id + " is not found in the " + species_id +" perc list.")
                    line.append("0.0")
            OUT.write("\t".join(line) + "\n")




"""
# Blast all bin fasta to reference genomes
len="2000"

echo "" > batch_cmd
sample_ids=( SWH-Cell_Y1 SWH-Cell_Y2 GZ-Cell_Y1 GZ-Cell_Y2 SWH-Xyl_Y2 SWH-Xyl_Y1 GZ-Xyl_Y2 GZ-Xyl_Y1 GZ-Seed_Y0 SWH-Seed_Y0 SWH-Cell55_Y2 )
for f in ~/MG/scaffolds_5000/mapping_to_references/*/*.fasta;do
        id=${f##*/}
        id=${id/.fasta/}
        if [ -d "./$id" ];then
                echo $id
                mkdir $id
                db_fn="./$id/$id.fasta"
                ln -s $f $db_fn
                ~/tools/blast/bin/makeblastdb -in $db_fn -dbtype nucl
                for sample_id in ${sample_ids[*]};do
                        sample_dir="/home/siukinng/MG/scaffolds_$len/"$sample_id"_"$len
                        for bin_fn in $sample_dir/*.fasta;do
                                bin_id=${bin_fn##*/}
                                bin_id=${bin_id/.fasta/}
                                cmd="~/tools/blast/bin/blastn -query $bin_fn -db $db_fn -outfmt 6 -perc_identity 89 -max_target_seqs 1 -evalue 1e-50 -out ./$id/$bin_id-$id.bla -num_threads 16"
                                echo $cmd >> batch_cmd
                        done
                done
        fi
done
"""




"""
import map_reads_to_reference

all_perc = map_reads_to_reference.process_all("Clostridium_thermocellum", 150)

"""
# def process_all(dir, cutoff_len):
#     import glob
#     import os
#     
#     excluding_ext = "-" + os.path.basename(dir) + ".fq"
#     
#     fq_perc = {}
#     fq_fns = glob.glob(dir + "/*.fq")
#     for fq_fn in fq_fns:
#         id = os.path.basename(fq_fn)
#         id = id.replace(excluding_ext, "")
#         if (os.stat(fq_fn)).st_size < 100:
#             fq_perc[id] = 0.0
#             continue
#         #fq_perc[id] = get_perc_cutoff(fq_fn, cutoff_len)
#         percs = get_perc_cutoff(fq_fn, cutoff_len)
#         fq_perc[id] = sum([percs[k] for k in percs.keys()])
#     return fq_perc

"""
Dir represents a tax 
"""
def process_all(dir, cutoff_len, mask_lower_case=False, fq_fn_ext="fq", delimiter="+"):
    import glob
    import os
    
    excluding_ext = delimiter + os.path.basename(dir) + "." + fq_fn_ext
    
    fq_perc = {}
    fq_fns = glob.glob(dir + "/*." + fq_fn_ext)
    for fq_fn in fq_fns:
        id = (os.path.basename(fq_fn)).split(delimiter)[0]
        #id = (os.path.basename(fq_fn)).replace(excluding_ext, "")
        #id = id.replace(excluding_ext, "")
        if (os.stat(fq_fn)).st_size < 100:
            fq_perc[id] = 0.0
            continue
        #fq_perc[id] = get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=mask_lower_case)
        percs = get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=mask_lower_case)
        fq_perc[id] = sum([percs[k] for k in percs.keys()])
        
    return fq_perc


"""

"""
#def process_fq(fq, cutoff_len, mask_lower_case=False, delimiter="+"):


def get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=False, washed_ids=False):
#def get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=False, combined_all_seq=True, delimiter="+", washed_ids=False):
    from Bio import SeqIO
    import re
    
    seqs = SeqIO.index(fq_fn, "fastq")
    
    if washed_ids:
        seqs = {get_ref_id(sid): str(seqs[sid].seq) for sid in seqs.keys()}
    else:
        seqs = {sid: str(seqs[sid].seq) for sid in seqs.keys()}
        
    percs = {}
    for seq_id in seqs.keys():
        seq_len = len(str(seqs[seq_id]))
        lens = estimate_coverage(seqs[seq_id],  mask_lower_case=mask_lower_case)
        
        if len(lens) == 0:
            percs[seq_id] -1.0    
    
        sum_len = sum([int(lens[l]) * int(l) for l in lens.keys() if int(l) >= cutoff_len])
        print_log("Length of sequence " + seq_id + " from " + fq_fn + ": " + str(seq_len) + ", Lens=" + str(len(lens)) + ", Total sum=" + str(sum_len) + ", Perc=" + str((sum_len / seq_len) * 100))
    
        percs[seq_id] = (sum_len / seq_len) * 100            

    return percs


# 
#      
#     seqs = [str(seqs[sid].seq) for sid in seqs.keys()]
#     seq_len = sum([len(s) for s in seqs])
#       
#     #print("Seq len from " + fq_fn + ": " + str(seq_len))  
#                 
#     lens = estimate_fq_coverage(fq_fn, mask_lower_case=mask_lower_case)
#     
#     if len(lens) == 0:
#         return -1.0
#     
#     #print("Lens=" + str(len(lens)))
#     
#     sum_len = sum([int(lens[l]) * int(l) for l in lens.keys() if int(l) >= cutoff_len])
#     print_log("Length of sequence from " + fq_fn + ": " + str(seq_len) + ", Lens=" + str(len(lens)) + ", Total sum=" + str(sum_len) + ", Perc=" + str((sum_len / seq_len) * 100))
#     
#     return (sum_len / seq_len) * 100


# def get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=False, combined_all_seq=True):
#     from Bio import SeqIO
#     import re
#     
#     
#     seqs = SeqIO.index(fq_fn, "fastq")    
#     seqs = [str(seqs[sid].seq) for sid in seqs.keys()]
#     seq_len = sum([len(s) for s in seqs])
#       
#     #print("Seq len from " + fq_fn + ": " + str(seq_len))  
#                 
#     lens = estimate_fq_coverage(fq_fn, mask_lower_case=mask_lower_case)
#     
#     if len(lens) == 0:
#         return -1.0
#     
#     #print("Lens=" + str(len(lens)))
#     
#     sum_len = sum([int(lens[l]) * int(l) for l in lens.keys() if int(l) >= cutoff_len])
#     print_log("Length of sequence from " + fq_fn + ": " + str(seq_len) + ", Lens=" + str(len(lens)) + ", Total sum=" + str(sum_len) + ", Perc=" + str((sum_len / seq_len) * 100))
#     
#     return (sum_len / seq_len) * 100



def print_log(msg):
    global LOG
    print(msg)
    if LOG is not None:
        LOG.write(msg + "\n")




"""

!!! 

"""
def estimate_coverage(seq, mask_lower_case=False):
    import re
    
    lens_count = {}
    
    if mask_lower_case:
        lens = re.split("[a-z]+", seq)
    else:
        lens = re.split("n+", seq)
        re
    for l in lens:
        l = len(l)
        if l == 0:
            continue
        
        if l not in lens_count.keys():
            lens_count[l] = 0
        lens_count[l] = lens_count[l] + 1
    return lens_count



"""

"""
def estimate_fq_coverage(fq_fn, mask_lower_case=False):
    from Bio import SeqIO
    
    seqs = SeqIO.index(fq_fn, "fastq")
    seqs = [str(seqs[sid].seq) for sid in seqs.keys()]
    
    lens_count = {}
    for seq in seqs:
        lens = estimate_coverage(seq, mask_lower_case)
        for l in lens:
            l = len(l)
            if l == 0:
                continue
            
            if l not in lens_count.keys():
                lens_count[l] = 0
            lens_count[l] = lens_count[l] + 1
    return lens_count



    
# def estimate_fq_coverage(fq_fn, mask_lower_case=False):
#     from Bio import SeqIO
#     import re
#     
#     seqs = SeqIO.index(fq_fn, "fastq")
#     seqs = [str(seqs[sid].seq) for sid in seqs.keys()]
#     
#     lens_count = {}
#     for seq in seqs:
#         #if mask_lower_case:
#         #    seq = "".join([s for s in seq if s.isupper() or s == 'n'])
#         #lens = re.split("n+", seq)
#         if mask_lower_case:
#             lens = re.split("[a-z]+", seq)
#         else:
#             lens = re.split("n+", seq)
#             re
#         for l in lens:
#             l = len(l)
#             if l == 0:
#                 continue
#             
#             if l not in lens_count:
#                 lens_count[l] = 0
#             lens_count[l] = lens_count[l] + 1
#     return lens_count



"""
import map_reads_to_reference

bla_fn = "blast_to_ref/" + "GZ-Cell_Y2.001+all_prokaryotes+BacterialDB.fasta.bla"
selected_sid_tbl = map_reads_to_reference.extract_tax_from_sid(bla_fn)




lens = map_reads_to_reference.compute_sid_match_len(bla_fn)
sorted_sids = sorted(lens.items(), key=lambda x:x[1], reverse=True)


lens[sorted_sids[0]]
tax_fn = "/home/siukinng/db/BioProject_Prokaryotes/prokaryotes.txt"
with open(tax_fn) as IN:
    tax = IN.read().splitlines()
tax1 = {t.split("\t")[12]:[t.split("\t")[4], t.split("\t")[5], t.split("\t")[0]] for t in tax}
tax2 = {t.split("\t")[8]:[t.split("\t")[4], t.split("\t")[5], t.split("\t")[0]] for t in tax}
tax1.update(tax2)

sid_tbl = {}
for sid in sorted_sids:
    sid = sid[0]
    if sid in tax1.keys():
        sid_tbl[sid] = tax1[sid]

selected_sid_tbl = {ss[0]:sid_tbl[ss[0]] for ss in sorted_sids if ss[1] > 40000 and ss[0] in sid_tbl.keys()}    

"""
# def extract_tax_from_sid(bla_fn, cutoff=40000, tax_fn="/home/siukinng/db/BioProject_Prokaryotes/prokaryotes.txt"):
#     lens = compute_sid_match_len(bla_fn)
#     sorted_sids = sorted(lens.items(), key=lambda x:x[1], reverse=True)
# 
#     with open(tax_fn) as IN:
#         tax = IN.read().splitlines()
#     tax1 = {t.split("\t")[12]:[t.split("\t")[4], t.split("\t")[5], t.split("\t")[0]] for t in tax}
#     tax2 = {t.split("\t")[8]:[t.split("\t")[4], t.split("\t")[5], t.split("\t")[0]] for t in tax}
#     tax1.update(tax2)
#     
#     sid_tbl = {}
#     for sid in sorted_sids:
#         sid = sid[0]
#         if sid in tax1.keys():
#             sid_tbl[sid] = tax1[sid]
#     selected_sid_tbl = {ss[0]:sid_tbl[ss[0]] for ss in sorted_sids if ss[1] > cutoff and ss[0] in sid_tbl.keys()}
#     
#     for sid in selected_sid_tbl.keys():
#         selected_sid_tbl[sid].append(lens[sid])
#     
#     return selected_sid_tbl
def extract_tax_from_sid(bla_fn, cutoff=30000, tax_fn="/home/siukinng/db/BioProject_Prokaryotes/prokaryotes.txt"):
    lens = compute_sid_match_len(bla_fn)
    sorted_sids = sorted(lens.items(), key=lambda x:x[1], reverse=True)

    sid_tbl = {}
    for sid in sorted_sids:
        sid = sid[0]
        
        tax_id = get_tax(sid)
        if tax_id is not None:
            sid_tbl[sid] = tax_id
            
        #if sid in tax1.keys():
        #    sid_tbl[sid] = tax1[sid]
    selected_sid_tbl = {ss[0]:sid_tbl[ss[0]] for ss in sorted_sids if ss[1] > cutoff and ss[0] in sid_tbl.keys()}
    
    #for sid in selected_sid_tbl.keys():
    #    selected_sid_tbl[sid].append(lens[sid])
    
    return selected_sid_tbl



"""

"""
def get_tax_from_bla_fns(bla_fns, cutoff=30000, tax_fn="/home/siukinng/db/BioProject_Prokaryotes/prokaryotes.txt"):
    combined_sid_tbl = {}
    for bla_fn in bla_fns:
        print("Assigning tax info to " + bla_fn)
        sid_tbl = extract_tax_from_sid(bla_fn, cutoff=cutoff, tax_fn=tax_fn)
        combined_sid_tbl.update(sid_tbl)
    return combined_sid_tbl


"""

"""
tax_db = None
def get_tax(sid, tax_fn="/home/siukinng/db/BioProject_Prokaryotes/prokaryotes.txt"):
    global tax_db
    
    if tax_db is None:
        print("Preparing tax_db")
        with open(tax_fn) as IN:
            tax = IN.read().splitlines()
        tax1 = {t.split("\t")[12]:[t.split("\t")[4], t.split("\t")[5], t.split("\t")[0]] for t in tax}
        tax2 = {t.split("\t")[8]:[t.split("\t")[4], t.split("\t")[5], t.split("\t")[0]] for t in tax}
        tax1.update(tax2)
        tax_db = tax1
    
    if sid in tax_db.keys():
        #print(sid+"="+",".join(tax_db[sid]))
        return tax_db[sid]
    else:
        print("No tax info available for " + sid)
        #return None
        return ["","",""]



"""

"""
def extract_species_from_tax(sid_tbl):
    species = []
    for sid in sid_tbl.keys():
        species.append(" ".join(sid_tbl[sid][2].split(" ")[0:2]))
    return species



"""

"""
def assign_contigs_to_sid(bla_fn, contig_fa_fn, contig2sid_tbl_fn=None, cutoff=40000, blast_len_cutoff=3000):
    print("")



"""

"""
def compute_sid_match_len(bla_fn):
    with open(bla_fn) as IN:
        bla = IN.read().splitlines()
    bla = [b.split("\t") for b in bla]
    
    sid_match_lens = {}
    for b in bla:
        sid = b[1]
        len = int(b[3])
        if sid not in sid_match_lens.keys():
            sid_match_lens[sid] = 0
        sid_match_lens[sid] = sid_match_lens[sid] + len
        
    return sid_match_lens
        


"""
If sequence id is in format "XXXX.X", this function extracts the text before "."
"""
def wash_ids(ids):
    washed_ids = [id.split(".")[0] for id in ids]
    return washed_ids


def wash_id(id):
    return id.split(".")[0]




"""
If the id is NCBI compliant format, we can extract the NC_ number
"""
def get_ref_id(sid):
    #washed_sids = []
    #for sid in sids:
    if "|" in sid:
        return sid.split("|")[3]
    #    washed_sids.append(sid.split("|")[3])
    else:
        return sid
    #    washed_sids.append(sid)
    #return washed_sids