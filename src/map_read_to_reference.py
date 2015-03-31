from __future__ import print_function
from __future__ import division


LOG = None

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

import map_read_to_reference

all_perc = map_read_to_reference.process_all_dirs(".", 150)


"""

def process_all_dirs(dir, cutoff_len=200, mask_lower_case=False, print_log=True):
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
        percs = process_all(dir, cutoff_len, mask_lower_case=mask_lower_case)
        if len(percs) > 0:
            #fq_perc.update(percs)
            fq_perc[os.path.basename(dir)] = percs
        
    export_perc(all_perc=fq_perc)
        
    if print_log:
        LOG.close()
        
    return fq_perc



def process_all_fq(fq_dir, cutoff_len, mask_lower_case=False, print_log=True):
    import glob
    import os
    global LOG
    
    if print_log:
        LOG = open("all_perc.log", "w")
    else:
        LOG = None
    
    fq_fns = glob.glob(fq_dir + "/*.fq")
    fq_perc = {}
    for fq_fn in fq_fns:
        print("Processing " + fq_fn)
        percs = get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=mask_lower_case)
        fq_perc[fq_fn.replace(".fq", "")] = percs
    
    
    
    if print_log:
        LOG.close()
        
    return fq_perc


"""
import map_read_to_reference
map_read_to_reference.export_perc()


"""
def expprt_perc2(all_perc, out_fn="all_perc.stat"):
    sample_ids = list(all_perc.keys())
    with open(out_fn, "w") as OUT:
        OUT.write("Species\t" + "\t".join(sample_ids) + "\n")
        for species_id in all_perc.keys():
            line = [species_id]
            perc = all_perc[species_id]
            for sample_id in sample_ids:
                if sample_id in perc.keys():
                    line.append(str(perc[sample_id]))
                else:
                    line.append("0.0")
            OUT.write("\t".join(line) + "\n")


def export_perc(out_fn="all_perc.stat", all_perc=None, dir=".", mask_lower_case=False, cutoff_len=150, sample_ids = ["GZ-Cell_Y2", "GZ-Cell_Y1", "GZ", "GZ-Xyl_Y1", "GZ-Xyl_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]):
    import map_read_to_reference

    if all_perc is None:
        print("Export_perc will extract info from the current directory.")
        all_perc = map_read_to_reference.process_all_dirs(".", 150)
        
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
import map_read_to_reference

all_perc = map_read_to_reference.process_all("Clostridium_thermocellum", 150)

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
#         fq_perc[id] = get_perc_cutoff(fq_fn, cutoff_len)
#         
#     return fq_perc


def process_all(dir, cutoff_len, mask_lower_case=False):
    import glob
    import os
    
    excluding_ext = "-" + os.path.basename(dir) + ".fq"
    
    fq_perc = {}
    fq_fns = glob.glob(dir + "/*.fq")
    for fq_fn in fq_fns:
        id = os.path.basename(fq_fn)
        id = id.replace(excluding_ext, "")
        if (os.stat(fq_fn)).st_size < 100:
            fq_perc[id] = 0.0
            continue
        fq_perc[id] = get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=mask_lower_case)
        
    return fq_perc




def get_perc_cutoff(fq_fn, cutoff_len, mask_lower_case=False):
    from Bio import SeqIO
    import re
    
    seqs = SeqIO.index(fq_fn, "fastq")
    seqs = [str(seqs[sid].seq) for sid in seqs.keys()]
    seq_len = sum([len(s) for s in seqs])
      
    #print("Seq len from " + fq_fn + ": " + str(seq_len))  
                
    lens = estimate_fq_coverage(fq_fn, mask_lower_case=mask_lower_case)
    
    if len(lens) == 0:
        return -1.0
    
    #print("Lens=" + str(len(lens)))
    
    sum_len = sum([int(lens[l]) * int(l) for l in lens.keys() if int(l) >= cutoff_len])
    print_log("Length of sequence from " + fq_fn + ": " + str(seq_len) + ", Lens=" + str(len(lens)) + ", Total sum=" + str(sum_len) + ", Perc=" + str((sum_len / seq_len) * 100))
    
    return (sum_len / seq_len) * 100



def print_log(msg):
    global LOG
    print(msg)
    if LOG is not None:
        LOG.write(msg + "\n")



def estimate_fq_coverage(fq_fn, mask_lower_case=False):
    from Bio import SeqIO
    import re
    
    seqs = SeqIO.index(fq_fn, "fastq")
    seqs = [str(seqs[sid].seq) for sid in seqs.keys()]
    
    lens_count = {}
    for seq in seqs:
        #if mask_lower_case:
        #    seq = "".join([s for s in seq if s.isupper() or s == 'n'])
        #lens = re.split("n+", seq)
        if mask_lower_case:
            lens = re.split("[a-z]+", seq)
        else:
            lens = re.split("n+", seq)
            re
        for l in lens:
            l = len(l)
            if l == 0:
                continue
            
            if l not in lens_count:
                lens_count[l] = 0
            lens_count[l] = lens_count[l] + 1
    return lens_count


"""
import map_read_to_reference

bla_fn = "blast_to_ref/" + "GZ-Cell_Y2.001+all_prokaryotes+BacterialDB.fasta.bla"
selected_sid_tbl = map_read_to_reference.extract_sid_tax(bla_fn)




lens = map_read_to_reference.compute_sid_match_len(bla_fn)
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
def extract_sid_tax(bla_fn, cutoff=40000, tax_fn="/home/siukinng/db/BioProject_Prokaryotes/prokaryotes.txt"):
    lens = compute_sid_match_len(bla_fn)
    sorted_sids = sorted(lens.items(), key=lambda x:x[1], reverse=True)

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
    selected_sid_tbl = {ss[0]:sid_tbl[ss[0]] for ss in sorted_sids if ss[1] > cutoff and ss[0] in sid_tbl.keys()}
    
    for sid in selected_sid_tbl.keys():
        selected_sid_tbl[sid].append(lens[sid])
    
    return selected_sid_tbl
   

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
        
