from __future__ import print_function
from __future__ import division

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

def process_all_dirs(dir, cutoff_len):
    import glob
    import os
    
    dirs = glob.glob(dir + "/*")
    dirs = [d for d in dirs if os.path.isdir(d)]
    fq_perc = {}
    for dir in dirs:
        percs = process_all(dir, cutoff_len)
        if len(percs) > 0:
            #fq_perc.update(percs)
            fq_perc[os.path.basename(dir)] = percs
    return fq_perc



def export_perc(out_fn="all_perc.stat"):
    import map_read_to_reference

    all_perc = map_read_to_reference.process_all_dirs(".", 150)
    sample_ids = ["GZ-Cell_Y2", "GZ-Cell_Y1", "GZ-Seed_Y0", "GZ-Xyl_Y1", "GZ-Xyl_Y2", "SWH-Xyl_Y2", "SWH-Xyl_Y1", "SWH-Seed_Y0", "SWH-Cell_Y1", "SWH-Cell_Y2", "SWH-Cell55_Y2"]
    
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
            OUT.write("\t".join(line))


"""
import map_read_to_reference

all_perc = map_read_to_reference.process_all("Clostridium_thermocellum", 150)

"""
def process_all(dir, cutoff_len):
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
        fq_perc[id] = get_perc_cutoff(fq_fn, cutoff_len)
        
    return fq_perc



def get_perc_cutoff(fq_fn, cutoff_len):
    from Bio import SeqIO
    import re
    
    seqs = SeqIO.index(fq_fn, "fastq")
    seqs = [str(seqs[sid].seq) for sid in seqs.keys()]
    seq_len = sum([len(s) for s in seqs])
      
    #print("Seq len from " + fq_fn + ": " + str(seq_len))  
                
    lens = estimate_fq_coverage(fq_fn)
    #print("Lens=" + str(len(lens)))
    
    sum_len = sum([int(lens[l]) * int(l) for l in lens.keys() if int(l) > cutoff_len])
    print("Length of sequence from " + fq_fn + ": " + str(seq_len) + ", Lens=" + str(len(lens)) + ", Total sum=" + str(sum_len) + ", Perc=" + str((sum_len / seq_len) * 100))
    
    return (sum_len / seq_len) * 100



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
            
        for l in lens:
            l = len(l)
            if l == 0:
                continue
            
            if l not in lens_count:
                lens_count[l] = 0
            lens_count[l] = lens_count[l] + 1
    return lens_count


