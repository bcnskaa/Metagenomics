from Bio import AlignIO
import glob

faa_fns = glob.glob("./*.aln")
for faa_fn in faa_fns:
    print("Processing " + faa_fn)
    out_fn = faa_fn + ".phy"

    alignments = AlignIO.read(open(faa_fn, "r"), "fasta")
    #nr_alignments = {a.id:a for a in alignments}
    #nr_alignments = [nr_alignments[k] for k in nr_alignments]
    
    #OUT = open(out_fn, "w")
    #for a in alignments:
    #    AlignIO.write(a, OUT, "phylip-relaxed")
        
    #OUT.close()
    
    OUT = open(out_fn, "w")
    AlignIO.write(nr_alignments, OUT, "phylip-relaxed")
    #AlignIO.convert(faa_fn, "fasta", out_fn, "phylip")
    #AlignIO.convert(faa_fn, "fasta", out_fn, "phylip-relaxed")
    OUT.close()




"""
records = SeqIO.parse(open("combined.faa.aln"), "fasta")

ids = []
seqs = []
for r in records:
    if r.id not in ids:
        seqs.apend(r)
        ids.append(r.id)
        




"""