from Bio import AlignIO
import glob
import sys
import os

def main(argv):
    faa_fns = glob.glob("./*.aln")
    for faa_fn in faa_fns:
        print("Processing " + faa_fn)
        out_fn = faa_fn + ".phy"
    
#        alignments = AlignIO.read(open(faa_fn, "r"), "fasta")
        
#        OUT = open(out_fn, "w")
#        AlignIO.write(alignments, OUT, "phylip")
#        OUT.close()
        
        AlignIO.convert(faa_fn, "fasta", out_fn, "phylip")


# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])


"""
for f in *.phy;do echo "Processing $f";echo -e "$f\nP\nP\nY\n" > current.cmd; ~/tools/phylip/exe/protdist < current.cmd > screenout; mv outfile $f.mtx;done

ls *.faa.renamed_faa.aln.phy | sed "s/.faa.renamed_faa.aln.phy//" | ~/tools/bin/parallel -j12 'mkdir {}; cd {}; echo "Processing {}"; echo -e "../{}.faa.renamed_faa.aln.phy\nP\nP\nY\n" > current.cmd; ~/tools/phylip/exe/protdist < current.cmd > screenout; mv outfile ../{}.faa.renamed_faa.aln.phy.mtx; cd ..; rm -r {}'





"""

"""

ggplot(contig_info2, aes(x=V2.y, y=V3.y, colour=V3.x)) + geom_point(aes(size=(V2.x), cex=2)) + geom_text(aes(label=V8), cex=3, hjust=-0.1,vjust=0) + xlim(10, 75) + ylim(50000,240000) + theme(panel.grid = element_blank(), panel.background = element_blank())

"""