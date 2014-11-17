from Bio import AlignIO
import glob
import sys
import os

def main(argv):
    faa_fns = glob.glob("./*.aln.renamed")
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

"""