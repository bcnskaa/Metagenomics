from Bio import AlignIO
import glob
import sys
import os
import random


def print_usage():
    print("Usage:")
    print("  convert_aln2phy.py -x FILE")
    print("")
    print("  Option (-x):")
    print("    -c  convert the specified FASTA file/files with specified extension into Phylip format.")
    print("    -p  perform Phylip Protdist analysis")
    print("    -d  perform Phylip DNAdist analysis")
    print("    -h  print usage")    
    print("")
    


def convert_aln2phy(aln_fn_pattern):
    faa_fns = glob.glob("./" + aln_fn_pattern)
    if len(faa_fns) == 0:
        print("No file is found. (pattern used: " + aln_fn_pattern)
    print("Converting ")
    

    for faa_fn in faa_fns:
        print("Processing " + faa_fn)
        out_fn = faa_fn + ".phy"
    
#        alignments = AlignIO.read(open(faa_fn, "r"), "fasta")
        
#        OUT = open(out_fn, "w")
#        AlignIO.write(alignments, OUT, "phylip")
#        OUT.close()
        
        AlignIO.convert(faa_fn, "fasta", out_fn, "phylip")



"""
    
seqboot_options = $f "\nR\n" $rep_n "\nY"  $randint 

"""
def do_phylip_dist(files, type="protein", prog_path="/home/siukinng/tools/phylip/exe", resample=False, replicate_n=100):
    if type == "protein":
        prog = "protdist"
        options = "\nP\nP\nY\n"
    else:
        prog = "dnadist"
        options = "\nD\nD\nG\nY\n"
    phy_fns = glob.glob(files)
    
    if len(phy_fns) == 0:
        print("No file specified, abort now.")
        return None
    
    print("Performing " + type + " distance analysis (" + prog + ")")
    
    #for f in *.phy;do echo "Processing $f";echo -e "$f\nP\nP\nY\n" > current.cmd; ~/tools/phylip/exe/protdist < current.cmd > screenout; mv outfile $f.mtx;done
    for phy_fn in phy_fns:

        print("processing " + phy_fn)
        tmp_fn = "current." + get_tmp_fn() + ".cmd"
        
        if os.path.isfile("outfile"):
            os.remove("outfile")
        
        mtx_fn = phy_fn + ".mtx"
        cmd = "echo -e '" + phy_fn + options + "' > " + tmp_fn + "; " + prog_path + "/" + prog + " < " + tmp_fn + " > screenout; mv outfile " + mtx_fn
        os.system(cmd)
        
        tree_fn = phy_fn + ".kitsch.tree"
        cmd = "echo -e '" + mtx_fn + "\nY\n' > " + tmp_fn + "; "  + prog_path + "/kitsch < " + tmp_fn + " > screenout; mv outtree " + tree_fn
        os.system(cmd)  
        
        if os.path.isfile(tmp_fn):
            os.remove(tmp_fn)
            
            
            

def get_tmp_fn(fn_len=8):
    return ''.join(random.choice('0123456789ABCDEF') for i in range(fn_len))
        
    


def main(argv):
    if len(argv) != 2:
        print_usage()
        return None
    
    option = argv[0]
    files = argv[1]
    
    if option == "-c":
        convert_aln2phy(files)
    elif option == "-p":
        do_phylip_dist(files)
    elif option == "-d":
        do_phylip_dist(files, type="nucleotide")
    else:
        print_usage()
        return None
        
        


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