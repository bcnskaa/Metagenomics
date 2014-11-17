
import re
import glob
from Bio import SeqIO
import string
import random
import sys
import os

def main(argv):
    wd_dir = "."
    aln_fn_ext = "aln"
    aln_fns = glob.glob(wd_dir + "/*." + aln_fn_ext)
    aln_format = "fasta"
    export_fn_ext = "renamed"
    map_fn = "map.lst"
     
    id_map = {}
    for aln_fn in aln_fns:
        print("Processing " + aln_fn)
        seqs = SeqIO.index(aln_fn, aln_format)
        out_fn = aln_fn + "." + export_fn_ext
        with open(out_fn, "w") as OUT:
            for id in seqs.keys():
                cflag_existed = False
                while not cflag_existed:
                    new_id = generate_id(id)
                    if new_id not in id_map.keys():
                        cflag_existed = True
                #print("Maping " + id + " to " + new_id)
                
                id_map[new_id] = id + "\t" + os.path.basename(aln_fn)
                seq = seqs[id]
                seq.id = new_id
                seq.name = new_id
                seq.description = ""
                
                SeqIO.write(seq, OUT, "fasta")
                
    with open(map_fn, "w") as OUT:
        for id in id_map:
            OUT.write(id + "\t" + id_map[id] + "\n")



def generate_id(id, n=10):
    digit_letters = (string.digits + string.letters)
    
    items = re.split(r"_|-|\.", id)
    id = ""
    if len(items) == 6:
        id = id + items[0][0] + items[1][0] + items[2][1] + "_" + str(int(items[3])) + "_"
        rlen = n - len(id)
        id = id + "".join([random.choice(digit_letters) for i in range(0, rlen)])
    else:
        id = "".join([random.choice(digit_letters) for i in range(0, n)])
    return id


# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    