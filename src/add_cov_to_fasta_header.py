import os
from Bio import SeqIO
import math
import sys


def main(argv):
    if len(argv) != 3:
        print("Missing argument")
        print("Usage: " + sys.argv[0] + " FASTA-INFILE ABUND-INFILE NEW-FASTA-OUTFILE")
        return
    
    fasta_fn = argv[0]
    abund_fn = argv[1]
    fasta_ofn = argv[2]
    
    with open(abund_fn) as IN:
        lines = IN.read().splitlines()
    abund_info = {l.split("\t")[0]:int(round(float(l.split("\t")[1]))) for l in lines}    
    
    export_count = 0
    skip_count = 0
    total_count = 0
    records = SeqIO.parse(fasta_fn, "fasta")
    with open(fasta_ofn, "w") as OUT:
        for r in records:
            total_count = total_count + 1
            cov = 0
            if r.id in abund_info.keys():
                cov = abund_info[r.id]
                r.id =  r.id + "_[cov=" + str(cov) + "]"
                r.description = ""
                r.name = ""
                SeqIO.write(r, OUT, "fasta")
                export_count = export_count + 1
            else:
                #print(r.id + " skipped (length=" + str(len(r.seq)) + ")")
                skip_count = skip_count + 1
    print("Number of sequence processed=" + str(total_count) + ", "+ str(export_count) + " exported, " + str(skip_count) + " skipped")
                          


# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])