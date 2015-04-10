# 2014 SKWoolf bcnskaa AT gmail DOT com
import os.path
import sys
import getopt
from collections import defaultdict
from Bio.SeqUtils.CheckSum import seguid


# Require biopython
# To run this script, path to biopython libraries has to be included in PYTHONPATH
from Bio import SeqIO



# Main 
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hi:o:n:")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)

    fasta_char_per_line = 70
    header = None
    fasta_ifn = None
    fasta_ofn = None

    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i"):
            fasta_ifn = arg
        elif opt in ("-o"):
            fasta_ofn = arg
        elif opt == "-n":
            header = arg
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)


    # Make sure everything is provided
    if(fasta_ifn is None or fasta_ofn is None or header is None):
        print("Missing argument")
        print_usage()
        sys.exit(0)
        
       
    if(fasta_ofn is None):
        fasta_ofn = fasta_ifn + ".merged.fa"
    
    [exported_n, seq] = merge(fasta_ifn, fasta_ofn, header)
    if exported_n > 0:
        export_fasta(seq[0], seq[1], fasta_ofn)
        print("Number of sequence merged: " + str(exported_n) + " (" + str(len(merged_seq)) + "bp)")
    else:
        print("Nothing to merge.")
         
         

def dedup(fa_fn, fa_ofn=None,dedup_by_id=True):
    if fa_ofn is None:
        fa_ofn = fa_fn + ".dedup.fa"
    
    nr_records = get_dedup_seqs(fa_fn, dedup_by_id)
    print("Exporting dedup sequences to " + fa_ofn)
    OUT = open(fa_ofn, "w")
    SeqIO.write(nr_records, OUT, "fasta")
    OUT.close()



def get_dedup_seqs(fa_fn, dedup_by_id=True):        
    records = SeqIO.parse(fa_fn, "fasta")
    
    if dedup_by_id:
        print("Checking sequence ids for duplication")
    else:
        print("Checking sequence for duplication")
    
    
    dup_n = 0
    seq_guids = set()
    for record in records:
        if dedup_by_id:
            seq_guid = seguid(record.id)
        else:
            seq_guid = seguid(record.seq)
        
        if seq_guid in seq_guids:
            dup_n = dup_n + 1
            continue
        
        seq_guids.add(seq_guid)
        yield record
    
    print(str(dup_n) + " sequences are duplicates.")



def merge(fasta_ifn, header):
#def merge(fasta_ifn, fasta_ofn, header, fasta_char_per_line = 70):
    # Call the Bio.SeqIO
    # Create a dict index on the fasta sequence file
    seqs = SeqIO.index(fasta_ifn, "fasta")
    merged_seq = ""
    exported_n = 0
    for seq_id in seqs.keys():
        merged_seq = merged_seq + str(seqs[seq_id].seq) 
        exported_n += 1
        
#     seq_lines = [merged_seq[i:i+fasta_char_per_line] for i in range(0, len(merged_seq), fasta_char_per_line)]
#     with open(fasta_ofn, "w") as OUT:
#         OUT.write(">" + header + "\n")
#         OUT.write("\n".join(seq_lines))
#         OUT.write("\n")

    return [exported_n, [header, merged_seq]]
    


def export_fastas(seqs, fa_ofn=None, OUT=None, fasta_char_per_line = 70):
    for seq_id in seqs.keys():
        export_fasta(seq_id, seqs[seq_id], fa_ofn=fa_ofn, OUT=OUT, fasta_char_per_line = fasta_char_per_line)


def export_fasta(header, seq, fa_ofn=None, OUT=None, fasta_char_per_line = 70):
    #print("Exporting " + header + " (" + str(len(seq)) + " bp)" )
    
    close_on_exit = False
    if OUT is None:
        if fa_ofn is not None:
            OUT = open(fa_ofn, "w")
        else:
            fa_ofn = header + ".fa"
            OUT = open(fa_ofn, "w")
        close_on_exit = True

    seq_lines = [seq[i:i+fasta_char_per_line] for i in range(0, len(seq), fasta_char_per_line)]
        
    #    print("Exporting " + header + " to " + fa_ofn)
    #else:
    #    print("Appending " + header + " to OUT handler")
        
    #with open(fasta_ofn, "w") as OUT:
    OUT.write(">" + header + "\n")
    OUT.write("\n".join(seq_lines))
    OUT.write("\n")
    
    if close_on_exit:
        OUT.close()
    
    return fa_ofn




# Print the usage of this script 
def print_usage():
    print("  python merge_fa.py -i FASTA-INFILE -o FASTA-OUTFILE -n HEADER")
    print("      -i STRING  Fasta infile")
    print("      -o STRING  Merged outfile")
    print("      -n STRING  Header of merged fasta sequence")
    print(" ")



# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    