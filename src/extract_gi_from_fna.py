import os.path
import sys


def main(argv):
    if len(argv) != 2:
        print_usage()
        sys.exit(1);
        
    infile = argv[0]
    outfile = argv[1]
    
    headers = extract_fasta_header(infile)

    print "Number of header=", len(headers)

    gi_headers = parse_headers(headers)
    #print [":".join([k, gi_headers[k], "\n"]) for k in gi_headers.keys()]

    # export the gi and name to an outfile
    with open(outfile, "w") as outfile:
        for gi in gi_headers.keys():
            outfile.write("\t".join([gi, gi_headers[gi], "\n"]))
        #outfile.write(["\t".join([gi, gi_headers[gi]]) for gi in gi_headers.keys()])
        outfile.close()
    

def print_usage():
    print("NA")

# FASTA header
def extract_fasta_header(fasta_fn):
    fasta_headers = []
    with open(fasta_fn, "r") as infile:
        fasta_headers = [line.replace("\n", "") for line in infile if line.startswith(">")]

    return fasta_headers


def parse_headers(fasta_headers):
    gi_headers = {}
    for header in fasta_headers:
        gi = header.split("|")[1]
        name = header.split("|")[4]
        gi_headers[gi] = name
        
    return gi_headers

# Invoke the main function
if __name__ == "__main__":
    
    main(sys.argv[1:])