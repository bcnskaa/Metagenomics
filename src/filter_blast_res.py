# 2014 SKWoolf bcnskaa AT gmail DOT com
import os.path
import sys
import getopt
import csv
import operator
from collections import defaultdict


# Main 
def main(argv):
    try:
        opts, args = getopt.getopt(argv,"hqsi:o:b:l:p:k:")
    except getopt.GetoptError:
        print_usage()
        sys.exit(2)
        
    infn = None
    outfn = None
    t_aln_len = 0
    t_bitscore = 0.0
    t_identity = 0.0
    flag_unique_query = False
    flag_unique_subject = False
    nr_sid_outfn = None
    # used with -k option
    sid_desc_fn = None


    for opt, arg in opts:
        if opt == '-h':
            print_usage()
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            infn = arg
        elif opt in ("-o", "--outfile"):
            outfn = arg
        elif opt == "-q":
            flag_unique_query = True
            print "Export the best hit for every query id"
        elif opt == "-s":
            #nr_sid_outfn = arg
            flag_unique_subject = True
            #print "Non-redundant subject id will be dumped to ", nr_sid_outfn
            #print "Non-redundant subject id will be dumped"
        elif opt == "-b":
            t_bitscore = float(arg)
            print "Bit-score Threshold: ", t_bitscore
        elif opt == "-p":
            t_identity = float(arg)
            print "Sequence Identity Threshold:", t_identity
        elif opt == "-l":
            t_aln_len = int(arg) 
            print "Alignment Length Threshold:", t_aln_len
        elif opt == "-k":
            sid_desc_fn = arg
            print "Subject ID map:", sid_desc_fn
        else:
            print("Unrecognized option", opt)
            print_usage()
            sys.exit(0)

    # Make sure infilename and outfilename are specified
    #if (infn is None and outfn is None):
    if (infn is None):
        print_usage()
        sys.exit(0)
    
    # Make sure outfilename is not same to infilename!
    if(infn == outfn):
        print "infile and outfile are the same!"
        sys.exit(0)
    
    
    # outputfilename is not specified, we name it! 
    if(outfn is None):
        outfn = infn + "."
        if t_bitscore > 0.0:
            outfn += "b" + `t_bitscore`
        if t_aln_len > 0:  
            outfn += "l" + str(t_aln_len)
        if t_identity > 0:  
            outfn += "p" + `t_identity`
        if flag_unique_query:  
            outfn += "q"
    
    # Specify an outfilename for dumping unique subject id
    if flag_unique_subject:
        nr_sid_outfn = outfn + ".sid"
        print "Non-redundant subject id will be dumped to", nr_sid_outfn
    
    # finalize the name of outfile
    outfn += ".bla"
    
    # import the blast results
    res = import_blast_res(infn, t_aln_len, t_bitscore, t_identity)
    
    # Export best hit for every query id
    if flag_unique_query:
        res = report_besthit(res)
    
    # Read description for subject ids if provided
    if sid_desc_fn is not None:
        sid_desc = import_subject_id_desc(sid_desc_fn)


    nr_sids = []
    
    # Export the filtered results to an outfile
    print "Exporting filtered results to", outfn
    export_n = 0
    with open(outfn, "w") as out:
        for qid, hits in res.items():
            
            # Deal with unique query id
            if flag_unique_query:
                export_n = export_n + 1
                # Replace sid with its description defined in description list
                if sid_desc_fn is not None:
                    hits[1] = sid_desc[hits[1]]
                
                # Chomp first whitespace character (if exists)
                hits[1] = hits[1].lstrip()
                    
                # Collect the sid and dump them later
                if nr_sid_outfn is not None:
                    if hits[1] not in nr_sids:
                        nr_sids.append(hits[1])
                        #print hits[1]
                            
                #print "\t".join(hits)
                out.write("".join(["\t".join(hits), "\n"]))
            else:  # Deal with all query id
                export_n = export_n + len(hits)
                for hit in hits:
                    # Replace sid with its description defined in description list
                    if sid_desc_fn is not None:
                        hit[1] = sid_desc[hit[1]]
                
                    # Chomp first whitespace character (if exists)
                    hit[1] = hit[1].lstrip()
                    
                    # Collect the sid and dump them later
                    if nr_sid_outfn is not None:
                        if hit[1] not in nr_sids:
                            nr_sids.append(hit[1])
                            #print hit[1]
                        
                    #print "\t".join(hit)
                    out.write("".join(["\t".join(hit), "\n"]))
                #print "\n".join(["\t".join(hit) for hit in hits])
                #print "\n".join([r for r in res[qid]])
            #[out.write("\t".join([r, "\n"])) for r in res[qid]]
    print export_n, "exported"


    # Sort and dump non-redundant to stdout
    if nr_sid_outfn is not None:
        # Sort the sids
        nr_sids.sort();
        #print "\n".join(nr_sids);
        with open(nr_sid_outfn, "w") as nr_sid_out:
            nr_sid_out.write("\n".join(nr_sids))
        


# Report the best hit for every query id 
def report_besthit(blast_results):    
    blast_besthits = {}
    
    for qid in blast_results.keys():
        blast_besthits[qid] = sorted(blast_results[qid], key=lambda x: x[11], reverse=True)[0]
    
    return blast_besthits


# Import information for subject id
# File format for the input file: sid\tdescript1...
def import_subject_id_desc(infilename):
    # 
    sid_desc = defaultdict(list)
    
    print "Importing from", infilename
    
    # infilename should be a tab-delimited csv file
    with open(infilename, "r") as infile:
        
        # read line
        for line in infile:
            line = line.replace("\n", "")
            [sid, desc] = line.split("\t", 1)
            #print sid, "=", desc
            sid_desc[sid] = desc
        
        # Close the input stream
        infile.close()
        
    # Exit
    return sid_desc


# Print the usage of this script
def print_usage():
    print("A simple python script for filtering blast result (generated with -outfmt 6 option) using threshold values")
    print(" ")
    print("Usage:")
    print("  python filter_blast_res.py -i BLAST-RESULT-INFILE -o FILTER-OUTFILE [-q] [-b BITSCORE-CUTOFF] [-l ALIGNMENT-LENGTH-CUTOFF] [-p PERCENTAGE-IDENTITY]")
    print("      -i STRING  Input file generated by BLAST with -m6 option")
    print("      -o STRING  Output file [optional]")
    print("      -h         Report the best hit for every query id [optional, default=disable]")
    print("      -b FLOAT   Bit-score threshold for filtering [optional, default=0.0]")
    print("      -l INT     Alignment length threshold for filtering [optional, default=0]")
    print("      -p FLOAT   Identity percentage threshold for filtering [optional, default=0.0]")
    print("      -s STRING  Dump non-redundant subject id")
    print("      -q STRING  Dump non-redundant query id")
    print(" ")
    print(" ") 
    print("Ver 0.2")


# Import blast results
def import_blast_res(infilename, thres_aln_len, thres_bitscore, thres_identity):
    # 
    blast_res = defaultdict(list)
    
    line_n = 0
    
    print "Importing from", infilename
    
    # infilename should be a tab-delimited csv file
    with open(infilename, "r") as infile:
        input_data = csv.reader(infile, delimiter='\t')
        
        # parse rows       
        for row in input_data:
            flag_discard = False
            
            # Track the current position
            line_n = line_n + 1
            
            try:
                # cast into corresponding data types
                bit_score = float(row[11])
                percentage_identity = float(row[2])
                alignment_length = float(row[3])
            except ValueError:
                print row[11], "is not a value. Discarded."
                flag_discard = True

            
            # Check if values exceed threshold values
            if (thres_aln_len is not None and alignment_length < thres_aln_len):
                flag_discard = True
            if (thres_identity is not None and percentage_identity < thres_identity):
                flag_discard = True
            if (thres_bitscore is not None and bit_score < thres_bitscore):
                flag_discard = True
                
            if not flag_discard:
                blast_res[row[0]].append(row)
        
        # Close the input stream
        infile.close()
    
    # Exit
    return blast_res  

 
# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])