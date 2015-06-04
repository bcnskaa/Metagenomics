from __future__ import print_function
from __future__ import division
import csv
from itertools import izip
import sys
import numpy

def print_usage():
    #print("Usage: python transpose_csv.py option INFILE OUTFILE [CUTOFF] [ROW_ID_TO_BE_DROPPED] [OPTION-SPECIFIC-ARGS]")
    print("Usage: python transpose_csv.py option INFILE OUTFILE [CUTOFF]")
    print("    Option:")
    print("            transpose   transpose the csv.")
    print("            filter      filter the row by a cutoff.")  
    print("            drop        drop row(s) by ids: ID;ID;ID;") 
    
    
def main(argv):
    if len(argv) < 2:
        print_usage()
        return 0
    
    option = argv[0]
    
    if option == "transpose":
        transpose(argv[1:])
    elif option == "filter":
        filter(argv[1:])
    elif option == "drop":
        drop_row(argv[1:])
    else:
        print_usage()
        return 0 
        

"""

"""
def filter(argv):
    if len(argv) < 3:
        print_usage()
        return 0
    infn = argv[0]
    outfn = argv[1]
    cutoff = float(argv[2])
    
    print("Filtering " + infn + " with cutoff=" + str(cutoff))
    
    # -c
    filtering_col = False
    # -n
    normalized = False

    # Parse options
    if len(argv) > 3:
        options = argv[3:]
        for option in options:
            if option == "-c":
                filtering_col = True
            elif option == "-n":
                normalized = True
            else:
                print("Unknown option: " + option)
    
    # Filter
    IN = open(infn)
    data = IN.read().splitlines()
    data = [d.split("\t") for d in data]
    header = data[0]
    del data[0]
    
    filtered_sums = [0.0 for d in header[1:]]
    filtered_total = 0
    
    # Row count
    filtered_n = 0
    processed_n = 0
    exported_n = 0
    OUT = open(outfn, "w")
    if not filtering_col:
        
        OUT.write("\t".join(header) + "\n")

        if normalized:
            # Calculate column sums
            mtx = [map(float, d[1:]) for d in data]
            #mtx = map(float, mtx) 
            # Calculate the col_sum
            col_sums = numpy.sum(mtx, axis=0)
        
        print("\t".join(map(str,col_sums)))

        # Normalized by row sum
        for d in data:
            processed_n += 1
            
            row_id = d[0]
            
            row_dat = map(float, d[1:])
            if normalized:
                row_dat = [rd / cs for rd, cs in zip(row_dat, col_sums)] 
            row_sum = sum(row_dat) 

            
            if row_sum < cutoff:
                filtered_n += 1
                filtered_total += sum(row_dat)
                #print(str(sum(row_dat)))
                filtered_sums = [a + b for (a, b) in zip(filtered_sums, row_dat)]
                continue
                
#             if normalized:
#                 row_dat_norm = [r / row_sum for r in row_dat]
#                 OUT.write(row_id + "\t" + "\t".join(map(str, row_dat_norm)) + "\n")
#                 exported_n += 1
#             else:

            OUT.write(row_id + "\t" + "\t".join(map(str, row_dat)) + "\n")
            exported_n += 1
        
        # Exported the filtered_sums
        OUT.write("Filtered" + "\t" + "\t".join(map(str, filtered_sums)) + "\n")
            
    print("Number of row processed: " + str(processed_n))
    print("Number of row filtered: " + str(filtered_n))
    print("Number of row exported: " + str(exported_n))
    print("Filtered total: " + str(filtered_total))
    
#     else:
#         mtx = [d[1:] for d in data]
#         mtx = map(float, mtx)
#         
#         # Calculate the col_sum
#         col_sums = numpy.sum(mtx, axis=1)    
#         
#         excluded_list = [i for i, cs in enumerate(col_sums) if cs < cutoff]
#         
#

    IN.close()
    OUT.close()
    

"""
Usage: 
 python transpose_csv.py drop INFILE OUTFILE "GZ-Seed_Y0+nr.m8;SWH-Seed_Y0+nr.m8"
 python transpose_csv.py drop samples-tax_id+go.samples.g.clustered.seed_drop.transposed samples-tax_id+go.samples.g.clustered.seed_drop.transposed "GZ-Seed_Y0+nr.m8;SWH-Seed_Y0+nr.m8"
 
 
"""
def drop_row(argv):
    if len(argv) < 3:
        print_usage()
        return 0
    infn = argv[0]
    outfn = argv[1]
    row_ids = argv[2].split(";")
    
    fixed = True
     
    if len(argv) > 3:
        options = argv[3:]
        for option in options:
            if option == "-w":
                fixed = False
            else:
                print("Unknown option: " + option)
    
    if fixed:
        print("The following row IDs will be dropped from " + infn)
        print("\n".join(row_ids))
    else:
        print("Row IDs with the following keywords will be dropped from " + infn)
        print("\n".join(row_ids))  

    IN = open(infn)
    OUT = open(outfn, "w")
    
    processed_n = 0
    discarded_n = 0
    for line in IN:
        processed_n += 1
        # Export comment lines
        if line.startswith("#"):
            OUT.write(line)
        else:
            # Check if the row_id specified in the list
            if fixed:
                if line.split("\t")[0] not in row_ids:
                    OUT.write(line)
                else:
                    discarded_n += 1
            else:
                matched = False
                id = line.split("\t")[0]
                for row_id in row_ids:
                    if row_id in id:
                        matched = True
                        break
                
                if not matched:
                    OUT.write(line)
                else:
                    discarded_n += 1

    print("Processed row: " + str(processed_n))
    print("Discarded row: " + str(discarded_n))
    
    IN.close()
    OUT.close()
    
     


def transpose(argv): 
    if len(argv) != 2:
        print_usage()
        return 0
    infn = argv[0]
    outfn = argv[1]
    
    print("Transposing " + infn)
    
    a = izip(*csv.reader(open(infn, "r"), delimiter='\t'))
    csv.writer(open(outfn, "w"), delimiter='\t').writerows(a)



# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])
    
    
"""

fn = "samples-tax_id+go.samples.transposed"
tax <- read.table(fn, header=T, comment.char="@", sep="\t", stringsAsFactors=F, row.names=1)
tax <- read.delim2(fn, header=T, comment.char="@", sep="\t", stringsAsFactors=F, row.names=1)

"""