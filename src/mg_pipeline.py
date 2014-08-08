#!/share/apps/Python-2.7.4/bin/python
import os
import sys
import getopt
import csv
import operator
from collections import defaultdict
import inspect
import subprocess
from subprocess import call
import glob
import os.path



# Global variables
HOME = "/home/siukinng"
TOOLS_HOME = HOME + "/tools"
IDBA_UD_HOME = TOOLS_HOME + "/idba_ud"
FASTQC_HOME = TOOLS_HOME + "/FastQC"
SEQTK_HOME = TOOLS_HOME + "/seqtk"
FASTQ_MCF_HOME = TOOLS_HOME + "/fastq-mcf"
MAXBIN_HOME = TOOLS_HOME + "/MaxBin"
BOWTIE2_HOME = TOOLS_HOME + "/"

FASTQ_EXT = "fq"
FASTA_EXT = "fa"

CUR_WD = "."

# Sequences of sequencing adaptors
ADAPTOR_SEQS = HOME + "/db/Sequencing/adaptors.fa"




# Main 
def main(argv):
    # Current working directory
    #CUR_WD = argv[0]
    #read_fns = glob.glob(CUR_WD + "/*." + FASTQ_EXT)
    #[preprocess(fn) for fn in read_fns]
    
    read_1_fn = argv[0]
    read_2_fn = argv[1]
    read_1_processed_fn = preprocess(read_1_fn)
    read_2_processed_fn = preprocess(read_2_fn)

    read_fns = run_fastq_mcf(read_1_processed_fn, read_2_processed_fn, ADAPTOR_SEQS)
    run_idba_ud(read_fns[0], read_fns[1])
    

######### Preprocessing 

# Preprocess FASTQ
def preprocess(read_fn): 
    print_status("Preprocessing " + read_fn)
    
    #run_FastQC(read_fn)
    
    processed_outfn = read_fn.replace("."+FASTQ_EXT, ".trimmed."+FASTQ_EXT)
    #run_seqtk(read_fn, processed_outfn)
    
    # Generate FastQC reports after trimming
    #run_FastQC(processed_outfn)    
    
    # Return my new name
    return processed_outfn
    
    
def run_FastQC(read_fn, report_outdir="fastqc_report"):
    print_status("Processing " + read_fn)    
    
    report_outdir_path = CUR_WD + "/" +  report_outdir
    if not os.path.exists(report_outdir_path):
        os.makedirs(report_outdir_path)
    
    cmd = FASTQC_HOME + "/fastqc " + read_fn + " --outdir=" + report_outdir_path
    
    print_status(cmd)
    
    #subprocess.Popen(cmd, shell=True)
    os.system(cmd)


def run_seqtk(read_fn, trimmed_read_fn):
    print_status("Processing " + read_fn)    
    
    cmd = SEQTK_HOME + "/seqtk trimfq " + read_fn + " > " + trimmed_read_fn
    print_status(cmd)
    os.system(cmd)


# http://onetipperday.blogspot.hk/2012/08/three-ways-to-trim-adaptorprimer.html
def run_fastq_mcf(read_1_fn, read_2_fn, adaptor_seq_fn, min_read_len=16, min_trim_quality=15, trim_win_len=4, N_percent=10):
    print_status("Processing" + read_1_fn + "and" + read_2_fn)    
    #"$FASTQ_MCF_HOME/fastq-mcf -o $outputfile_1 -o $outputfile_2 -l 16 -q 15 -w 4 -x 10 $ADAPTOR_SEQS $inputfile_1 $inputfile_2"
    
    read_1_outfn = read_1_fn.replace("."+FASTQ_EXT, ".cleaned." + FASTQ_EXT)
    read_2_outfn = read_2_fn.replace("."+FASTQ_EXT, ".cleaned." + FASTQ_EXT)
    
    #cmd = FASTQ_MCF_HOME + "/fastq-mcf -o " + read_1_outfn + " -o " + read_2_outfn + " -l " + min_read_len + " -q " + min_trim_quality + " -w " + trim_win_len + " -x " + N_percent + " " + adaptor_seq_fn + " " + read_1_fn + " " + read_2_fn
    cmd = FASTQ_MCF_HOME + "/fastq-mcf -o " + str(read_1_outfn) + " -o " + str(read_2_outfn) + " -l " + str(min_read_len) + " -q " + str(min_trim_quality) + " -w " + str(trim_win_len) + " -x " + str(N_percent) + " " + adaptor_seq_fn + " " + read_1_fn + " " + read_2_fn
    
    print_status(cmd)
    os.system(cmd)
    
    return [read_1_outfn, read_2_outfn]
    

####### Assemble stage #########
# IDBA-UD

def run_idba_ud(read_1_fn, read_2_fn, outdir="idba_ud", min_contig=1200, mink=20, maxk=80, step=10, num_threads=16):
    print_status("Initializing IDBA-UD")
    
    print_status("Merging mate files")
    
    outfn = read_1_fn.replace("." + FASTQ_EXT, "") + read_2_fn.replace("." + FASTQ_EXT, "") + "." + FASTQ_EXT
    cmd = IDBA_UD_HOME + "/bin/fq2fa --merge --filter " + read_1_fn + " " + read_2_fn + " " + outfn
    print_status(cmd)
    os.system(cmd)
    
    if not os.path.isfile(outfn):
        print_status("Problem on merging mate files")
        return False
    
    print_status("Provoking IDBA_UD")
    cmd = IDBA_UD_HOME + "/bin/idba_ud --long_read " + outfn + " -o " + outdir + " --mink " + str(mink) + " --maxk " + str(maxk) + " --num_threads " + str(num_threads) + " --min_contig " + str(min_contig) + " --pre_correction --step " + str(step)
    print_status(cmd)
    os.system(cmd)

    return True
    
    
    
####### Binning #############
def run_MaxBin(contig_fn, merged_read_fn, thread=16):
    print_status("Initializing MaxBin")    
    outfn = contig_fn.replace("."+FASTA_EXT, ".maxbin")
    
    cmd = MAXBIN_HOME + "/run_MaxBin.pl -contig " + contig_fn + " -out " + outfn + " -thread " + str(thread) + " -plotmarker -reads " + merged_read_fn
    print_status(cmd)
    os.system(cmd)
    
    return outfn


def run_ESOM():
    print_status("Initializing ESOM")    
    


####### Functional annotation stage #########
# IDBA-UD
def blastn():
    print_status("Initializing for blastn")

def blastp():
    print_status("Initializing for blastp")

    
def blastx():
    print_status("Initializing blastx()")    


def filter_blast():
    print_status("filtering blast results")    
  

def prepare_reference_genome():
    print_status("Preparing reference genomes")        


####### CAZy Specific stage #########






####### Auxillary routines #########


# 
def parse_opt():
    print_status("calling parsing opts")


# Print the usage of this script 
def print_usage():
    print("An integrated pipeline for processing metagenomic sequencing reads.")
    print(" ")
    print("Usage:")
    print("  python -i BLAST-RESULT-INFILE -o FILTER-OUTFILE [-q] [-b BITSCORE-CUTOFF] [-l ALIGNMENT-LENGTH-CUTOFF] [-p PERCENTAGE-IDENTITY]")
    print("      -i STRING  Input file generated by BLAST with -m6 option")
    print(" ")
    print(" ") 
    print("Ver 0.2")



# 
def sumarize_read(read_fn):
    print_status("Summarising Read file")


# Assert 
def assert_process():
    print_status("assert_process()")


# 
def import_settings():
    print_status("importing settings")


def log(log_msg):
    print_status(log_msg)


# Print status
def print_status(msg):
    caller_name = inspect.stack()[1][3]
    msg = "[ " + caller_name + " ]  " + msg
    print msg
    

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])