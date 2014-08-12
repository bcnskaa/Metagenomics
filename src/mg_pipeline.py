#!/share/apps/Python-2.7.4/bin/python

"""


"""

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
import time
from time import strftime
from datetime import datetime
import math


# Global variables
HOME = "/home/siukinng"
TOOLS_HOME = HOME + "/tools"
SCRIPTS_HOME = TOOLS_HOME + "/scripts"
IDBA_UD_HOME = TOOLS_HOME + "/idba_ud"
FASTQC_HOME = TOOLS_HOME + "/FastQC"
SEQTK_HOME = TOOLS_HOME + "/seqtk"
FASTQ_MCF_HOME = TOOLS_HOME + "/fastq-mcf"
MAXBIN_HOME = TOOLS_HOME + "/MaxBin"
BOWTIE2_HOME = TOOLS_HOME + "/"
PRODIGAL_HOME = TOOLS_HOME + "/Prodigal"
BLAST_HOME = TOOLS_HOME + "/blast"
EMIRGE_HOME = TOOLS_HOME + "/EMIRGE"
BBMAP_HOME = TOOLS_HOME + "/BBMap"


FASTQ_EXT = "fq"
FASTA_EXT = "fa"
FASTA_PROT_EXT = "faa"
FASTA_DNA_EXT = "fna"

CUR_WD = "."

# Sequences of sequencing adaptors
DB_HOME = HOME + "/db"
ADAPTOR_SEQS = DB_HOME + "/Sequencing/adaptors.fa"

# NCBI 16S database
NCBI16S_DB = DB_HOME + "/Markers/ncbi16s/16SMicrobial.idx.fa"
NCBI16S_GI_IDX = DB_HOME + "Markers/ncbi16s/16SMicrobial.idx"

EMIRGE_16S_DB = EMIRGE_HOME + "/db/SSURef_111_candidate_db.fasta"

CAZY_DB = DB_HOME + "/Markers/CAZy/CAZy_id.lst.retrieved.faa"


VERBOSE_ONLY = True

LOG = True
LOG_FN= "pipeline.log"
LOG_FILE = open(LOG_FN, "w")

# Output directories
MARKER_OUTDIR = "Markers"
RESULT_OUTDIR = "Reports"


# Fastq Statistics


# Main 
def main(argv):
    # Current working directory
    #CUR_WD = argv[0]
    #read_fns = glob.glob(CUR_WD + "/*." + FASTQ_EXT)
    #[preprocess(fn) for fn in read_fns]
    
    #if LOG:
    #    LOG_FILE = open(LOG_FN, "w")
    
    CUR_WD = os.getcwd()
    print_status("Current Working Directory: " + CUR_WD)
    
    read_1_fn = argv[0]
    read_2_fn = argv[1]
    read_1_processed_fn = preprocess(read_1_fn)
    read_2_processed_fn = preprocess(read_2_fn)

    read_fns = run_fastq_mcf(read_1_processed_fn, read_2_processed_fn, ADAPTOR_SEQS)
    
    if len(read_fns) != 2:
        print_status("Problem on calling fastq_mcf results.")
        raise OSError, "Problem on processing fastq_mcf outputs, abort now."
    
    merged_read_fn = merge_paired_end_seq(read_fns[0], read_fns[1])
    
    if merged_read_fn is None:
        print_status("Problem on reading from merged_read_fn.")
        raise OSError, "Problem at merging reads files, abort now."
    
    if not run_idba_ud(merged_read_fn):
        print_status("Problem at completing IDBA_UD stage")
        raise OSError, "Problem at IDBA_UD stage, abort now."
    
    contig_fn = postprocess_idba_ud()  
    
    if contig_fn is None:
        print_status("Unable to locate contig.fa")
        raise OSError, "Contig file is not avaiable from IDBA_UD, abort now."
    

    maxbin_outdir = run_MaxBin(contig_fn, merged_read_fn)
    #run_idba_ud(read_fns[0], read_fns[1])
    
    prodigal_info = run_Prodigal(contig_fn)
    # {"source_fn":fna_infn, "outdir":prodigal_outdir, "prodigal_outfn":outfn, "protein_outfn":faa_outfn, "nucleotide_outfn":fna_outfn}


    # blastp -query ../../Prodigal/contig.fa.prodigal.faa -db ~/db/Markers/CAZy/CAZy_id.lst.retrieved.faa -outfmt 6 -out contig.fa.prodigal-CAZy_id.lst.retrieved.bla  -num_threads 16
    if not os.path.isfile(prodigal_info["protein_outfn"]):
        blastp(prodigal_info["protein_outfn"], CAZY_DB, outdir=MARKER_OUTDIR + "/CAZy")
    else:
        print_status("Warning: Predictive protein sequence is not available from Prodigal, step of mapping to CAZy database is skipped.")
        
    # Blast any 16s sequence fragments existed in newly assembled contig sequences
    blastn(contig_fn, NCBI16S_DB, outdir=MARKER_OUTDIR + "/ncbi16s")
    
    
    running_EMIRGE(read_fns[0], read_fns[1])
    
    


    # Close the log stream
    if LOG:
        LOG_FILE.close()
    

######### Preprocessing 

# Preprocess FASTQ
def preprocess(read_fn): 
    print_status("Preprocessing " + read_fn)
    
    run_FastQC(read_fn)
    
    processed_outfn = read_fn.replace("."+FASTQ_EXT, ".trimmed."+FASTQ_EXT)
    
    run_seqtk(read_fn, processed_outfn)
    
    # Generate FastQC reports after trimming
    run_FastQC(processed_outfn)    
    
    # Return my new name
    return processed_outfn
    
    
def run_FastQC(read_fn, report_outdir="fastqc_report"):
    print_status("Processing " + read_fn)    
    
    report_outdir_path = CUR_WD + "/" +  report_outdir
    if not os.path.exists(report_outdir_path):
        os.makedirs(report_outdir_path)
    
    cmd = FASTQC_HOME + "/fastqc " + read_fn + " --outdir=" + report_outdir_path
    
    print_status("command: " + cmd)
    
    #subprocess.Popen(cmd, shell=True)
    if not VERBOSE_ONLY:
        os.system(cmd)


def run_seqtk(read_fn, trimmed_read_fn):
    print_status("Processing " + read_fn)    
    
    cmd = SEQTK_HOME + "/seqtk trimfq " + read_fn + " > " + trimmed_read_fn
    print_status("command: " + cmd)
    
    if not VERBOSE_ONLY:
        os.system(cmd)


# http://onetipperday.blogspot.hk/2012/08/three-ways-to-trim-adaptorprimer.html
def run_fastq_mcf(read_1_fn, read_2_fn, adaptor_seq_fn, min_read_len=16, min_trim_quality=15, trim_win_len=4, N_percent=10):
    print_status("Processing" + read_1_fn + "and" + read_2_fn)    
    #"$FASTQ_MCF_HOME/fastq-mcf -o $outputfile_1 -o $outputfile_2 -l 16 -q 15 -w 4 -x 10 $ADAPTOR_SEQS $inputfile_1 $inputfile_2"
    
    read_1_outfn = read_1_fn.replace("."+FASTQ_EXT, ".cleaned." + FASTQ_EXT)
    read_2_outfn = read_2_fn.replace("."+FASTQ_EXT, ".cleaned." + FASTQ_EXT)
    
    #cmd = FASTQ_MCF_HOME + "/fastq-mcf -o " + read_1_outfn + " -o " + read_2_outfn + " -l " + min_read_len + " -q " + min_trim_quality + " -w " + trim_win_len + " -x " + N_percent + " " + adaptor_seq_fn + " " + read_1_fn + " " + read_2_fn
    cmd = FASTQ_MCF_HOME + "/fastq-mcf -o " + str(read_1_outfn) + " -o " + str(read_2_outfn) + " -l " + str(min_read_len) + " -q " + str(min_trim_quality) + " -w " + str(trim_win_len) + " -x " + str(N_percent) + " " + adaptor_seq_fn + " " + read_1_fn + " " + read_2_fn
    
    print_status("command: " + cmd)
   
    if not VERBOSE_ONLY:
        os.system(cmd)
    
    return [read_1_outfn, read_2_outfn]
    
    

####### Statistics stage #########
"""
 
"""
#def running_bbmap(read_1_fn, read_2_fn, ihist_fn="bbamp.ihist", outdir=RESULT_OUTDIR, read_n=8000000):
def summarize_statistics(read_1_fn, read_2_fn, ihist_fn="bbamp.ihist", outdir=RESULT_OUTDIR, read_n=8000000):
    print_status("Estimating summary statistics for " + read_1_fn + " and " + read_2_fn)
    
    # "~/tools/BBMap/bbmerge.sh in1=GZ-cell_S1_R1.trimmed.cleaned.fq in2=GZ-cell_S1_R2.trimmed.cleaned.fq ihist=ihist.txt reads=1000000 &> test.out"
    
    if not os.path.exists(RESULT_OUTDIR):
        print_status("Output directory, " + RESULT_OUTDIR + ", does not exist, we will create it.")
        os.makedirs(RESULT_OUTDIR)
    
    print_status("Estimating maximum read length from " + read_1_fn)
    cmd = "wc -L " + read_1_fn + " | cut -d ' ' -f1"
    max_read_len = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    if len(max_read_len) > 0:
        max_read_len = int(max_read_len)
    else:
        raise OSError, "\nError: unable to obtain the maximum read length from " + read_1_fn
        
    print_status("Estimating maximum read length from " + read_2_fn)
    cmd = "wc -L " + read_2_fn + " | cut -d ' ' -f1"
    max_read_len2 = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    if len(max_read_len2) > 0:
        max_read_len2 = int(max_read_len2)
        if max_read_len2 > max_read_len:
            max_read_len = max_read_len2
    else:
        raise OSError, "\nError: unable to obtain the maximum read length from " + read_2_fn
        
    print_status("Estimating insert size")
    summary_fn = outdir + "/" + read_1_fn.replace("."+FASTQ_EXT, "") + "-" + read_2_fn.replace("."+FASTQ_EXT, "")
    
    cmd = BBMAP_HOME + "/bbmerge.sh in1=" + read_1_fn + " in2=" + read_2_fn + " ihist=" + outdir + "/" + ihist_fn + " reads=" + str(read_n) + " &> " + summary_fn
    print_status("command: " + cmd)
    
    #if not VERBOSE_ONLY:    
    #    os.system(cmd)
    os.system(cmd)

    if not os.path.isfile(summary_fn):
        print_status("Problem on estimating files")        
        raise OSError, "\nError: problem on processing read data. (Please check BBMap settings)"
    
    insert_size_mean = float(subprocess.Popen("cat " + summary_fn + "  | grep -e \"Avg Insert\" | cut -f3", shell=True, stdout=subprocess.PIPE).stdout.read())
    insert_size_sd = float(subprocess.Popen("cat " + summary_fn + "  | grep -e \"Standard\" | cut -f3", shell=True, stdout=subprocess.PIPE).stdout.read())
    read_1_count = int(subprocess.Popen("cat " + read_1_fn + "  | grep -c \"@\"", shell=True, stdout=subprocess.PIPE).stdout.read()) 
    read_2_count = int(subprocess.Popen("cat " + read_2_fn + "  | grep -c \"@\"", shell=True, stdout=subprocess.PIPE).stdout.read())

    if read_1_count != read_2_count:
        print_status("Warning: reads count between paired-end reads is not consistent.(" + str(read_1_count) + " =/= " + str(read_2_count))
    
    #return summary_fn
    return {"read_1_fn":read_1_fn, "read_2_fn":read_2_fn, "read_1_count":read_1_count, "read_2_count":read_2_count, "insert_size_mean":insert_size_mean, "insert_size_sd":insert_size_sd, "max_read_length":max_read_len}



####### Assemble stage #########
# For using IDBA-UD, it is required to merge both pair-end files into one single files 
def merge_paired_end_seq(read_1_fn, read_2_fn, outdir="idba_ud", merged_read_fn=None):
    print_status("Merging mate files")
    
    if merged_read_fn is None:
        merged_read_fn = read_1_fn.replace("." + FASTQ_EXT, "") + read_2_fn.replace("." + FASTQ_EXT, "") + "." + FASTA_EXT
    
    cmd = IDBA_UD_HOME + "/bin/fq2fa --merge --filter " + read_1_fn + " " + read_2_fn + " " + merged_read_fn
    print_status("command: " + cmd)
   
    if not VERBOSE_ONLY:    
        os.system(cmd)
    
    if not os.path.isfile(merged_read_fn):
        print_status("Problem on merging pair-end mate files")
        return None
    
    return merged_read_fn

       
# IDBA-UD
def run_idba_ud(merged_read_fn, idba_ud_outdir="idba_ud", min_contig=1200, mink=20, maxk=80, step=10, num_threads=16):
    print_status("Initializing IDBA-UD")
    
    print_status("Provoking IDBA_UD")
    
    # Check sure that the constant kShortSequence in short_sequence.h is modified to some value larger than the default 128
    cmd = IDBA_UD_HOME + "/bin/idba_ud --read " + merged_read_fn + " -o " + idba_ud_outdir + " --mink " + str(mink) + " --maxk " + str(maxk) + " --num_threads " + str(num_threads) + " --min_contig " + str(min_contig) + " --pre_correction --step " + str(step)
    print_status("command: " + cmd)
    
    # Count the lapsed time
    start_time = time.time()
    
    if not VERBOSE_ONLY:
        os.system(cmd)

    print_status("Time used: " + str(time.time() - start_time) + " seconds")
    
    print_status("Checking if idba_ud is completed properly (" + idba_ud_outdir + "/end" + ")")
    
    #return True
    if os.path.isfile(idba_ud_outdir + "/end") and os.path.isfile(idba_ud_outdir + "/contig.fa"):
        return True
    else:
        return False


"""
 This routine verifies the 
"""
def postprocess_idba_ud(idba_ud_outdir="idba_ud"):
    print_status("Processing outputs from IDBA_UD stage")
    
    # If everything goes well, idba_ud produces a dummy file, "end"
    if not os.path.isfile(idba_ud_outdir + "/end"):
        raise OSError, "IDBA_UD stage is not completed properly."

    if os.path.isfile(idba_ud_outdir + "/contig.fa"):
        cmd = SCRIPTS_HOME + "/fasta_stat.pl " + idba_ud_outdir + "/contig.fa"
        contig_summary = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
        print_status("Contig Summary:")
        print_status("\t" + contig_summary)
        return idba_ud_outdir + "/contig.fa"
    else:
        raise OSError, "Contig file does not exist at " + idba_ud_outdir + " . IDBA_UD stage is not completed properly."
        return None


####### Binning #############
"""
 Check sure that bash, "source /home/siukinng/.bashrc", "set PERL5LIB=/usr/lib/perl5/site_perl/5.8.8:"$PERL5LIB" is included in the qsub script  
"""
def run_MaxBin(contig_fn, merged_read_fn, maxbin_outdir="MaxBin", maxbin_outprefix="contig", thread=16):
    print_status("Initializing MaxBin")  
    
    cmd = MAXBIN_HOME + "/run_MaxBin.pl -contig " + contig_fn + " -out " + maxbin_outdir + "/" + maxbin_outprefix + " -thread " + str(thread) + " -plotmarker -reads " + merged_read_fn
    print_status("command: " + cmd)
    
    if not os.path.exists(maxbin_outdir):
        print_status("Output directory, " + maxbin_outdir + ", does not exist, we will create it.")
        os.makedirs(maxbin_outdir)
    
    if not VERBOSE_ONLY:
        os.system(cmd)

    return maxbin_outdir


"""
 
"""
def postprocess_MaxBin(maxbin_outdir):
    print_status("Processing outputs from MaxBin")  
    
    

# ORFs prediction
def run_Prodigal(fna_infn, p_opt="meta", prodigal_outdir="Prodigal", out_prefix=None, outfn=None, faa_outfn=None, fna_outfn=None):
    print_status("Initializing Prodigal");
    
    #"~/tools/Prodigal/prodigal -p meta -i ../P3_contig.fa  -o P3_contig-Prodigal.out -a P3_contig-Prodigal.faa  -d P3_contig-Prodigal.fna"
    
    if out_prefix is None:
        if len(fna_infn.split("/")) > 0:
            tmp = fna_infn.split("/")
            out_prefix = tmp[len(tmp) - 1]
            out_prefix = (out_prefix[::-1].split(".", 1)[1])[::-1]

    if outfn is None:
        outfn = prodigal_outdir + "/" + out_prefix + ".prodigal.out"
    if faa_outfn is None:
        faa_outfn = prodigal_outdir + "/" + out_prefix + ".prodigal." + FASTA_PROT_EXT
    if fna_outfn is None:
        fna_outfn = prodigal_outdir + "/" + out_prefix + ".prodigal." + FASTA_DNA_EXT
  
    cmd = PRODIGAL_HOME + "/prodigal -p " + p_opt + " -i " + fna_infn + " -a " + faa_outfn + " -d " + fna_outfn + " -o " + outfn
    print_status("command: " + cmd)
    
    if not os.path.exists(prodigal_outdir):
        print_status("Output directory, " + prodigal_outdir + ", does not exist, we will create it.")
        os.makedirs(prodigal_outdir)
    
    if not VERBOSE_ONLY:
        os.system(cmd)
    #os.system(cmd)
    
    #return prodigal_outdir
    return {"source_fn":fna_infn, "outdir":prodigal_outdir, "prodigal_outfn":outfn, "protein_outfn":faa_outfn, "nucleotide_outfn":fna_outfn}


# Contig binning using ESOM 
def run_ESOM():
    print_status("Initializing ESOM")    
    


####### Functional annotation stage #########
# BLASTN
def blast(query_fn, db_fn, outdir=".", outfn=None, outfmt=6, num_threads=16, evalue=0.01, blast_program="blastn"):
    print_status("Initializing for " + blast_program)

    # Name the outfile
    if outfn is None:
        subject_id = db_fn
        
        if len(subject_id.split("/")) > 1:
            subject_id = subject_id.split("/")[len(subject_id.split("/")) - 1]
        
        query_id = query_fn
        if len(query_fn.split("/")) > 0:
            tmp = query_id.split("/")
            query_id = tmp[len(tmp) - 1]
            query_id = (query_id[::-1].split(".", 1)[1])[::-1]
                  
        outfn = query_id + "-" + subject_id + ".bla"

    # We have to check if the output directory existed, if not, we make it
    if outdir != ".":
        if not os.path.exists(outdir):
            print_status("Output directory, " + outdir + ", does not exist, we will create it.")
            os.makedirs(outdir)     


    cmd = BLAST_HOME + "/bin/" + blast_program + " -db " + db_fn + " -outfmt " + str(outfmt) + " -num_threads " + str(num_threads) + " -evalue " + str(evalue) + " -query " + query_fn + " -out " + outdir + "/" + outfn
    print_status("command: " + cmd)
      
    #if not VERBOSE_ONLY:
    #    os.system(cmd)
    os.system(cmd)
    
    return outfn
 
    
def blastp(query_fn, db_fn, outdir=".", outfn=None, outfmt=6, num_threads=16):
    # blastp -query ../../Prodigal/contig.fa.prodigal.faa -db ~/db/Markers/CAZy/CAZy_id.lst.retrieved.faa -outfmt 6 -out contig.fa.prodigal-CAZy_id.lst.retrieved.bla  -num_threads 16
    blast(query_fn, db_fn, outdir, outfn, outfmt, num_threads, blast_program="blastp")
    
def blastn(query_fn, db_fn, outdir=".",outfn=None, outfmt=6, num_threads=16):
    # blastp -query ../../Prodigal/contig.fa.prodigal.faa -db ~/db/Markers/CAZy/CAZy_id.lst.retrieved.faa -outfmt 6 -out contig.fa.prodigal-CAZy_id.lst.retrieved.bla  -num_threads 16
    blast(query_fn, db_fn, outdir, outfn, outfmt, num_threads, blast_program="blastn")
    
    
        
def blastx():
    print_status("Initializing blastx()")    

"""
 Based on given threshold values, this routine will invoke filter_blast_res.py.
"""
def filter_blast_results(blast_fn, identity=80.0, length=100, subject_id_desc_fn=None, nr_query=True): 
    print_status("Filtering BLAST result") 
    print_status("Filtering thresholds: " + " Identity=" + str(identity) + " Alignment Length=" + str(length)) 
    
    

"""
 If EMIRGE is installed, please check emirge.py can run properly. In case it complains missing _emirge module, 
 you have to put "sys.path.append('path_to__emirge.so')" into _emirge.py before "import _emirge"
 
 If the sequencing reads were generated by an platform other than Illumina, please set platform_illumina to False
"""
def running_EMIRGE(read_1_fn, read_2_fn, DIR=MARKER_OUTDIR+"/EMIRGE", ssu_fn=EMIRGE_16S_DB, ssu_bowtie_fn=EMIRGE_16S_DB, platform_illumina=True):
    print_status("Initializing EMIRGE")
 
    read_summary = summarize_statistics(read_1_fn, read_2_fn)
    
    ## {"read_1_fn":read_1_fn, "read_2_fn":read_2_fn, "read_1_count":read_1_count, "read_2_count":read_2_count, "insert_size_mean":insert_size_mean, "insert_size_sd":insert_size_sd, "max_read_length":max_read_len}

    print_status("< Statistics of Paired-end Short Reads >")
    print_status("   Mean Insert Size: " + str(read_summary["insert_size_mean"]) + " bp")
    print_status("   Insert Size SD: " + str(read_summary["insert_size_sd"]) + " bp")
    print_status("   Maximum Read Length: " + str(read_summary["max_read_length"]) + " bp")
    
    cmd = EMIRGE_HOME + "/bin/emirge.py " + DIR + " -1 " + read_1_fn + " -2 " + read_2_fn + " -f " + ssu_fn + " -b " + ssu_bowtie_fn + " -l " + str(read_summary["max_read_length"]) + " -i " + str(int(math.ceil(read_summary["insert_size_mean"]))) + " -s " + str(int(math.ceil(read_summary["insert_size_sd"])))
    
    if platform_illumina:
        cmd = cmd + " --phred33"
    
    print_status("command: " + cmd)
      
    #if not VERBOSE_ONLY:
    #    os.system(cmd)  
    os.system(cmd)
    
 
"""
 
"""
def running_HMMER(query_fn, hmm_fn):
    print_status("Initializing HMMER3")
    
    
    
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



# Assert 
def assert_process(fn_to_be_asserted):
    print_status("assert_process()")


# 
def import_settings():
    print_status("importing settings")


def log(log_msg):
    print_status(log_msg)


# Print status
def print_status(msg):
    caller_name = inspect.stack()[1][3]
    msg = "[ " + caller_name + " ] " + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "  "  + msg
    
    print msg
    
    if LOG:
        LOG_FILE.write(msg + "\n")
    

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])