#!/share/apps/Python-2.7.4/bin/python

# 2014 SKWoolf bcnskaa AT gmail DOT com

"""
 This pipeline script is in a highly integrated and automatic manner designed for analyzing paired-end metagenomic data generated
 by next-generated sequencing platforms. After installing all dependencies (first-time only), the only input from end-user is two
 paired-end read files. Currently, this pipeline is tested under Linux environment only.
 
 The following dependencies are needed for this pipeline:
     1. BBMap
     2. NCBI BLAST+
     3. BOWTIE
     4. EMIRGE
     5. ESOM
     6. Fastq-MCF
     7. FastQC
     8. HMMER3
     9. IDBA_UD
    10. MaxBin
    11. Perl (bioperl)
    12. Prodigal
    13. Python (biopython)
    14. R (ggplot2, reshape2, ??)
    15. Samtool
    16. Seqtk
    17. Tetra-ESOM

 In addition to the above dependencies, this pipeline is also expected the following databases available:
     1. NCBI bacterial genomes (genomic + protein sequences)
     2. NCBI 16s ribosomal database
     3. CAZy database
     4. dbCAN
     5. EMIRGE 16s database
     6. (Single Copy Marker Gene database) specI 
    
  
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
import numpy
#from __future__ import print_function

# To run this script, path to biopython libraries has to be included in PYTHONPATH
from Bio import SeqIO
from Bio import SeqUtils
from Bio import Seq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW


####################### Apps Path #######################
# Global variables
HOME = "/home/siukinng"
TOOLS_HOME = HOME + "/tools"
SCRIPTS_HOME = TOOLS_HOME + "/scripts"

# Paths to individual packages' home directories
BBMAP_HOME = TOOLS_HOME + "/BBMap"
BLAST_HOME = TOOLS_HOME + "/blast"
BOWTIE2_HOME = TOOLS_HOME + "/"
EMIRGE_HOME = TOOLS_HOME + "/EMIRGE"
ESOM_HOME = TOOLS_HOME + "/ESOM"
FASTQ_MCF_HOME = TOOLS_HOME + "/fastq-mcf"
FASTQC_HOME = TOOLS_HOME + "/FastQC"
HMMER_HOME = TOOLS_HOME + "/hmmer"
IDBA_UD_HOME = TOOLS_HOME + "/idba_ud"
MAXBIN_HOME = TOOLS_HOME + "/MaxBin"
PRODIGAL_HOME = TOOLS_HOME + "/Prodigal"
SEQTK_HOME = TOOLS_HOME + "/seqtk"
SPECI_HOME = TOOLS_HOME + "/specI"
TETRA_ESOM_HOME = TOOLS_HOME + "/tetra-ESOM"
TRNA_SCAN_HOME = TOOLS_HOME + "/tRNAscan-SE"


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
DBCAN_HOME = DB_HOME + "/Markers/dbCAN"
DBCAN_HMM = DBCAN_HOME + "/dbCAN-fam-HMMs.txt.v3.txt"  # CAZy HMM profile 

# Folder contains all prokaryotic genomes available in NCBI
NCBI_BACTERIAL_GENEOMES_DB = DB_HOME + "/BacterialDB/all_fna"


LOG = True
LOG_FN= "pipeline.log"
LOG_FILE = open(LOG_FN, "w")

# Output directories
MARKER_OUTDIR = "Markers"
RESULT_OUTDIR = "Reports"

HMMER_OUTDIR = MARKER_OUTDIR + "/HMMER"

VERBOSE_ONLY = False


######### Main subroutine 
def main(argv):

    #if LOG:
    #    LOG_FILE = open(LOG_FN, "w")
    
    # Get current working directory
    CUR_WD = os.getcwd()
    print_status("Current Working Directory: " + CUR_WD)
    
    read_1_fn = argv[0]
    read_2_fn = argv[1]
    
    if len(argv) == 3:
        run_id = argv[2]
    else:
        run_id = read_1_fn.split("_1.",1)[0]

    print_status("Output Prefix = " + run_id)    
    
    # Preprocess raw reads and generate quality reports 
    read_1_processed_fn = preprocess(read_1_fn)
    read_2_processed_fn = preprocess(read_2_fn)
    
    
    # Filter read files
    read_fns = run_fastq_mcf(read_1_processed_fn, read_2_processed_fn, ADAPTOR_SEQS)
    
    # Go to hell if thing gets nasty after filtering
    if len(read_fns) != 2:
        print_status("Problem on calling fastq_mcf results.")
        raise OSError, "Problem on processing fastq_mcf outputs, abort now."
    
    # Generate a summary of read data
    summarize_statistics(read_fns[0], read_fns[1])
    
    # We need to merge paired-end files (fastq) into one single file (fasta) before calling IDBA-UD
    merged_read_fn = merge_paired_end_seq(read_fns[0], read_fns[1])
    
    # Die if anything goes wrong 
    if merged_read_fn is None:
        print_status("Problem on reading from merged_read_fn.")
        raise OSError, "Problem at merging reads files, abort now."
    
    # Alright, we are about to do assembly 
    if not run_idba_ud(merged_read_fn, min_contig=800):
        print_status("Problem at completing IDBA_UD stage")
        raise OSError, "Problem at IDBA_UD stage, abort now."
    
    # Double check and clean up output files generated by assembling steps
    contig_fn = postprocess_idba_ud()  
    
    # Sorry, I am kicking you off as something goes wrong
    if contig_fn is None:
        print_status("Unable to locate contig.fa")
        raise OSError, "Contig file is not avaiable from IDBA_UD, abort now."
    
    # Initial stage of inferring candidate taxonomical identities in the current dataset
    blast_outfn = run_id+"-ncbi16s.bla"
    blastp(contig_fn, NCBI16S_DB, outdir=MARKER_OUTDIR, outfn=blast_outfn)
    

    # Great, we are moving to do binning
    maxbin_outdir = run_MaxBin(contig_fn, merged_read_fn)
    
    # Assure, filter and manage output files from MaxBin
    postprocess_MaxBin(maxbin_outdir)
    
    
    # Do ESOM classification
    run_tetraESOM()
    
    
    # 
    bin_group_fasta_fns = glob.glob(maxbin_outdir + "/*.fasta")
    if len(bin_group_fasta_fns) == 0:
        print_status("No bin group is found in " + maxbin_outdir)
        raise OSError, "No bin group is found in " + maxbin_outdir
    
    
    # For each binned group
    for bin_group_fasta_fn in bin_group_fasta_fns:
        bin_group_id = (bin_group_fasta_fn.replace(maxbin_outdir + "/", "")).replace(".fasta", "")
        prodigal_info = run_Prodigal(bin_group_fasta_fn, out_prefix=bin_group_id)

 
    # For all combined result  
    prodigal_info = run_Prodigal(contig_fn)
    # {"source_fn":fna_infn, "outdir":prodigal_outdir, "prodigal_outfn":outfn, "protein_outfn":faa_outfn, "nucleotide_outfn":fna_outfn}


    # blastp -query ../../Prodigal/contig.fa.prodigal.faa -db ~/db/Markers/CAZy/CAZy_id.lst.retrieved.faa -outfmt 6 -out contig.fa.prodigal-CAZy_id.lst.retrieved.bla  -num_threads 16
    if not os.path.isfile(prodigal_info["protein_outfn"]):
        blastp(prodigal_info["protein_outfn"], CAZY_DB, outdir=MARKER_OUTDIR + "/CAZy")
    else:
        print_status("Warning: Predictive protein sequence is not available from Prodigal, step of mapping to CAZy database is skipped.")
    
    
    # Scan predicted protein sequences for CAZy domains 
    if(running_HMMER_search(DBCAN_HMM, prodigal_info["protein_outfn"])):
        # Parse the output files from HMMER3
        post_process_HMMER_search()
    else:
        print_status("Unable to process outputs from HMMER search")
    
     
    # Blast any 16s sequence fragments existed in newly assembled contig sequences
    blastn(contig_fn, NCBI16S_DB, outdir=MARKER_OUTDIR + "/ncbi16s")
    
    
    # Based on the 16s blast result, we construct a reference genome library
    prepare_reference_genome()  
    
    # 
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


  
# run_FastQC   
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



# Seqtk 
def run_seqtk(read_fn, trimmed_read_fn):
    print_status("Processing " + read_fn)    
    
    cmd = SEQTK_HOME + "/seqtk trimfq " + read_fn + " > " + trimmed_read_fn
    print_status("command: " + cmd)
    
    if not VERBOSE_ONLY:
        os.system(cmd)



###############################################
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
  Generate some statistics about the read data 
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
        merged_read_fn = read_1_fn.replace("." + FASTQ_EXT, "") + "-" + read_2_fn.replace("." + FASTQ_EXT, "") + "." + FASTA_EXT
    
    cmd = IDBA_UD_HOME + "/bin/fq2fa --merge --filter " + read_1_fn + " " + read_2_fn + " " + merged_read_fn
    print_status("command: " + cmd)
   
    if not VERBOSE_ONLY:    
        os.system(cmd)
    
    if not os.path.isfile(merged_read_fn):
        print_status("Problem on merging pair-end mate files")
        return None
    
    return merged_read_fn



###### 
# IDBA-UD
"""
 
"""
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
 This routine verifies if IDBA_UD completed successfully, and generates 
 statistics on assembled contig sequences.
"""
def postprocess_idba_ud(idba_ud_outdir="idba_ud", report_outdir=RESULT_OUTDIR):
    print_status("Processing outputs from IDBA_UD stage")
    
    # If everything goes well, idba_ud produces a dummy file, "end"
    if not os.path.isfile(idba_ud_outdir + "/end"):
        raise OSError, "IDBA_UD stage is not completed properly."

    if os.path.isfile(idba_ud_outdir + "/contig.fa"):
        cmd = SCRIPTS_HOME + "/fasta_stat.pl " + idba_ud_outdir + "/contig.fa"
        contig_summary = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
        print_status("Contig Summary:")
        print_status("\t" + contig_summary)
        
        # Summarize statistics of every sequence inside contig.fa in a tabular format to a file under RESULT_OUTDIR
        #cmd = SCRIPTS_HOME + "/fasta_len_spectrum.pl "
        
        # Compress large interim files
        cmd = "bzip2 -9 " + idba_ud_outdir + "/contig-* " + idba_ud_outdir + "/align-* " + idba_ud_outdir + "/graph-* " + idba_ud_outdir + "/local-* " + idba_ud_outdir + "/kmer"   
        print_status("Compressing IDBA_UD interim files")
        
        if not VERBOSE_ONLY:
            os.system(cmd)
        
        
        return idba_ud_outdir + "/contig.fa"
    else:
        raise OSError, "Contig file does not exist at " + idba_ud_outdir + " . IDBA_UD stage is not completed properly."
        return None




####### Binning ###########
"""
 Generate a summary of fasta sequences
"""
def summarize_fasta(fasta_fn, outtbl_fn=None):
    print_status("Summarizing " + fasta_fn) 
    
    # if outfile does not specify, we name it
    if outtbl_fn is None:
        outtbl_fn = fasta_fn + ".summary"

    seq_n = 0
    
    # Go through every fasta sequences
    with open(fasta_fn, "r") as infile, open(outtbl_fn, "w") as outfile:
        outfile.write("ID\tLEN\tGC\n")
        for seq in SeqIO.parse(infile, "fasta"):
            print("Processing " + str(seq.id))
            txt = seq.id + "\t" + str(len(seq.seq)) + "\t" + str(SeqUtils.GC(seq.seq)) + "\n"
            outfile.write(txt)
            seq_n += 1
            
            
    print_status("Sequence processed: " + str(seq_n))



"""
 Check sure that bash, "source /home/siukinng/.bashrc", "set PERL5LIB=/usr/lib/perl5/site_perl/5.8.8:"$PERL5LIB" is included in the qsub script  
"""
def run_MaxBin(contig_fn, merged_read_fn, maxbin_outdir="MaxBin", maxbin_outprefix="contig", thread=16):
    print_status("Initializing MaxBin")  
    
    cmd = MAXBIN_HOME + "/run_MaxBin.pl -contig " + contig_fn + " -out " + maxbin_outdir + "/" + maxbin_outprefix + " -thread " + str(thread) + " -plotmarker -reads " + merged_read_fn
    print_status("command: " + cmd)
    
    # 
    if not os.path.exists(maxbin_outdir):
        print_status("Output directory, " + maxbin_outdir + ", does not exist, we will create it.")
        os.makedirs(maxbin_outdir)
    
    # 
    if not VERBOSE_ONLY:
        os.system(cmd)

    return maxbin_outdir



"""
 Filter binning groups
 maxbin_outdir: Directory containing files generated by MaxBin
 min_total_contig_len: Minimum value of accumulative length of a bin group
 marker_gene_yield: Number of marker gene observed in a bin group
"""
def postprocess_MaxBin(maxbin_outdir, min_total_contig_len=1000000, marker_gene_yield=50.0):
    print_status("Processing outputs from MaxBin")  
    
    maxbin_summary_fn = glob.glob(maxbin_outdir + "/*.summary")
    
    if len(maxbin_summary_fn) != 1:
        print_state("Unable to process maxbin results at " + maxbin_outdir)
        return False
    
    return True

    
    
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
def run_tetraESOM():
    print_status("Initializing ESOM")    
    
    
    
    
####### Taxonomical assingment stage #########
"""
for f in ../*.fasta; do
    id=${f/.fasta/};
    id=${id/..\/};
    cmd="~/tools/Prodigal/prodigal -p meta -i ../"$id".fasta -a "$id"-prodigal.faa -d "$id"-prodigal.fna -o "$id"-prodigal.out";
    echo $cmd;
    eval $cmd;
    cmd="python ~/db/Markers/specI/specI.py "$id"-prodigal.fna "$id"-prodigal.faa 4 "$id"-prodigal /home/siukinng/samples/lab/GZ-xyl_S2/MaxBin/specI";
    echo $cmd;     
    eval $cmd;
done

    glsearch36 are needed for using specI
    http://faculty.virginia.edu/wrpearson/fasta/fasta36/
"""
def run_specI(gene_fn, protein_fn, base_dir, runName=None, num_of_thread=16):
    print_status("Initializing SpecI for assigning taxonomical identity to " + gene_fn)
    
    # No name provided, then we create it
    if runName == None:
        runName = (gene_fn[::-1].split(".", 1))[1]
        runName = runName[::-1]
    
    # Assume a result generated by the specI
    specI_result_fn = base_dir + "/" + runName + ".results"
    print_status("Results will be exported to " + specI_result_fn)
    
        
    # Formulate the command
    cmd = "python " + SPECI_HOME + "/specI.py " + gene_fn + " " + protein_fn + " " + str(num_of_thread) + " " + runName + " " + base_dir
    print_status("command: " + cmd)
    
    # Execute the command 
    if not VERBOSE_ONLY:
        os.system(cmd)
    
    # Check if a results outfile is generated by SpecI
    if not os.path.isfile(specI_result_fn):
        print_status("No result is avialable for " + gene_fn)
        return None

    # Read the result
    with open(specI_result_fn) as input:
        lines = input.read() #splitlines('\t')
        lines = lines.splitlines()
    #IN.close()
     
    # Prepare a dict
    specI_results = {}
    for result in lines:
        result = result.split("\t")
        
        if(result[3] == "-"):
            specI_results["No Similar Found"] = 0.0
        else:
            specI_results[result[3]] = result[4]  
    # 
    return specI_results


     
"""
    This routine will search 16s sequence from the input fna sequence file. The identified 16s sequence will then
    blast against NCBI web blast. Result will be parsed, and subject id of top hit will be returned. 
"""
def run_16s_mapping(fna_fn, outdir=None, outprefix=None, blast_bitscore_threshold=200, wwwblast_evalue=1e-10, wwwblast_hitlist_size=5):
    print_status("Searching 16s sequence from " + fna_fn)

    if outdir is None:
        outdir = "."  
    
    if outprefix is None:
        outprefix = fna_fn[::-1].split("/", 1)
        outprefix = outprefix[0][::-1]

    # Blast the fna against 16s database
    blast_outfn = outprefix + "-ncbi16s.bla"
    print_status("Blast result will be exported to " + blast_outfn)

    
    blastn(query_fn=fna_fn, db_fn=NCBI16S_DB, outdir=outdir, outfn=blast_outfn, best_hit_score_edge=None, best_hit_overhang=None, perc_identity=None, max_target_seqs=None, evalue=None)
    
    # Verify if the blast result is available
    if not os.path.isfile(blast_outfn):
        return None
        
    # Check out the hit with highest bit-score
    mapped_16s = []
    bres = open(blast_outfn).read().splitlines()
    if len(bres) > 0:
        items = bres[0].split("\t")
        bit_score = float(items[11])
        if bit_score > blast_bitscore_threshold:
            print outprefix + "=" + items[0] + ": " + str(bit_score)
            # query, q_spos, q_epos, subject, s_spos, s_epos, length, identity, bit-score 
            mapped_16s = [items[0], int(items[6]), int(items[7]), items[1], int(items[8]), int(items[9]), int(items[3]), float(items[4]), bit_score] 
    
    #return mapped_16s
    
    if len(mapped_16s) == 0:
        print_status("No 16s sequence found in " + blast_outfn)
        return None
    
    # Pick the 16s sequence from the contig sequences
    seq_16s = pick_seq(fna_fn, mapped_16s[0], mapped_16s[1], mapped_16s[2])
    
    # Problem on retrieving the 16s sequence
    if seq_16s is None:
        print_status("Cannot extract the 16s sequence from " + fna_fn)
        return None
    
    # Reverse complement if it is in minus strand
    if mapped_16s[5] - mapped_16s[4] < 0:
        seq_16s = str(Seq(seq_16s).reverse_complement())

    # Submit the sequence to NCBI wwwblast and check if there is any hit.
    #handler = NCBIWWW.qblast("blastn", "nr", seq_16s, expect=wwwblast_evalue, hitlist_size=wwwblast_hitlist_size)
    records = run_ncbi_wwwblast(seq=seq_16s, database="nr", blast_program="blastn", expect=wwwblast_evalue, hitlist_size=wwwblast_hitlist_size)
    
    #species = None
    candidate_species = []
    if records is not None:
        if len(records[0].alignments) > 0:    
            for alignment in records[0].alignments:
                if alignment.hsps[0].score > blast_bitscore_threshold:
                    species = str((alignment.title).split("|")[4])
                    species = species.replace("16S ribosomal RNA gene, complete sequence", "")
                    species = species.replace("16S rRNA gene, ", "")
                    species = species.replace("gene for 16S rRNA, partial sequence, ", "")
                    species = species.replace("16S ribosomal RNA gene, partial sequence", "")
                    species = species.strip()
                    candidate_species.append(species)

#             alignment = records[0].alignments[0]
#              
#             if alignment.hsps[0].score > blast_bitscore_threshold:
#                 species = str((alignment.title).split("|")[4])
#                 species = species.replace("16S ribosomal RNA gene, complete sequence", "")
#                 species = species.replace("16S rRNA gene, ", "")
#                 species = species.replace("gene for 16S rRNA, partial sequence, ", "")
#                 species = species.replace("16S ribosomal RNA gene, partial sequence", "")
#                 species = species.strip()
                
    return candidate_species
    


"""
 Invoke a NCBI web blast search
"""
def run_ncbi_wwwblast(seq, database, blast_program="blastn", expect=10, hitlist_size=50):
    print_status("Connecting to NCBI Web BLAST: " + blast_program + " e-value=" + str(expect) + " hitlist_size=" + str(hitlist_size))
    
    result_handler = NCBIWWW.qblast(program=blast_program, database=database, sequence=seq, expect=expect, hitlist_size=hitlist_size)
    
    # Store the wwwblast results to local memory
    records = list(NCBIXML.parse(result_handler))
    
    return records



"""
    
""" 
def pick_seq(seq_fn, seq_id, spos=-1, epos=-1, format="fasta"):
    print_status("Picking " + seq_id + " from " + seq_fn)

    # Build an index of the sequences available in the file
    seq_idx = SeqIO.index(seq_fn, format)
    
    seq = None
    # Check if the seq_id is existed in the index
    if seq_id in seq_idx.keys():
        if spos != -1 and epos != -1:
            seq = str(seq_idx[seq_id][spos:epos].seq)
        else:
            seq = str(seq_idx[seq_id].seq)
             
    return seq




####### Functional annotation stage #########
# BLASTN
# -evalue 1e-10 -best_hit_score_edge 0.05 -best_hit_overhang 0.25 -perc_identity 80 -max_target_seqs 2
def blast(query_fn, db_fn, outdir=".", outfn=None, outfmt=6, num_threads=16, evalue=1e-10, best_hit_score_edge=0.05, best_hit_overhang=0.25, perc_identity=80, max_target_seqs=2, blast_program="blastn"):
    print_status("Initializing " + blast_program)

    # Name the outfile
    if outfn is None:
        subject_id = db_fn
        
        if len(subject_id.split("/")) > 1:
            subject_id = subject_id.split("/")[len(subject_id.split("/")) - 1]
        
        query_id = query_fn
        if len(query_fn.split("/")) > 0:
            tmp = query_id.split("/")
            query_id = tmp[len(tmp) - 1]
            query_id = (query_id[::-1].split(".", 1)[1])[::-1]  # Extract the first part
                  
        outfn = query_id + "-" + subject_id + ".bla"

    # We have to check if the output directory existed, if not, we make it
    if outdir != ".":
        if not os.path.exists(outdir):
            print_status("Output directory, " + outdir + ", does not exist, we will create it.")
            os.makedirs(outdir)     
            
    cmd = BLAST_HOME + "/bin/" + blast_program + " -db " + db_fn + " -outfmt " + str(outfmt) + " -num_threads " + str(num_threads)
    
    if evalue is not None:
        cmd = cmd + " -evalue " + str(evalue)
        
    if best_hit_score_edge is not None:
        cmd = cmd + " -best_hit_score_edge " + str(best_hit_score_edge)
        
    if best_hit_overhang is not None: 
        cmd = cmd + " -best_hit_overhang " + str(best_hit_overhang)
        
    if perc_identity is not None: 
        cmd = cmd + " -perc_identity " + str(perc_identity)  
        
    if max_target_seqs is not None: 
        cmd = cmd + " -max_target_seqs " + str(max_target_seqs)     
    
    cmd = cmd + " -query " + query_fn + " -out " + outdir + "/" + outfn
    
    print_status("command: " + cmd)
      
    #if not VERBOSE_ONLY:
    #    os.system(cmd)
    os.system(cmd)
    
    return outfn




"""
"""    
def blastp(query_fn, db_fn, outdir=".", outfn=None, outfmt=6, num_threads=16, best_hit_score_edge=0.05, best_hit_overhang=0.25, perc_identity=80, max_target_seqs=2, evalue=1e-10):
    # blastp -query ../../Prodigal/contig.fa.prodigal.faa -db ~/db/Markers/CAZy/CAZy_id.lst.retrieved.faa -outfmt 6 -out contig.fa.prodigal-CAZy_id.lst.retrieved.bla  -num_threads 16
    blast(query_fn=query_fn, db_fn=db_fn, outdir=outdir, outfn=outfn, outfmt=outfmt, num_threads=num_threads, best_hit_score_edge=best_hit_score_edge, best_hit_overhang=best_hit_overhang, perc_identity=perc_identity, max_target_seqs=max_target_seqs, evalue=evalue, blast_program="blastp")
    

    
"""
 
"""
def blastn(query_fn, db_fn, outdir=".",outfn=None, outfmt=6, num_threads=16, best_hit_score_edge=0.05, best_hit_overhang=0.25, perc_identity=80, max_target_seqs=2, evalue=1e-10):
    # blastp -query ../../Prodigal/contig.fa.prodigal.faa -db ~/db/Markers/CAZy/CAZy_id.lst.retrieved.faa -outfmt 6 -out contig.fa.prodigal-CAZy_id.lst.retrieved.bla  -num_threads 16
    blast(query_fn=query_fn, db_fn=db_fn, outdir=outdir, outfn=outfn, outfmt=outfmt, num_threads=num_threads, best_hit_score_edge=best_hit_score_edge, best_hit_overhang=best_hit_overhang, perc_identity=perc_identity, max_target_seqs=max_target_seqs, evalue=evalue, blast_program="blastn")
    

 
"""
 
"""       
def blastx():
    print_status("Initializing blastx()")
    


"""
 Based on given threshold values, this routine will invoke filter_blast_res.py.
"""
def filter_blast_results(blast_fn, identity=80.0, length=100, subject_id_desc_fn=None, nr_query=True): 
    print_status("Filtering BLAST result") 
    print_status("Filtering thresholds: " + " Identity=" + str(identity) + " Alignment Length=" + str(length)) 
    
    cmd = SCRIPTS_HOME + "/filter_blast_res.py -i " + blast_fn + " -p " + str(identity) + " -l " + str(length)
    # If the option -q is needed
    if nr_query:
        print_status("Best hit will be exported for each query id") 
        cmd = cmd + " -q "
        
    # If a list of subject ids is provided
    if subject_id_desc_fn is not None:
        print_status("Subject ID will be prepared and exported.") 
        cmd = cmd + " -k " + subject_id_desc_fn + " -s"
    
    print_status("command: " + cmd)
      
    if not VERBOSE_ONLY:
        os.system(cmd)
        
    
   

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
    
    # Special care is needed for processing quality values of the data generated from Illumina platforms
    if platform_illumina:
        cmd = cmd + " --phred33"
    
    print_status("command: " + cmd)
      
    if not VERBOSE_ONLY:
        os.system(cmd)  
    


"""
 This routine scans protein sequences against protein HMM profiles
"""
def running_HMMER_scan(prot_faa_fn, hmm_profile_fn, outdir=HMMER_OUTDIR):
    print_status("Initializing HMMER3 hmmscan")
    

"""
 This routine processes output files from HMMER3 scan
  
"""
def postprocess_HMMER_scan(outdir=HMMER_OUTDIR, mean_posterior_prob=0.8):
    print_status("Parsing HMMER3 hmmscan outfiles")
    
    # Select HMM hits with mean posterior probability higher than the threshold
    # Estimate the domain number of each protein sequence
    # Concatenate all hit domains
    # 
    


# Generate a random sequence of length (default=100bp) for estimating the significance of HMMsearch results
def generate_random_seq(length=100, base = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z']):
    base_len = len(base) - 1
    seq = "".join([base[random.randint(0, base_len)] for i in range(1, length)])
    
    return seq



"""
 This routine searches protein HMM profiles against a protein sequence database. 
 ~/tools/hmmer/bin/hmmsearch 
 -o contig.fa.prodigal-dbCAN.hmm.out 
 -A contig.fa.prodigal-dbCAN.hmm.aln 
 --tblout contig.fa.prodigal-dbCAN.hmm.seq.tbl 
 --domtblout contig.fa.prodigal-dbCAN.hmm.dom.tbl 
 --pfamtblout contig.fa.prodigal-dbCAN.hmm.pfam.tbl 
 ~/db/Markers/dbCAN/dbCAN-fam-HMMs.txt.v3.txt contig.fa.prodigal.faa
 
"""
def running_HMMER_search(hmm_profile_fn, prot_faa_fn, outdir=HMMER_OUTDIR, outfn_prefix=None, tblout_fn=None, domtblout_fn=None, pfamtblout_fn=None, alignment_outfn=None):
    print_status("Initializing HMMER3 hmmsearch")
 
    if not os.path.isfile(hmm_profile_fn):
        print_status("Unable to read from " + hmm_profile_fn + ", hmmsearch is skipped.")
        return None

    # Number of profiles in the file
    cmd = "cat " + hmm_profile_fn + " | grep -c \"NAME\""
    hmm_profile_count = int(subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read())
    
    # Length of profiles
    cmd = "cat " + hmm_profile_fn + " | grep -c \"LENG\" | cut -d ' ' -f3"
    hmm_lens = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    hmm_lens = [int(len.strip()) for len in hmm_lens.splitlines()]
    
    # Number of sequences recruited in each profile
    cmd = "cat " + hmm_profile_fn + " | grep \"NSEQ\" | cut -d ' ' -f3"
    hmm_nseqs = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
    hmm_nseqs = [int(nseq.strip()) for nseq in hmm_nseqs.splitlines()] 
    
    
    print_status("< Statistics of Query HMM Profiles [" + hmm_profile_fn + "] >")
    print_status("   Number of HMM profiles: " + str(hmm_profile_count))
    print_status("   Length of Profiles (residue): Mean=" + str(numpy.mean(hmm_lens)) + " Min=" + str(min(hmm_lens)) + " Max=" + str(max(hmm_lens)))
    print_status("   NSEQ in Profile: Mean=" + str(numpy.mean(hmm_nseqs)) + " Min=" + str(min(hmm_nseqs)) + " Max=" + str(max(hmm_nseqs)))
 

    if outfn_prefix is None:
        outfn_prefix = prot_faa_fn.replace("."+FASTA_PROT_EXT, "")
   
    if alignment_outfn is None:
        alignment_outfn = outdir + "/" + outfn_prefix + ".aln"        
    if tblout_fn is None:
        tblout_fn = outdir + "/" + outfn_prefix + ".tbl"
    if domtblout_fn is None:
        domtblout_fn = outdir + "/" + outfn_prefix + ".dom.tbl"
    if pfamtblout_fn is None:
        pfamtblout_fn = outdir + "/" + outfn_prefix + ".pfam.tbl"
    
    outfn = outdir + "/" + outfn_prefix + ".out"         
    
    cmd = HMMER_HOME + "/bin/hmmsearch " + " -o " + outfn + " -A " + alignment_outfn + " --tblout " + tblout_fn + " --domtblout " + domtblout_fn + " --pfamtblout " + pfamtblout_fn + " " + hmm_profile_fn + " " + prot_faa_fn
    print_status("command: " + cmd)
      
    #if not VERBOSE_ONLY:
    #    os.system(cmd)  
    os.system(cmd)

    # Make sure 
    if assert_proc(outfn) and assert_proc(alignment_outfn) and assert_proc(tblout_fn) and assert_proc(domtblout_fn) and assert_proc(pfamtblout_fn):
        print_status("HMMsearch completed!")
        return True
    else:
        return False


# Return the length of overlapping between two regions
def getOverlap(x, y):
    #print "x[1]=", x[1], ", y[1]=", y[1], ", x[0]=", x[0], ", y[0]=", y[0]
    return max(0, min(x[1], y[1]) - max(x[0], y[0]))



"""
 This routine processes output files from HMMER3 scan
 cat contig.fa.prodigal-dbCAN.hmm.dom.tbl | grep -v "^#" | sed 's/\s\s*/ /g' | cut -d ' ' -f22
"""
def postprocess_HMMER_search(hmm_dir=HMMER_OUTDIR, mean_posterior_prob=0.8, hmm_score_threshold=60.0, dom_overlapping_threshold=20):
    print_status("Parsing HMMER3 hmmsearch outfiles")
    
    # Check path exists
    if not os.path.exists(hmm_dir):
        print_status("Unable to read from" + hmm_dir)
        return None   
    
    # 1. Consolidate the queries sharing a same query id and having mean posterior probability higher than a threshold (default=0.8)
    # 2. Summary statistics of domain counts
    # 3. Domain models of every query sequence


    # Import .dom.tbl outfile
    #line_n = 0
    #for file in glob.glob("*.dom.tbl"):
    
    file = glob.glob(hmm_dir + "/*.dom.tbl")
    if len(file) != 1:
        print_status("Does not find any .dom.tbl in the folder \"" + hmm_dir + "\"")
        return False
    
    file = file[0]    
    
    # ORF with HMM hits
    hmm_orf_dict = {}
    hmm_scores = []
    skipped_dom_n = 0
    processed_dom_n = 0
    discard_dom_n = 0
    with open(file, "r") as IN:
        for line in IN:
            
            #if line_n >= 3 and not line.startswtih("#"):
            if not line.startswith("#"):
                dom = line.split()
                  
                hmm_score = float(dom[7])
                hmm_scores.append(hmm_score)
                
                hmm_dom_score = float(dom[13])
                
                processed_dom_n += 1
                if hmm_score >= hmm_score_threshold:
                    tid = dom[0]
                    tlen = int(dom[2])
                    aln_spos = int(dom[17])
                    aln_epos = int(dom[18])
                    hmm_id = dom[3].replace(".hmm","")
                    hmm_len = int(dom[5])
                    hmm_spos = int(dom[15])
                    hmm_epos = int(dom[16])
                                  
                    # 
                    if tid not in hmm_orf_dict.keys():
                        #hmm_orf_dict[dom[0]] = ["".join([hmm_id, "@", aln_spos, "-", aln_epos, "=", dom[7]])]
                        hmm_orf_dict[tid] = [[tid, tlen, aln_spos, aln_epos, hmm_score, hmm_id, hmm_len, hmm_spos, hmm_epos, hmm_dom_score]]
                                                
                        #hmm_orf_dict[dom[0]] = [hmm_id]
                    else:
                        #hmm_orf_dict[dom[0]].append("".join([hmm_id, "@", aln_spos, "-", aln_epos, "=", dom[7]]))
                        hmm_orf_dict[tid].append([tid, tlen, aln_spos, aln_epos, hmm_score, hmm_id, hmm_len, hmm_spos, hmm_epos, hmm_dom_score])

                else:
                    #print dom[0], hmm_id, "=", dom[7], "skipped"
                    skipped_dom_n += 1
    
    
    # Check if there are overlapping HMM domains, keep the one with highest score    
    for k,v_arr in hmm_orf_dict.items():
        if len(v_arr) > 0:
            # Sorted by HMM dom score
            sorted_a = sorted(v_arr, key=lambda v: v[9], reverse=True)
            
            flags = [0 for i in range(len(sorted_a))]
            for idx, v in enumerate(sorted_a):
                if flags[idx] == 1:
                    continue
                
                for idx2, v2 in enumerate(sorted_a):
                    if idx2 > idx:
                        if getOverlap(v[2:4], v2[2:4]) > dom_overlapping_threshold:
                            flags[idx2] = 1
                            discard_dom_n += 1
                            #print v2, "discarded by", v, " overlap_len=", getOverlap(v[2:4], v[2:4])
            
            # Discard tagged items             
            sorted_a = [sorted_a[idx] for idx, flag in enumerate(flags) if flag == 0]
            # Sorted by position
            sorted_a = sorted(v_arr, key=lambda v: v[1], reverse=False)
            
            hmm_orf_dict[k] = sorted_a
    
    print_status("Processed HMM domains = " + str(processed_dom_n))   
    print_status( "Skipped HMM domains = "+str(skipped_dom_n))
    print_status( "Discarded HMM domains = "+str(discard_dom_n))
    
    return hmm_orf_dict  
    



# Mapping identified HMM domains to bin_groups
# Keys from hmm_orf_dict should be ORFs names generated by Prodigal. They are contig prefix + "_" + number.
# The contig prefixes are expected to be found in binning groups by MaxBin
# Get Current Dir: os.path.dirname(os.path.abspath("."))
def map_hmm2maxbin(hmm_orf_dict, maxbin_dir): 
    # Check path exists
    if not os.path.exists(maxbin_dir):
        print_status("Unable to read from \"" + maxbin_dir + "\"")
        return None


    ##################################
    # Extracting contig abundances
    contig_abunds = {}
    abund_fn = glob.glob(maxbin_dir + "/*.abund") 
    abund_fn = abund_fn[0]
    with open(abund_fn, "r") as IN:
        for idx, line in enumerate(IN):
            [contig_name, contig_abund] = (line.replace("\n", "")).split("\t")
            contig_abunds[contig_name] = float(contig_abund)
            
    
    ##################################
    #orf2contig_dict = {orf:(orf[::-1].split("_", 1))[1][::-1] for orf in hmm_orf_dict}
    contig2orf_dict = {}
    for orf in hmm_orf_dict:
        contig_id = (orf[::-1].split("_", 1))[1][::-1]
        if contig_id not in contig2orf_dict.keys():
            contig2orf_dict[contig_id] = [orf]
        else:
            contig2orf_dict[contig_id].append(orf)
    

    ##################################
    # Mapping predicted ORFs with HMM domains info to binned groups
    # Read the summary file
    summary_fn = glob.glob(maxbin_dir + "/*.summary")
    
    if len(summary_fn) != 1:
        print_state("Unable to read *.summary from"+maxbin_dir)
        return False
    
    summary_fn = summary_fn[0]
    
    # Contig summary
    contig_n = 0
    bin_groups = {}
    bin_group_n = 0
    mapped_dom_n = 0
    with open(summary_fn, "r") as IN:
        for idx, line in enumerate(IN):
            if idx > 0:
                items = (line.replace("\n", "")).split("\t")
                
                bin_group_n += 1
                bin_id = items[0].replace(".fasta", "")
               
                bin_fn = maxbin_dir + "/" + items[0]
                
                print_status("Importing from" + bin_fn)
                
                print_status("++++++++++++++++++++++++++++++++++++++++")
                print_status( "ID: " + bin_id)
                print_status( "++++++++++++++++++++++++++++++++++++++++")
                
                bin_groups[bin_id] = []
                
                # Get all the header of the current bin group
                with open(bin_fn, "r") as IN:
                    for line in IN:
                        if line.startswith(">"):
                            contig_name = (line.replace(">", "")).replace("\n", "")
                            
                            print(contig_name)
                            
                            if contig_name in contig2orf_dict.keys():
                                contig_n += 1
                                orf_names = contig2orf_dict[contig_name]
                                
                                for orf in orf_names:
                                    if orf in hmm_orf_dict.keys():
                                        mapped_dom_n += len(hmm_orf_dict[orf])
                                        bin_groups[bin_id].extend(hmm_orf_dict[orf])
                                
#                             if bin_id not in bin_groups.keys():
#                                 bin_groups[bin_id] = [contig_name]
#                             else:
#                                 bin_groups[bin_id].append(contig_name)
                                
#                             #bin_groups[bin_id].append((line.replace(">", "")).replace("\n", ""))
#                             contig_name = (line.replace(">", "")).replace("\n", "")
#                             if contig_name not in bin_groups.keys():
#                                 bin_groups[contig_name] = bin_id
#                                 #print (line.replace(">", "")).replace("\n", "")
#                                 contig_n += 1
#                             else:
#                                 print ("Duplicate found.("+contig_name+")")

    print_status("Bin Group# = " + str(bin_group_n) + "  Total Contig# = " + str(contig_n) + "  Mapped Domain# = " + str(mapped_dom_n))

    if bin_group_n == 0 or contig_n == 0:
        print_status("Something wrong with MaxBin data. Abort now.")
        return None
    else:
        
        bin_groups = {group: sorted(bin_groups[group], key=lambda v:v[4], reverse=True) for group in bin_groups.keys()}
        
        return bin_groups

  

"""
 Export the bin groups into a tabular file
"""
def generate_map_table(bin_groups, outfn_prefix):
    
    dom_list = {}
    #dom_list["__GROUP_NAME__"] = sorted(bin_groups.keys())
    dom_list_header = sorted(bin_groups.keys())
    # Construct a new list of dom ids
    for key in bin_groups.keys():
        group = bin_groups[key]
        for item in group:
            dom = item[5]
            if dom not in dom_list.keys():
                #dom_list[dom] = []
                dom_list[dom] = [0 for i in range(len(bin_groups))]
    
    
    # Populate the dom_list
    for idx, key in enumerate(dom_list_header):
    #for idx, key in enumerate(dom_list["__GROUP_NAME__"]):
        group = bin_groups[key]
        
        # Create a domain list with unique orf id
        unique_group = list(set([v[5]+"<>"+v[0] for v in group]))
        
        
        for items in unique_group:
            dom = (items.split("<>"))[0]
            dom_list[dom][idx] += 1
            
    # Export the list
    outfn = outfn_prefix + ".tbl"
    print_status("Exporting to "+outfn)
    out = open(outfn, "w")
    out.write("\t".join( ["HMM_DOM", ("\t".join(dom_list_header)), "\n"]))
              
    #print >>out, "\t".join( ["HMM_DOM", ("\t".join(dom_list_header)) ])
    for key in sorted(dom_list.keys()):
        items = dom_list[key]
        out.write("\t".join([key, "\t".join([str(v) for v in items]), "\n"]))
        #print >>out, "\t".join([key, "\t".join([str(v) for v in items])])
        
    out.close()
    
    
    
   
"""

"""
def prepare_reference_genome(sid_list, output_prefix=None, bacterialdb_path=NCBI_BACTERIAL_GENEOMES_DB):
    print_status("Preparing reference genomes")
    



"""
 This routine detects and estimates the percentage of reads originated from a host genome (human)
"""
def detect_host_genome(read_fn, host="human"):
    print_status("Preparing reference genomes")
    
       

####### CAZy Specific stage #########
####### Auxillary routines #########
def parse_opt():
    print_status("Calling parsing opts")


# Print the usage of this script 
def print_usage():
    print("An integrated pipeline for processing metagenomic sequencing reads (paired-end).")
    print(" ")
    print("Usage:")
    print(" python " + sys.argv[0] + " READ_1.fq READ_2.fq [OUTPUT_PREFIX]")
    #print("  python -i BLAST-RESULT-INFILE -o FILTER-OUTFILE [-q] [-b BITSCORE-CUTOFF] [-l ALIGNMENT-LENGTH-CUTOFF] [-p PERCENTAGE-IDENTITY]")
    #print("      -i STRING  Input file generated by BLAST with -m6 option")
    print(" ")
    print(" ") 
    print("Ver 0.2b")


def create_link():
    print_status("create_link")


# Assert 
def assert_proc(fn_to_be_asserted):
    return os.path.isfile(fn_to_be_asserted) and os.stat(fn_to_be_asserted)[stat.ST_SIZE] != 0
        

# 
def import_settings():
    print_status("Importing settings not yet implemented")

# Output message to a log file
def log(log_msg):
    if LOG:
        LOG_FILE.write(log_msg + "\n")


# Print status
def print_status(msg):
    caller_name = inspect.stack()[1][3]
    msg = "[ " + caller_name + " ] " + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "  " + msg
    print(msg)
    
    log(msg)

    

# Invoke the main function
if __name__ == "__main__":
    main(sys.argv[1:])