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
DECONSEQ_HOME = TOOLS_HOME + "/deconseq"
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
NCBI16S_GI_IDX = DB_HOME + "/Markers/ncbi16s/16SMicrobial.idx"

EMIRGE_16S_DB = EMIRGE_HOME + "/db/SSURef_111_candidate_db.fasta"

CAZY_DB = DB_HOME + "/Markers/CAZy/CAZy_id.lst.retrieved.faa"
DBCAN_HOME = DB_HOME + "/Markers/dbCAN"
DBCAN_HMM = DBCAN_HOME + "/dbCAN-fam-HMMs.txt.v3.txt"  # CAZy HMM profile 
STRING_DB_HOME = DB_HOME + "Markers/STRING"


# Folder contains all prokaryotic genomes available in NCBI
NCBI_BACTERIAL_GENEOMES_DB = DB_HOME + "/BacterialDB/all_fna"


LOG = False
LOG_FN= "pipeline.log"
LOG_FILE = open(LOG_FN, "w")

# Output directories
MARKER_OUTDIR = "Markers"
RESULT_OUTDIR = "Reports"
BINNING_OUTDIR = "Binning"
TMP_OUTDIR = "temp"

HMMER_OUTDIR = MARKER_OUTDIR + "/HMMER"

VERBOSE_ONLY = False



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
    


"""
 Mapping identified HMM domains to bin_groups
 Keys from hmm_orf_dict should be ORFs names generated by Prodigal. They are contig prefix + "_" + number.
 The contig prefixes are expected to be found in binning groups by MaxBin
 Get Current Dir: os.path.dirname(os.path.abspath("."))
"""
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
        print_status("Unable to read *.summary from"+maxbin_dir)
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
                bin_id = items[0].replace("." + FASTA_DNA_EXT, "")
               
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
    

# Return the length of overlapping between two regions
def getOverlap(x, y):
    #print "x[1]=", x[1], ", y[1]=", y[1], ", x[0]=", x[0], ", y[0]=", y[0]
    return max(0, min(x[1], y[1]) - max(x[0], y[0]))



# Print status
def print_status(msg):
    caller_name = inspect.stack()[1][3]
    msg = "[ " + caller_name + " ] " + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "  " + msg
    print(msg)
    
    
"""

 import process_HMM
 maxbin_path="../../../../MaxBin"
 out_fn = "test.tbl"
 hmm_orf_dict = process_HMM.postprocess_HMMER_search(hmm_dir=".")
 bin_groups = process_HMM.map_hmm2maxbin(hmm_orf_dict, maxbin_path)
 process_HMM.generate_map_table(bin_groups, out_fn)
 
"""
    
    
"""
id=P4_contig-Prodigal
~/tools/hmmer/bin/hmmsearch -o $id-dbCAN.out -A $id-dbCAN.aln --tblout $id-dbCAN.tbl --domtbl $id-dbCAN.dom.tbl --pfamtblout $id-dbCAN.pfam.tbl ~/db/Markers/dbCAN/dbCAN-fam-HMMs.txt.v3.txt ../../../../Prodigal/$id.all.faa 


swh_cell <- read.table("SWH-cell35_S1-dbCAN.HMM.tbl", sep="\t", header=T, stringsAsFactors=F)
swh_cell <- swh_cell[,-ncol(swh_cell)]


swh_xyl <- read.table("SWH-xyl35_S3.all-dbCAN.HMM.tbl", sep="\t", header=T, stringsAsFactors=F)
swh_xyl <- swh_xyl[,-ncol(swh_xyl)]


gz_xyl <- read.table("GZ-xyl_S2-dbCAN.HMM.out.tbl", sep="\t", header=T, stringsAsFactors=F)
gz_xyl <- gz_xyl[,-ncol(gz_xyl)]


gz_cell <- read.table("GZ-cell_S1-dbCAN.HMM.tbl", sep="\t", header=T, stringsAsFactors=F)
gz_cell <- gz_cell[,-ncol(gz_cell)]

p1 <- read.table("P1_contig-Prodigal-dbCAN.HMM.tbl", sep="\t", header=T, stringsAsFactors=F)
p1 <- p1[,-ncol(p1)]

p2 <- read.table("P2_contig-Prodigal-dbCAN.HMM.tbl", sep="\t", header=T, stringsAsFactors=F)
p2 <- p2[,-ncol(p2)]

p3 <- read.table("P3_contig-Prodigal-dbCAN.HMM.tbl", sep="\t", header=T, stringsAsFactors=F)
p3 <- p3[,-ncol(p3)]

p4 <- read.table("P4_contig-Prodigal-dbCAN.HMM.tbl", sep="\t", header=T, stringsAsFactors=F)
p4 <- p4[,-ncol(p4)]


samples <- list()
samples[["SWH-Cell-Y2"]] <- swh_cell;
samples[["SWH-Xyl-Y2"]] <- swh_xyl;
samples[["GZ-Cell-Y2"]] <- gz_cell;
samples[["GZ-Xyl-Y2"]] <- gz_xyl;
samples[["SWH-Xyl-Y1"]] <- p1;
samples[["GZ-Xyl-Y1"]] <- p2;
samples[["SWH-Cell-Y1"]] <- p3;
samples[["GZ-Cell-Y1"]] <- p4;


hmm_dom_list <- c(swh_cell$HMM_DOM, swh_xyl$HMM_DOM, gz_xyl$HMM_DOM, gz_cell$HMM_DOM, p1$HMM_DOM, p2$HMM_DOM, p3$HMM_DOM, p4$HMM_DOM)
hmm_dom_list <- as.character(unique(hmm_dom_list))

#
hmm_dom_list <- hmm_dom_list[grep("GH|CBM|GT|CE", hmm_dom_list)]
hmm_dom_list <- hmm_dom_list[grep("GH", hmm_dom_list)]


pdf_fn <- "test.pdf"
bin_n <- 5;
gene_n_threshold <- 1;
max_gene_n <- 8;
min_value_per_row_threshold <- 4;

hmm_tbl <- matrix(0, length(hmm_dom_list), length(samples));
rownames(hmm_tbl) <- hmm_dom_list;
colnames(hmm_tbl) <- names(samples);
for(hmm in hmm_dom_list)
{
    for(sample in names(samples))
    {
        idx = which(samples[[sample]][,1] == hmm);
        if(length(idx) == 1)
        {
            hmm_tbl[hmm, sample] <- sum(samples[[sample]][idx, 2:(bin_n + 1)]);
        }
    }
}

sum_per_row <- sapply(1:nrow(hmm_tbl), function(i) sum(hmm_tbl[i,2:ncol(hmm_tbl)]) )
min_per_row <- sapply(1:nrow(hmm_tbl), function(i) min(hmm_tbl[i,2:ncol(hmm_tbl)]) )


plot_data <- melt(hmm_tbl[-which(sum_per_row < 6 | min_per_row > 6), ]);
plot_data$Var2 <- factor(plot_data$Var2, level=c("SWH-Cell-Y2", "SWH-Cell-Y1", "GZ-Cell-Y2", "GZ-Cell-Y1", "SWH-Xyl-Y2", "SWH-Xyl-Y1", "GZ-Xyl-Y2", "GZ-Xyl-Y1"))

#plot_data$Var1 <- as.character(plot_data$Var1);
#sorted_labels <- c(read.table("Sorted.txt", header=F, sep="\n", quote="", strip.white=TRUE, stringsAsFactors=F))[[1]]
sorted_labels <- data.frame(label=rownames(hmm_tbl), diff=sapply(1:nrow(hmm_tbl), function(i) sum(hmm_tbl[i,c("GZ-Cell-Y1", "GZ-Cell-Y2", "SWH-Cell-Y1", "SWH-Cell-Y2")]) - sum(hmm_tbl[i,c("GZ-Xyl-Y1", "GZ-Xyl-Y2", "SWH-Xyl-Y1", "SWH-Xyl-Y2")])))
sorted_labels <- sorted_labels[order(sorted_labels$diff), ]
plot_data$Var1 <- factor(plot_data$Var1, level=as.character(sorted_labels$label))

plot_data <- cbind(plot_data, label=plot_data$value);
plot_data$value[which(plot_data$value > max_gene_n)] <- max_gene_n;

color_scale <- colorRampPalette(c("steelblue", "yellow", "red"))(n = 20)
library(reshape2)
library(ggplot2)

legend_title = "CAZy Domain#";
pdf(pdf_fn, width=20, height=4);
ggplot(plot_data, aes(x=Var1, y=Var2, fill=value)) + geom_tile(aes(height=0.97, width=0.97)) +
  scale_fill_gradientn(colours=color_scale, name=legend_title) +
  labs(x="CAZy Class", y="Sample") +
  geom_text(aes(label=label, size=1)) + 
  theme(panel.background=element_blank(), axis.text.x=element_text(angle=315, hjust=0, vjust=1));
dev.off()



sapply(1:nrow(hmm_tbl), function(i) sum(hmm_tbl[i,c("GZ-Cell-Y1", "GZ-Cell-Y2", "SWH-Cell-Y1", "SWH-Cell-Y2")]) - sum(hmm_tbl[i,c("GZ-Xyl-Y1", "GZ-Xyl-Y2", "SWH-Xyl-Y1", "SWH-Xyl-Y2")]))


"""
