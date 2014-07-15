#!/bin/bash

########################################################################################
#
# This is a pipeline script for reconstruction of biological pathways from NGS-based
# metagenomic data. Data from Illumina and SOLiD can be handled in this pipeline.
#
# Flow of the Pipeline
# 	Stage 0 - Preprocessing
# 		1. Filter raw fastq datasets
# 
#	Stage 1 - Taxonomic Mappings and Sequence Assembly
# 		1. Construction of a preliminary taxonomic profile using 16s rRNA sequences
# 		2. Inference of ORFs/genes
# 		3. Mapping and assemblies of short reads
#			a. de novo
#			b. guided by 16s taxonomic profiles  
#
#	Stage 2 - Functional Analysis
#		1. Taxonomic analysis
#		2. Functional analysis of proteins
#	
#	Stage 3 - Pathway Analysis
#		1. Recontruction of genomic model
#		2. Mapping of biological pathways
#		3. Stoichiometric modelling
#
#	Stage 4 - Community Analysis
#		
########################################################################################



# User input argument
# Sample id
SAMPLE_ID=$1


# paired-end sequence id
SEQ_FQ1=$1
SEQ_FQ2=$2
#SEQ_FQ1="$SAMPLE_ID\_1.fq"
#SEQ_FQ2="$SAMPLE_ID\_2.fq"


#############################
# Environmental variables
#############################
TOOLS_HOME=$HOME/sk/tools
SEQTK_HOME=$TOOLS_HOME/seqtk
METAPHLAN_HOME=$TOOLS_HOME/Metaphlan


##########################################################
# Stage 0 - Preprocessing
##########################################################

inputfile_1=$SEQ_FQ1
inputfile_2=$SEQ_FQ2
outputfile_1=${inputfile_1/.fq/trimmed.fq}
outputfile_2=${inputfile_2/.fq/trimmed.fq}

# Filtering fastq files using seqtk
# seqtk will trim low-quality bases from both ends using the Phred algorithm
$SEQTK_HOME/seqtk trimfq $inputfile_1 > $outputfile_1
$SEQTK_HOME/seqtk trimfq $inputfile_2 > $outputfile_2


# Generate FASTQC report


##########################################################
# Stage 1 - Taxonomic Mappings and Sequence Assembly
##########################################################


# Using MetaPhlAn for preliminarily profiling the composition of microbial community
inputfile_1=$outputfile_1
outputfile_1="$inputfile_1\.BM.txt"
inputfile_2=$outputfile_2
outputfile_2="$inputfile_2\.BM.txt"
output_folder="profile_samples"

# Read_1
eval "cat $input_filename | $METAPHLAN_HOME/metaphlan.py --bowtie2db $METAPHLAN_HOME/db/bowtie2db/mpa --bt2_ps very-sensitive --input_type multifastq > ./$output_folder/$output_filename"

# Read_2
input_filename="$SEQ_FQ1"
output_filename=$SEQ_FQ1\.BM.txt"
eval "cat $input_filename | $METAPHLAN_HOME/metaphlan.py --bowtie2db $METAPHLAN_HOME/db/bowtie2db/mpa --bt2_ps very-sensitive --input_type multifastq > ./$output_folder/$SAMPLE_ID_1.BM.txt"


# Generate 


##########################################################
# Stage 2 - Functional Analysis
##########################################################







