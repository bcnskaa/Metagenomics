#!/bin/bash

########################################################################################
#
# This is a pipeline script for reconstruction of biological pathways from NGS-based
# metagenomic datasets. Data from Illumina and SOLiD can be handled in this pipeline.
#
# Flow of the Pipeline
# 	Stage 0 - Preprocessing
# 		1. Filter raw fastq datasets
#		2. Quality checking
#		3. Low complexity filter (LCF)
# 
#	Stage 1 - Taxonomic Mappings and Sequence Assembly
# 		1. Construction of a preliminary taxonomic profile using 16s amplicon dataset
# 		2. Prediction of ORFs/genes directly from reads
# 		3. Mapping and assemblies of short reads
#			a. de novo (idba_ud, metavelvet)
#			b. guided by reference library based on 16s taxonomic profiles (metavelvet-SL)
#		4. Binning (contigs with length >800bp)
#			a. Sequence similarity of contigs to reference genomes (blastn)
#			b. Yield contigs with best hits to reference genomes
#			c. Classify contigs to 
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

# Working directory
WD=

# paired-end sequence id
#SEQ_FQ1=$1
#SEQ_FQ2=$2
SEQ_FQ1="$SAMPLE_ID\_1.fq"
SEQ_FQ2="$SAMPLE_ID\_2.fq"



#############################
# Environmental variables
#############################
TOOLS_HOME=$HOME/tools
SEQTK_HOME=$TOOLS_HOME/seqtk
METAPHLAN_HOME=$TOOLS_HOME/Metaphlan
VELVET_HOME=$TOOLS_HOME/velvet
FASTQC_HOME=$TOOLS_HOME/FastQC
FASTQ_MCF_HOME=$TOOLS_HOME/fastq-mcf

ADAPTOR_SEQS="/home/siukinng/db/Sequencing/adaptors.fa"

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
$FASTQC_HOME/fastqc $inputfile_1 --outdir=fastqc_reports
$FASTQC_HOME/fastqc $inputfile_2 --outdir=fastqc_reports
$FASTQC_HOME/fastqc $outputfile_1 --outdir=fastqc_reports
$FASTQC_HOME/fastqc $outputfile_2 --outdir=fastqc_reports


# Quality filter using fastq-mcf
inputfile_1=$outputfile_1
inputfile_2=$outputfile_2
outputfile_1=${inputfile_1/.fq/cleaned.fq}
outputfile_2=${inputfile_2/.fq/cleaned.fq}
$FASTQ_MCF_HOME/fastq-mcf -o $outputfile_1 -o $outputfile_2 -l 16 -q 15 -w 4 -x 10 $ADAPTOR_SEQS $inputfile_1 $inputfile_2

# Generate FASTQC report
$FASTQC_HOME/fastqc $outputfile_1 --outdir=fastqc_reports
$FASTQC_HOME/fastqc $outputfile_2 --outdir=fastqc_reports


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


# Using MetaVelvet-SL to reconstruct scaffolds
$VELVET_HOME/velveth $SAMPLE_ID\_MetaVelvet_SL 100 -shortPaired -fastq [read file] 



# ************************************************
# For MetaVelvet-SL
# ************************************************
set TOOLS_HOME = $HOME"/sk/tools"
set VELVET_HOME = $TOOLS_HOME"/velvet"
set METAVELVET_HOME = $TOOLS_HOME"/MetaVelvet-SL"
set WD = "/disk/rdisk08/jiapchen/sk/pipelines/assembling/lab"
set SVMLIB_HOME=$TOOLS_HOME/libsvm
set MODEL="SoilFeatureAll"

#Must 'cd' to the working directory of the script as shown below to run your job
cd $WD


# Go through all the lab samples, P1, P2, P3 and P4
foreach f (`ls *_1.trimmed.fq`)
#foreach f (`ls P4_1.trimmed.fq`)
  set id = "`echo $f | sed -e "s/_1.trimmed.fq//g"`"
  set ifn1 = $id"_1.trimmed.fq"
  set ifn2 = $id"_2.trimmed.fq"
  set outdir = $WD"/metavelvet/"$id

  set kmer = 31
  set insert_size = 160
  
  
  # Create de Bruijn graph using velvet 
  set cmd = "$VELVET_HOME/velveth $outdir $kmer -shortPaired -fastq $ifn1 $ifn2 "
  echo "************************************************"
  echo $cmd
  eval $cmd
  
  
  echo "************************************************"
  set cmd = "$VELVET_HOME/velvetg $outdir -ins_length $insert_size -read_trkg yes -exp_cov auto"
  echo $cmd
  eval $cmd
  
  
  #  Extract candidates of chimera nodes
  echo "************************************************"
  set cmd = "$METAVELVET_HOME/meta-velvete $outdir -ins_length $insert_size"
  echo $cmd
  eval $cmd


  # Extract feature 
  echo "************************************************"
  set cmd = "perl $METAVELVET_HOME/LearningModelFeatures/FeatureExtractPredict.perl $outdir/meta-velvetg.subgraph__TitleChimeraNodeCandidates $outdir/Features $outdir/Features3Class $outdir/ChimeraNodeName"
  echo $cmd
  eval $cmd
  
  # Doing classification using pretrained model 
  echo "************************************************"
  set cmd = "$SVMLIB_HOME/svm-scale -r $SVMLIB_HOME/tools/LibraryPretrainedClassification/$MODEL.range $outdir/Features3Class > $outdir/Features3Class.scale"
  echo $cmd
  eval $cmd
  
  # Doing classification using pretrained model 
  echo "************************************************"
  set cmd = "$SVMLIB_HOME/svm-predict $outdir/Features3Class.scale $SVMLIB_HOME/tools/LibraryPretrainedClassification/$MODEL.model $outdir/P3/ClassificationPredict"
  echo $cmd
  eval $cmd
  
  # Doing assembly using pretrained model 
  echo "************************************************"
  set cmd = "$METAVELVET_HOME/meta-velvetg $outdir/meta-velvetg.subgraph__TitleChimeraNodeCandidates $outdir/Features $outdir/Features3Class $outdir/ChimeraNodeName"
  echo $cmd
  eval $cmd
  
end

# perl /home/jiapchen/sk/tools/MetaVelvet-SL/LearningModelFeatures/BLAST_map/eval.pl -i /disk/rdisk08/jiapchen/sk/pipelines/assignments/lab/models/P2_model/meta-velvetg.subgraph__ChimeraNodeCandidates -n P2_model -d P2_model_all_refgenomes.fna -p P2_model -L 120


##########################################################
# Stage 2 - Functional Analysis
##########################################################


# Using Prodigal to predict ORFs from the contig data






# Binning using ESOM
set WORKING_DIR=""
set FASTA_DIR=""
set FASTA_EXT="fa"
set TETRA_ESOM_HOME=""


# Filter the sequence shorter than 2500bp


# Calculate the ESOM
perl $TETRA_ESOM_HOME"/"esomWrapper.pl -path $FASTA_DIR -ext $FASTA_EXT -scripts $TETRA_ESOM_HOME





# Extract marker gene ids
# cat P1_contig-MaxBin.marker | grep -e "Total" | cut -f19- | sed -e "s/\t/\n/g"

# Extract marker gene description from the corresponding HMM file
# for id in `cat marker_ids.lst`;do desc=`grep "DESC" /home/siukinng/db/Markers/TIGRFAM/HMM/$id.HMM`;echo -e "$id\t${desc/DESC /}"; done 







