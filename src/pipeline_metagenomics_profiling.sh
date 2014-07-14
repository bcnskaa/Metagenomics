#!/bin/bash

#############################
# Environmental variables
#############################
SEQTK_HOME=$HOME/tools/seqtk

# paired-end sequence id
SEQ_FQ1=$1
SEQ_FQ2=$2


# Filtering fastq files using seqtk
# seqtk will trim low-quality bases from both ends using the Phred algorithm
$SEQTK/seqtk trimfq $SEQ_FQ1 > ${SEQ_FQ1/.fq/trimmed.fq}
$SEQTK/seqtk trimfq $SEQ_FQ2 > ${SEQ_FQ2/.fq/trimmed.fq}




