from __future__ import division

import glob
import os
import numpy
import sys
from Bio import SeqIO


sys.path.append(os.path.abspath("/home/siukinng/tools/scripts"))


import mg_pipeline
import pick_seq


"""
#https://www.biostars.org/p/45654/

Count the number of unmapped reads
samtools view -f4 whole.bam | cut -f10 | sort | uniq -c | sort -nr > unmapped_unique.count


http://seqanswers.com/forums/showthread.php?t=23389
Count the number of reads mapped mutliple times
grep 'XS:' your_alignment_file.sam | wc -l

Remove lines with "XS" tag from a sam file
sed '/XS:/d' your_alignment_file.sam > your_alignment_file_1alignmentonly.sam



"""



