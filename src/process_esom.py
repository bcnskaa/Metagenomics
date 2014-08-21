#!/share/apps/Python-2.7.4/bin/python

# 2014 SKWoolf bcnskaa AT gmail DOT com

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



"""
 R code:
 
 dom <- read.table("tmp.dom", header=F, sep="\t", stringsAsFactors=F);
 dom <- sapply(seq(1,nrow(dom)), function(i) strsplit(dom$V1[i], "\\s+"));
 dom_n <- length(dom);
 
 dom_df <- data.frame(contig_name=character(dom_n), bin_name=);
 for(i in 1 : dom_n)
 {
     dom[[i]] <- 
 } 
 
"""