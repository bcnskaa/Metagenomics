#!/bin/bash

infile=$1
outfile=$2

# Command to extract GI from the fasta header
cat $infile | awk -F "\t" '{print $2}' | awk -F "|" '{print $2}' > $outfile


