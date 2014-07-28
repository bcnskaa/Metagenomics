#!/bin/bash


# If subject ids have patterns like these:
# gi|123467|NC_12345.7|
# This script will extract the gi id and replace the original subject id

infile=$1
outfile=$2



echo "Reading from $infile"


#IFS='\t' read -a ARR << $infile
#while IFS=' ' read -ra ARR; do
#  echo "${ARR[2]}"
#for i in "${ARR[@]}"; do
#    echo "$i"
# done

while read line;
do
  IFS=$'\t'
  read -a ITEMS <<< "$line"
  
  IFS=$'|'
  read -a SID <<< "${ITEMS[1]}"
  
  # Replace the original sid
  ITEMS[1]=${SID[1]}

  IFS=$'\t'
  cat <<< "${ITEMS[*]}" >> "$outfile"

done < "$infile"


function join { local IFS="$1"; shift; echo "$*"; }
