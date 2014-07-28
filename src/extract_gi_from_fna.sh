#!/bin/bash

infile=$1

# Command to extract GI from the fasta header
# cat $infile | awk -F "\t" '{print $2}' | awk -F "|" '{print $2}' > $outfile



GI_LIST='gi_list.list'

if [ -e "$GI_LIST" ];then
	rm "$GI_LIST"
fi 

while read line;
do
	if [[ "$line" =~ ^\> ]];then
		IFS='|' read -ra ARRAY <<< "$line"
		# Print header
		echo ">${ARRAY[1]}"

		IFS='[' read -ra DESC <<< "${ARRAY[4]}"
		species=${DESC[1]/]/}
		
		# Export GI info
		echo -e "${ARRAY[1]}\t${ARRAY[3]}\t${DESC[0]}\t$species" >> "$GI_LIST" 
	else
		# Print sequence
		echo $line
	fi
done < "$infile"



