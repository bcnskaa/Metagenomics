#!/bin/bash

infile=$1

# Command to extract GI from the fasta header
# cat $infile | awk -F "\t" '{print $2}' | awk -F "|" '{print $2}' > $outfile

# Extract GI# from original fasta headers, and use them as new fasta headers
# cat all_ref.fa | cut -d'|' -f1,2 | sed 's/gi|//g' > all_ref.gi.fa
# Convert Fasta GI header into tabular format
# cat all_ref.fa | grep -e "^>"  | sed 's/\,\ complete genome//g' | sed 's/, complete sequence//g' | sed 's/>gi|//' | sed 's/ref|//'| sed 's/|/\t/g


# Generate a GI info list of all_faa.faa
#cat all_faa.faa | grep ">" | sed 's/ \[/|/g' | sed 's/\]//g' | sed 's/| /\|/g' | cut -d'|' -f2,4,5,6 | sed 's/|/\t/g' > all_faa.lst
# Convert all_faa.faa headers into gi headers
#cat all_faa.faa | cut -d'|' -f1,2 | sed 's/gi|//g' > all_faa.gi.faa &



# Convert the fasta header to GI number only
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



