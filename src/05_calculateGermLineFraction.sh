#!/bin/bash

## Script adds max liklihood information of case and control classes in the input file obtained from 
## "Reads_likelihood_InSpermCellTypes_ByChunk.R" script output  
## Calculate germline fraction at different LL differences  ##


usage() { echo -e "\n Usage: $0 [options] \n\n \
Options: \n \
  [-I <Input file (String), output of 04_Reads_likelihood_InSpermCellTypes_ByChunk.R >]\n \
  [-O <Output Directory ]\n \
  [-N <Number of CPUs (Int)>]\n \
  [-h <Help>]" 1>&2; exit 1; }

cpus=3 # default CPU counts

while getopts ":I:O:N:" flag; do
    case "${flag}" in
 	    I) INPUT_FILE=${OPTARG};;
      O) OUTDIR=${OPTARG};;
      N)
         # Validate that the argument is a number
         if ! [[ "$OPTARG" =~ ^-?[0-9]+$ ]]; then
             echo "Error: -N option requires a valid integer argument." >&2
             usage
         fi
         ncpus=${OPTARG};;
	    *) usage;;
    esac
done


### Default output file ##

if [[ -z "$OUTDIR" ]]; then
    OUTDIR=$(pwd)
    echo "output directory is $OUTDIR"
fi

# Ensure the mandatory argument was provided
if [[ -z "$INPUT_FILE" ]]; then
    echo -e "\nError: The -I option is mandatory and was not provided." >&2
    usage
    exit 1
fi


##

#INPUT_FILE=$LL_allChrFile
#OUTDIR=$(pwd)
OUTPUT_FILE=${OUTDIR}/$(basename "${INPUT_FILE%.*}"_withLLDiff.tsv)
tempDir=${OUTDIR}/TEMP_$(openssl rand -hex 12)
mkdir $tempDir

## Awk function to calulate max likelihood value and its index in corresponding case and control samples
max_Ind_n_Val_awk_body='NR >1 {
    max_val = $1   # Assume the first column is the maximum initially
    max_idx = 1    # Index of the maximum value
    for (i = 2; i <= NF; i++) {   # Loop through each field from the second onwards
        if ($i > max_val) {
            max_val = $i
            max_idx = i
        }
    }
    #print "Row", NR, "Max value:", max_val, "at index:", max_idx
    printf "%s,%s\n", max_val, max_idx
}'


# Calculate maxLL_Germ
cat "$INPUT_FILE" | cut -f 4-8  | awk "$max_Ind_n_Val_awk_body" > ${tempDir}/tmp_maxLL_Germ.txt
cat "$INPUT_FILE" | cut  -f 9-10  | awk "$max_Ind_n_Val_awk_body" > ${tempDir}/tmp_maxLL_noGerm.txt

# Join the columns and save as a TSV file
paste -d $'\t' ${tempDir}/tmp_maxLL_Germ.txt ${tempDir}/tmp_maxLL_noGerm.txt | tr ',' '\t' > ${tempDir}/tmp.tsv # replcae coma by tab seperator
#cat ${tempDir}/tmp.tsv | awk '{print $1-$3}' > ${tempDir}/tmp_LL_maxGerm_minus_maxNonGerm.tsv
cat ${tempDir}/tmp.tsv | awk -F'\t' '{diff=$1-$3; if (diff < 0) {  diff_abs=-diff } else { diff_abs=diff }; print diff "\t" diff_abs }' >  ${tempDir}/tmp_LL_maxGerm_minus_maxNonGerm.tsv
#paste -d $'\t' ${tempDir}/tmp.tsv ${tempDir}/tmp_LL_maxGerm_minus_maxNonGerm.tsv > ${tempDir}/tmp_2.tsv
paste -d $'\t' ${tempDir}/tmp_maxLL_Germ.txt ${tempDir}/tmp_maxLL_noGerm.txt ${tempDir}/tmp_LL_maxGerm_minus_maxNonGerm.tsv > ${tempDir}/tmp_2.tsv

head ${tempDir}/tmp_2.tsv

col_names="maxLL_Germ\tmaxLL_noGerm\tLL_maxGerm_minus_maxNonGerm\tabs_LL_maxGerm_minus_maxNonGerm"
echo -e "$col_names" > ${tempDir}/tmp_ll_data.tsv
cat ${tempDir}/tmp_2.tsv >> ${tempDir}/tmp_ll_data.tsv
paste -d $'\t' $INPUT_FILE ${tempDir}/tmp_ll_data.tsv > $OUTPUT_FILE

## add germ or nogerm column #
awk -F'\t' 'BEGIN {OFS="\t"} NR==1 {print $0, "germ"; next} {print $0, ($(NF-1) > 0 ? "TRUE" : "FALSE")}' "$OUTPUT_FILE" > "${OUTPUT_FILE}2"
mv "${OUTPUT_FILE}2" "$OUTPUT_FILE"
head $OUTPUT_FILE | column -t
#head $OUTPUT_FILE | cat -t



# Clean up temporary files
rm -r ${tempDir}
#rm ${tempDir}/tmp_maxLL_Germ.txt ${tempDir}/tmp_maxLL_noGerm.txt ${tempDir}/tmp_LL_maxGerm_minus_maxNonGerm.tsv ${tempDir}/tmp.tsv


## Calculate germline frac at different likelihood diffrence between germ and non-germ cells.

germ_frac() {
  local ll_th="$1"
  local file="$2"

  # Check if OUTPUT_FILE exists and is readable
  if [[ ! -r "$file" ]]; then
    echo "Error: OUTPUT_FILE does not exist or cannot be read."
    return 1
  fi

  all_reads=$(awk -v ll="$ll_th" 'NR >1 && $(NF-1) > ll {count++} END {print count}' "$file")
  totalreadsInFile=$(cat "$file" | wc -l)
  
  # Handle case where all_reads could be zero
  if [[ "$all_reads" == "0" ]]; then
    echo "No reads above threshold."
    return
  fi

  #all_reads=$(cat $file | awk -v ll="$ll_th" 'NR >1 {if($NF > ll) print $NF}' | wc -l)
  germ_reads=$(cat $file | awk -v ll="$ll_th" 'NR >1 {if($(NF-1) > ll && $NF=="TRUE") print $(NF-1)}' | wc -l)
  germ_frac=$(echo "scale=3;($germ_reads/$all_reads)*100" | bc)
  reads_remained=$(echo "scale=3;($all_reads/$totalreadsInFile)*100" | bc)
  echo "germFrac at ll> $ll_th = $germ_frac , readsRemained= $reads_remained"
}

export -f  germ_frac

parallel -k  germ_frac {} "$OUTPUT_FILE"  ::: $(seq 0 25) > ${OUTDIR}/germFrac.txt

