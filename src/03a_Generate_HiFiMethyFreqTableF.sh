#!/bin/bash

## script produces frequency of methylation values for each read along rows,
# followed by splitting file germline and non-germline reads ##


usage() { echo -e "\n Usage: $0 [options] \n\n \
Options: \n \
  [-B <BAM file (String)>]\n \
  [-O <Output File ]\n \
  [-D <Data Name(Sample ID)>]\n \
  [-N <Number of CPUs (Int)>]\n \
  [-h <Help>]" 1>&2; exit 1; }

while getopts ":B:O:D:N:" flag; do
    case "${flag}" in
 	      B) bam=${OPTARG};;
        O) outFile=${OPTARG};;
        D) dataName=${OPTARG};;
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



###

awk_body1='{for (i=12; i<=NF; ++i) { if ($i ~ "^ML:B:"){ td[substr($i,1,2)] = substr($i,7,length($i)); } }; print td["ML"] }'
awk_body2='{print $2, $1}'
#awk_body3='{sum+=$1} {print $1/sum, $2/255}'


s=(1 $(seq 2 2 48))
s2=( "${s[@]/#/\$}" )
printf '%s ' "${s2[@]}"
s3=$(printf '%s,' "${s2[@]}")
s4=${s3%?}
awk_body4="{print $s4}"

echo -e "\n"
##

source /home/vinodsingh/miniforge3/bin/activate "samtools_1.16.1"



outDir2=$(dirname $outFile)
echo ">> Output Dir: $outDir2"

parallel -j $ncpus -v "samtools view  $bam  chr{} |  awk '$awk_body1'  | awk '{print}' ORS='' |  tr ',' '\n' | sed 1,1d |  sort | uniq -c | sort -nk2 | awk '$awk_body2' > ${outDir2}/tmp_${dataName}_HiFiMethyFreqTable_chr{}.txt " ::: $(seq 1 22) X Y
paste -d' ' $(ls ${outDir2}/tmp_${dataName}_HiFiMethyFreqTable_chr*.txt) >  ${outDir2}/tmp_${dataName}_HiFiMethyFreqTable.txt
awk -f <(echo "$awk_body4") ${outDir2}/tmp_${dataName}_HiFiMethyFreqTable.txt > ${outDir2}/tmp_${dataName}_HiFiMethyFreqTable2.txt
cat ${outDir2}/tmp_${dataName}_HiFiMethyFreqTable2.txt | awk '{for(i=2;i<=NF;i++) t+=$i; print $1, t; t=0}' > "$outFile"
echo -e "MethyLevel\tFreq" | cat - "$outFile" | sed 's/ \+/\t/g'  > temp && mv temp "$outFile"  # add colnames
rm ${outDir2}/tmp_${dataName}_*.txt


#GeneralisedScript="${codeDir}/Generate_HiFiMethyFreqTableF.sh"
#sbatch -J chunk_${dataName}_${MethyType} --mem 140g --cpus-per-task=25 --account mutationalscanning --time=04:00:00 $GeneralisedScript -B $bamFile -D $dataName -O $MethyOutDir -C 25
