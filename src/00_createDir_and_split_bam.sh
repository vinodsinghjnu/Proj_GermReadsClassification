#!/bin/bash


usage() { echo -e "\n Usage:  $0 [options] \n\n \
Options: \n \
  [-B <BAM file (String)>]\n \
  [-S <sample Name(Sample ID)>]\n \
  [-G <GenomeAssembly name (hg38) >]\n \
  [-O <Output Dir (Dir Path)>]\n \
  [-N <Number of ncpus (Int), Defualt: 2 >]\n \
  [-h <Help>]" 1>&2; exit 1; }

ncpus=4 # default ncpus


while getopts B:S:G:O:N: flag
do
    case "${flag}" in
        B) bamPath=${OPTARG};;
        S) SampleName=${OPTARG};;
        G) GenomeAssembly=${OPTARG};;
        O) outDir=${OPTARG};;
        N)
          # Validate that the argument is a number
          if ! [[ "$OPTARG" =~ ^-?[0-9]+$ ]]; then
              echo "Error: -C option requires a valid integer argument." >&2
              usage
          fi
            ncpus=${OPTARG};;
	    *) usage;;
    esac
done

# Ensure the mandatory argument was provided
if [[ -z "$bamPath" ]]; then
    echo -e "\nError: The -B option is mandatory and was not provided." >&2
    usage
    exit 1
elif [[ -z "$SampleName" ]]; then
    echo -e "\nError: The -S option is mandatory and was not provided." >&2
    usage
    exit 1
elif [[ -z "$GenomeAssembly" ]]; then
    echo -e "\nError: The -G option is mandatory and was not provided." >&2
    usage
    exit 1
elif [[ -z "$outDir" ]]; then
    echo -e "\nError: The -O option is mandatory and was not provided." >&2
    usage
    exit 1
fi


echo ">== Parameters =="
echo "bamPath: $bamPath";
echo "SampleName: $SampleName";
echo "GenomeAssembly: $GenomeAssembly";
echo "OutDir: $outDir";
echo "ncpus avalaiable: $ncpus";
echo "====<"

nproc --all
echo sleeping20...
sleep 20



##

# Function to handle errors
error_handler() {
    echo "An error occurred. Stopping the script."
    exit 1  # Exit with a non-zero status
}

# Set trap for ERR signal
trap 'error_handler' ERR

##

source /home/vinodsingh/miniforge3/bin/activate samtools_1.16.1


#outDir=$(realpath $outDir)
#mkdir $outDir
# #bamPath=$1
# #SampleName=$2
# #GenomeAssembly=$3
# #outDir=$4

# bamPath=/faststorage/project/mutationalscanning/DerivedData/bam/HiFi/human/reference/${GenomeAssembly}
# wd=/home/vinodsingh/mutationalscanning/Workspaces/vinod/${SampleName}_chrBAMs_${GenomeAssembly}/
 ref='/home/vinodsingh/mutationalscanning/Workspaces/vinod/reference/hg38/hg38.fasta'

# source $HOME/miniforge3/etc/profile.d/conda.sh
# conda activate "samtools_1.16.1"

cd $outDir

#ln -s ${bamPath}/${SampleName}_${GenomeAssembly}.bam ${SampleName}_${GenomeAssembly}.bam
#ln -s ${bamPath}/${SampleName}_${GenomeAssembly}.bam.bai ${SampleName}_${GenomeAssembly}.bam.bai
#ls -l
#filename=${outDir}/${SampleName}_${GenomeAssembly}.bam

filename=${bamPath}


mkdir -p BAM_reports
prefixBam_with_outDirPath=${outDir}/$(basename -- ${filename%.bam})


#split BAM by chrs
echo '> Split BAM by Chr: COMMAND' 
parallel -v -j "$ncpus" "samtools view -bh $filename chr{} > ${prefixBam_with_outDirPath}_chr{}.bam" ::: $(seq 1 22) X Y
#for chrom in `seq 1 22` X Y; do echo chr$chrom; samtools view -bh $file chr${chrom} > $prefixBam_with_outDirPath${filename}_chr${chrom}.bam  done
echo  -e "\n"

## Check if BAM files are sorted ##
echo '> Check if BAM file is sorted?..:  COMMAND'
parallel -j "$ncpus" -v -k  "samtools stats ${prefixBam_with_outDirPath}_chr{}.bam | grep 'is sorted:'" ::: $(seq 1 22) X Y > BAM_reports/${SampleName}_Sort_StatusOfBAMfiles.txt
#for file in *.bam; do echo $file; samtools stats $file | grep "is sorted:" done > BAM_reports/${SampleName}_Sort_StatusOfBAMfiles.txt
echo  -e "\n"


## make index file
echo '> Making index..: COMMAND'
parallel -j "$ncpus" -v  "samtools index ${prefixBam_with_outDirPath}_chr{}.bam"  ::: $(seq 1 22) X Y
#for file in *.bam; do echo $file; samtools index $file  done
echo  -e "\n"


## add NM and MD tag to bam
echo '> Adding NM and MD tag..: COMMAND'
parallel -j "$ncpus" -v  "samtools calmd -b ${prefixBam_with_outDirPath}_chr{}.bam $ref > ${prefixBam_with_outDirPath}_MD_NMtag_chr{}.bam" ::: $(seq 1 22) X Y
rm ${prefixBam_with_outDirPath}_chr*.bam
rm ${prefixBam_with_outDirPath}_chr*.bam.bai


echo ">rename to original..: COMMAND" 
parallel -j "$ncpus" -v "mv ${prefixBam_with_outDirPath}_MD_NMtag_chr{}.bam  ${prefixBam_with_outDirPath}_chr{}.bam" ::: $(seq 1 22) X Y
echo  -e "\n"

## make index file after adding NM and MD tag.
echo '> Making index again after adding NM and MD tag..: COMMAND'
parallel -j "$ncpus" -v  "samtools index ${prefixBam_with_outDirPath}_chr{}.bam"  ::: $(seq 1 22) X Y
#for file in *.bam; do echo $file; samtools index $file  done
echo  -e "\n"

#validate bam file
echo '> Picard Validation of BAM files ..: command' 
mkdir -p BAM_reports/bam_PicardValidation_report
#parallel -j "$ncpus" -v "picard ValidateSamFile I=${SampleName}_${GenomeAssembly}_chr{}.bam -R $ref MODE=SUMMARY > BAM_reports/bam_PicardValidation_report/${SampleName}_${GenomeAssembly}_chr{}_Picard_BAMValidation.txt" ::: $(seq 1 22) X Y
parallel -j "$ncpus" -v "picard ValidateSamFile I=${prefixBam_with_outDirPath}_chr{}.bam R=$ref MODE=SUMMARY > BAM_reports/bam_PicardValidation_report/${SampleName}_${GenomeAssembly}_chr{}_Picard_BAMValidation.txt" ::: $(seq 1 22) X Y
echo  -e "\n"


## coverage of BAM files ##
echo '> Coverage of BAM files ..: COMMAND'
awk_body='{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'
parallel -j "$ncpus" -k -v "samtools depth ${prefixBam_with_outDirPath}_chr{}.bam  |  awk '$awk_body'" ::: $(seq 1 22) X Y > BAM_reports/${SampleName}_CoverageOf_BAMfiles.txt
echo  -e "\n"


## Read counts per chromosome ##
echo '> Read Counts per Chromosome of BAM files ..: COMMAND'
parallel -j "$ncpus" -k samtools view -c ${SampleName}_${GenomeAssembly}_chr{}.bam ::: $(seq 1 22) X Y > BAM_reports/${SampleName}_readCounts_perChr.txt
echo  -e "\n"


## ML amd MM tags counts ##
echo '> MMtag Counts of BAM files .. : COMMAND'

awk_body1='{for (i=12; i<=NF; ++i) { if ($i ~ "^MM:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; print td["MM"] }'
awk_body2='{sum+=$1} END{print sum}'
methy_pat="C+m"
parallel -j "$ncpus" -k -v "samtools view ${SampleName}_${GenomeAssembly}_chr{}.bam | awk '$awk_body1' | cut -d, -f1 | sort | uniq -c | grep -F '$methy_pat' | awk '$awk_body2' " ::: $(seq 1 22) X Y > BAM_reports/${SampleName}_MMtag_Counts.txt
cat BAM_reports/${SampleName}_MMtag_Counts.txt | grep -E -v '^samtools' > BAM_reports/${SampleName}_MMtag_Counts2.txt
echo  -e "\n"

## make chunks of 10k  reads for CHR bam files ##
echo '> Making chunks of BAM files ..'
echo "making dirs..: COMMAND"
parallel -j "$ncpus" -v "mkdir -p ${prefixBam_with_outDirPath}_bamChunks_chr{}" ::: $(seq 1 22) X Y
echo  -e "\n"


awk_body='BEGIN { i=0; line=0; out=outPath"/chunk"i+1".bam"}  
/^@/ {header = header $0 "\n"; next;}  
{ 
  if (line >= n && $1 != last_read) {
    close(out); 
    i++; 
    out = outPath "/chunk" i+1 ".bam"; 
    line = 0;
  }
  print (line == 0 ? header $0 : $0) | "samtools view -b -o " out;
  last_read = $1; 
  line++;
}'

echo "making chunks..: COMMAND"
parallel -j "$ncpus"  -v "samtools sort -n  -T tmpChr{}_ -O SAM "${prefixBam_with_outDirPath}_chr{}.bam" | awk -v n=10000 -v outPath="${prefixBam_with_outDirPath}_bamChunks_chr{}/" '$awk_body'" ::: $(seq 1 22) X Y
echo  -e "\n"

# 
# prefixBam_with_outDirPath='/home/vinodsingh/mutationalscanning/Workspaces/vinod/da1_chrBAMs_hg38//da1_BinoClassificationResults//BAMChunks_00//da1_hg38'
# parallel  -v "samtools sort -n  -T {} -O SAM "${prefixBam_with_outDirPath}_chr{}.bam" | awk -v n=10000 -v outPath="${prefixBam_with_outDirPath}_bamChunks_chr{}/" '$awk_body'" ::: 20
# parallel  -v "samtools sort -n  -T tmpChr{}_ -O SAM "${prefixBam_with_outDirPath}_chr{}.bam" | awk -v n=10000 -v outPath="${prefixBam_with_outDirPath}_bamChunks_chr{}/" '$awk_body'" ::: 20


#parallel -v "samtools sort my.sam > my_sorted.bam
#parallel -v "samtools index ${prefixBam_with_outDirPath}_bamChunks_chr{}/*.bam" ::: 20
#parallel -v "samtools sort {} -T {/.}_  -o {//}/sorted.{/.}.bam && samtools index {//}/sorted.{/.}.bam" :::  ${prefixBam_with_outDirPath}_bamChunks_chr20/chunk*.bam

echo "Sorting and indexing chunks..: COMMAND"
for chr in $(seq 1 22) X Y; do
    parallel -j "$ncpus" -v "samtools sort {} -T {/.}_  -o {//}/sorted.{/.}.bam && samtools index {//}/sorted.{/.}.bam" ::: \
     ${prefixBam_with_outDirPath}_bamChunks_chr${chr}/chunk*.bam
done
echo  -e "\n"

soreted_outbamsCount=$(ls  ${prefixBam_with_outDirPath}_bamChunks_chr${chr}/sorted.chunk*.bam | wc -l)
soreted_outbamsIndCount=$(ls  ${prefixBam_with_outDirPath}_bamChunks_chr${chr}/sorted.chunk*.bam.bai | wc -l)
outbamsCount=$(ls  ${prefixBam_with_outDirPath}_bamChunks_chr${chr}/chunk*.bam | wc -l)

if [[ "$soreted_outbamsCount" == "$outbamsCount" ]] & [[ "$soreted_outbamsIndCount" == "$outbamsCount" ]]; then
    echo "deleting unncessary files"
    rm ${prefixBam_with_outDirPath}_bamChunks_chr${chr}/chunk*.bam
else
    echo "Some chunks are not sorted or indexed properly "
    exit 1  # Exit with a non-zero status
fi

echo "===Task Completed=="

###

#working
# BEGIN { i=0; line=0; out=outPath"/chunk"i".bam"} /^@/ {header = header $0 "\n"; next;} { if (line >= n && $1 != last_read) { close(out); i++; out = outPath "/chunk" i ".bam"; line = 0; } print (line == 0 ? header $0 : $0) | "samtools view -b -o " out; last_read = $1; line++; }



