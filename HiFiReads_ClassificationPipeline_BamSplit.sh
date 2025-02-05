## This pipeline is designed for the SLURM workload manager.
## Please modify the pipeline to suit your specific scheduling system or requirements.

 


#source $HOME/.bash_profile
HPC_accountName='mutationalscanning'
GenomeAssembly='hg38'  # Supports only Hg38 genome assembly
codeDir='/faststorage/project/mutationalscanning/Workspaces/vinod/codes/Codes_HiFiReadsClassification/src/'
dataName='marlock'
baseDir='/home/vinodsingh/mutationalscanning/Workspaces/vinod/'
DataDir="${baseDir}/${dataName}_chrBAMs_${GenomeAssembly}/"
bam=${DataDir}/${dataName}_${GenomeAssembly}.bam  # bam file should be named as "HT45_chrBAMs_hg38"
OUTDIR="${DataDir}/${dataName}_BinoClassificationResults/"
mkdir -p $OUTDIR
cd $OUTDIR



#### split bams ##
script_00='00_createDir_and_split_bam.sh'
OUTDIR_00="${OUTDIR}/00_BAMChunks/"
mkdir $OUTDIR_00
logs_00=$OUTDIR_00/logs
mkdir $logs_00
NCPUS_00=24

sbatch -J ${dataName}_splitBAM --time 10:00:00 --cpus-per-task $NCPUS_00 --mem=100GB --account $HPC_accountName \
  $codeDir/${script_00} \
  -B $bam \
  -S $dataName \
  -G $GenomeAssembly \
  -O $OUTDIR_00 \
  -N $NCPUS_00

ls -halt $OUTDIR_00
ls -halt $OUTDIR_00/${dataName}_${GenomeAssembly}_bamChunks_chr1
rm $OUTDIR_00/${dataName}_${GenomeAssembly}_bamChunks_chr*/chunk*.bam

grep -E  "===Task Completed==" $(ls  xx-${dataName}_splitBAM*.out) 

mv xx-${dataName}_splitBAM*.out $logs_00

## germFrac screen session for da1 ##

#######################################################
## Convert  splitted BAM files into chunks of RDS objects
# with binary stirng for ML and MM tag on the read ##
#######################################################
script_01='01_BAM2RDS_ByChunk.R'
OUTDIR_01="${OUTDIR}/01_BAM_RDS_ByChunk/"
mkdir $OUTDIR_01
logs_01=$OUTDIR_01/logs
mkdir $logs_01
NCPUS_01=10

parallel -v sbatch -J ${dataName}_BAM2RDS_ByChunk_chr{} --time 08:00:00 --cpus-per-task $NCPUS_01 --mem=100GB --account $HPC_accountName  \
  ${codeDir}/runR_masterCluster.sh  \
  ${codeDir}/${script_01}  \
  -I ${OUTDIR_00}/${dataName}_${GenomeAssembly}_bamChunks_chr{}  \
  -N $NCPUS_01 \
  -O ${OUTDIR_01}/BAMRDS_chr{} ::: $(seq 1 22) X Y

# varify run#
grep -E  ">> " $(ls  xx-${dataName}_BAM2RDS_ByChunk_chr*.out) | sort -V
grep -E  ">> all MM overlap with CpGs:" $(ls  xx-${dataName}_BAM2RDS_ByChunk_chr*.out) | sort -V
grep -E  "===Task Completed==" $(ls  xx-${dataName}_BAM2RDS_ByChunk_chr*.out) | sort -V
grep -E  "===Task Completed==" $(ls  xx-${dataName}_BAM2RDS_ByChunk_chr*.out) | wc -l
grep -E  "error|Error" $(ls  xx-${dataName}_BAM2RDS_ByChunk_chr*.err)
# save logs
mv xx-${dataName}_BAM2RDS_ByChunk_chr* $logs_01/


#######################################################
## chunks of RDS objects are converted to bed files ##
##  bed file has CpG methylation for each reads ##
## Approx Runtime for 1 million long reads with 10 CPUs: 8 MIN
#######################################################

script_02='02_MethyperCpG_forReads_fromHiFi_ByChunk.R'
OUTDIR_02="${OUTDIR}/02_MethyperCpG_forReads_ByChunk/"
mkdir $OUTDIR_02
logs_02=$OUTDIR_02/logs
mkdir $logs_02
NCPUS_02=12

parallel  sbatch -J ${dataName}_MethyperCpG_forReads_ByChunk_chr{} --time 15:00:00 --cpus-per-task $NCPUS_02 --mem=80GB --account $HPC_accountName  \
  ${codeDir}/runR_masterCluster.sh  \
  ${codeDir}/${script_02}  \
  -G $GenomeAssembly \
  -I ${OUTDIR_01}/'BAMRDS_'chr{}  \
  -C chr{} \
  -D $dataName \
  -O ${OUTDIR_02}/'MethyperCpG_forReads_'chr{} \
  -N $NCPUS_02 ::: $(seq 1 22) X Y

# varify run#
grep -E "error|Error" $(ls xx-${dataName}_MethyperCpG_forReads_ByChunk_chr*.*.err)
grep -E ">> completed" $(ls xx-${dataName}_MethyperCpG_forReads_ByChunk_chr*.*.out) | sort -V
grep -E ">> completed" $(ls xx-${dataName}_MethyperCpG_forReads_ByChunk_chr*.*.out) | sort -V | wc -l

# save logs
mv xx-${dataName}_MethyperCpG_forReads_ByChunk_chr* $logs_02/

#######################################################
## Reads Likelihood calculation in each celltype and Binomial classification ##
##   1. Calculate Prior component
##   2. Reads likelihood calculation
#######################################################


## calculate methylated and unmethylated component prior ##

# find methylated table #
script_03a="03a_Generate_HiFiMethyFreqTableF.sh"
MethyType='CpG'
OUTDIR_03=${OUTDIR}/03_${MethyType}_priorData/
mkdir $OUTDIR_03
logs3=$OUTDIR_03/logs
mkdir $logs3
NCPUS_03=24

sbatch -J xx-${dataName}_${MethyType}_freqTable --mem 140g --cpus-per-task=$NCPUS_03  --time=04:00:00 --account $HPC_accountName \
  ${codeDir}/${script_03a}  \
  -B $bam  \
  -D $dataName  \
  -O ${OUTDIR_03}/${dataName}_HiFiMethyFreqTable.tsv \
  -N $NCPUS_03


# estimate prior parameters from methylated table # 
script_03b="03b_GA_ParametersEstimation.R"
OUTDFILE_03="${OUTDIR_03}/Parameters.tsv"
NCPUS_03=1

sbatch -J ${dataName}_calPrior_parameters --mem 50GB --cpus-per-task=$NCPUS_03 --account $HPC_accountName --time 1:00:00 \
  ${codeDir}/runR_masterCluster.sh  \
  ${codeDir}/${script_03b} \
  -I ${OUTDIR_03}/${dataName}_HiFiMethyFreqTable.tsv \
  -D $dataName \
  -O $OUTDFILE_03 \
  -N $NCPUS_03

# save logs
mv xx-${dataName}_${MethyType}_freqTable* $logs3/
mv xx-${dataName}_calPrior_parameters* $logs3/
  

## Calculate Likelihood ##
trainData="${baseDir}/CpG_binomial_trainData/sorted_Cov.g5_Binomial_parameters_allCellTypes.rds"
#rs3_1='Reads_likelihood_InSpermCellTypes_ByChunk_intApproach.R'
script_04='04_Reads_likelihood_InSpermCellTypes_ByChunk.R'
OUTDIR_04="${OUTDIR}/04_HiFiReadsLL_InCellTypes_ByChunk/"
OUTFILE_04=${OUTDIR_04}/${dataName}_${GenomeAssembly}_allchr_LL_of_readInCellType.tsv # final output file witl LL of all chromsomes.
mkdir $OUTDIR_04
logs_04=$OUTDIR_04/logs
mkdir $logs_04
NCPUS_04=10


parallel  sbatch -J ${dataName}.ReadsLL_ByChunk_chr{} --time 08:00:00 --cpus-per-task $NCPUS_04 --mem=150GB --account $HPC_accountName  \
  ${codeDir}/runR_masterCluster.sh  \
  ${codeDir}/${script_04}  \
  -G $GenomeAssembly \
  -I ${OUTDIR_02}/'MethyperCpG_forReads_'chr{}  \
  -C chr{} \
  -D $dataName \
  -T $trainData \
  -P $OUTDFILE_03 \
  -O ${OUTDIR_04}/'Reads_likelihood_'chr{} \
  -N $NCPUS_04 ::: $(seq 1 22) X Y

# varify run#
grep -E "error|Error" $(ls xx-${dataName}.ReadsLL_ByChunk_chr*.*.err)
grep -E ">> completed" $(ls xx-${dataName}.ReadsLL_ByChunk_chr*.*.out) | sort -V
grep -E ">> completed" $(ls xx-${dataName}.ReadsLL_ByChunk_chr*.*.out) | sort -V | wc -l
cat xx-${dataName}.ReadsLL_ByChunk_chr*.out | grep -E ">>germ frac:"

# save logs
mv xx-${dataName}.ReadsLL_ByChunk_chr*.*.* $logs_04/


# merge all chromosomes LL datasets #


LL_ChrFiles=($(ls -1 ${OUTDIR_04}/${dataName}_${GenomeAssembly}_chr*_LL_of_readInCellType.tsv | sort -V)) # chromosome LL files 
if [ ${#LL_ChrFiles[@]} == 24 ]; then
  echo true
  { head -n 1 ${LL_ChrFiles[0]}; tail -n +2 -q ${LL_ChrFiles[*]}; } > $OUTFILE_04
else
  echo false
fi

head $OUTFILE_04


## add LLdiff columns to the data  ##
## Calculate germline frac at different likelihood diffrence between germ and non-germ cells.

script_05=05_calculateGermLineFraction.sh
OUTDIR_05="${OUTDIR}/05_GermFrac/"
mkdir $OUTDIR_05
logs_05=$OUTDIR_05/logs
mkdir $logs_05
NCPUS_05=10

sh ${codeDir}/${script_05} \
  -I ${OUTDIR_04}/${dataName}_${GenomeAssembly}_allchr_LL_of_readInCellType.tsv \
  -O $OUTDIR_05 \
  -N $NCPUS_05

# varify out put ##
cat ${OUTDIR_05}/germFrac.txt
(
  echo -e "LLdiff\tgermFrac\treadsRemained"
  cat ${OUTDIR_05}/germFrac.txt | awk -F'[ ,=]+' '/germFrac/ {print $4 "\t" $5 "\t" $7}'
) > ${OUTDIR_05}/germFrac2.txt
cat ${OUTDIR_05}/germFrac2.txt | column -t

##--------------xxx-------------##


