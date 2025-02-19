# HiFi Reads Classification Pipeline

This method was trained on published methylation data from sorted cell types of human spermatogenesis and somatic proxies (blood and neurons). We can reliably classify the HiFi reads to somatic versus germline.

This pipeline is designed for the SLURM workload manager.
Please modify it to suit your specific scheduling system or requirements.


## Step 0: Build conda environments to run the script
 
- **GermReadsClassification_Renv:** conda env create -f test_GermReadsClassification_Renv.yml
- **GermReadsClassification_NGSenv:** conda env create -f test_GermReadsClassification_NGSenv.yml

- **Packages used to build conda environment**
  - **Linux Packages:** GNU parallel 20240722, samtools 1.16.1, picard 3.0.0, R 4.2.3 (GermReadsClassification_NGSenv)
  - **R-packages:**  optparse 1.7.5,  parallel 4.2.3, plyranges 1.18.0 , EnvStats 3.0.0, tibble 3.2.1, gtools 3.9.4, patchwork 1.3.0, GA 3.2.3, ggplot2 3.5.1 , BSgenome.Hsapiens.UCSC.hg38 1.4.5, BSgenome 1.66.3, Biostrings 2.66.0, Rsamtools, GenomicAlignments 1.34.1, rtracklayer 1.58.0, (Detailed in "Rsession_info.txt") (GermReadsClassification_Renv)

## Step 1: Setup and Variables

- Sets the HPC account, genome assembly, and paths for input and output directories.
- Ensures the output directory exists.

```bash
# Define HPC account name and genome assembly version
HPC_accountName='mutationalscanning'
GenomeAssembly='hg38'  # Supports only Hg38 genome assembly
Ref='/home/vinodsingh/mutationalscanning/Workspaces/vinod/reference/hg38/hg38.fasta'


# Define directories
projDir='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/Proj_GermReadsClassification/'
codeDir="${projDir}/src/"
dataName='ob2v'
baseDir='/home/vinodsingh/mutationalscanning/Workspaces/vinod/'
DataDir="${baseDir}/${dataName}_chrBAMs_${GenomeAssembly}_test/"

# Train data 
trainData="${projDir}/test_train.bed.gz"


# Define BAM file path
bam=${DataDir}/${dataName}_${GenomeAssembly}.bam

# Define output directory and create it
OUTDIR="${DataDir}/${dataName}_BinoClassificationResults/"
mkdir -p $OUTDIR
cd $OUTDIR
```

## Step 2: Split BAM files

- Runs the BAM splitting script as a SLURM job.
- Lists output files and removes temporary chunk BAMs.
- Verifies completion and moves logs to the log directory.

```bash
CONDAenv_00='GermReadsClassification_NGSenv'
script_00='00_createDir_and_split_bam.sh'
OUTDIR_00="${OUTDIR}/00_BAMChunks/"
mkdir $OUTDIR_00
logs_00=$OUTDIR_00/logs
mkdir $logs_00
NCPUS_00=24

sbatch -J ${dataName}_splitBAM --time 10:00:00 --cpus-per-task $NCPUS_00 --mem=100GB --account $HPC_accountName \
  ${codeDir}/runOn_masterCluster.sh $CONDAenv_00 \
  $codeDir/${script_00} \
  -B $bam \
  -S $dataName \
  -G $GenomeAssembly \
  -R $Ref \
  -O $OUTDIR_00 \
  -N $NCPUS_00

# Check output and cleanup
ls -halt $OUTDIR_00
ls -halt $OUTDIR_00/${dataName}_${GenomeAssembly}_bamChunks_chr1

# Verify job completion and move logs
grep -E  "===Task Completed==" $(ls  xx-${dataName}_splitBAM*.out) 
mv xx-${dataName}_splitBAM*.out $logs_00
rm $OUTDIR_00/${dataName}_${GenomeAssembly}_bamChunks_chr*/chunk*.bam


```

## Step 3: Convert BAM Chunks (1 million reads) to RDS Objects

- Converts BAM chunks into RDS objects with binary strings for ML and MM tags.
- Submits parallel jobs for each chromosome.
- Verifies task completion and checks for errors.

```bash
CONDAenv_01='GermReadsClassification_Renv'
script_01='01_BAM2RDS_ByChunk.R'
OUTDIR_01="${OUTDIR}/01_BAM_RDS_ByChunk/"
mkdir $OUTDIR_01
logs_01=$OUTDIR_01/logs
mkdir $logs_01
NCPUS_01=10

parallel -v sbatch -J ${dataName}_BAM2RDS_ByChunk_chr{} --time 08:00:00 --cpus-per-task $NCPUS_01 --mem=100GB --account $HPC_accountName \
  ${codeDir}/runOn_masterCluster.sh $CONDAenv_01 \
  ${codeDir}/${script_01} \
  -I ${OUTDIR_00}/${dataName}_${GenomeAssembly}_bamChunks_chr{} \
  -N $NCPUS_01 \
  -O ${OUTDIR_01}/BAMRDS_chr{} ::: $(seq 1 22) X Y
  
# Verify logs and errors
grep -E "===Task Completed==" $(ls xx-${dataName}_BAM2RDS_ByChunk_chr*.out) | wc -l # should be 24
grep -E "error|Error" $(ls xx-${dataName}_BAM2RDS_ByChunk_chr*.err)

# Move logs
mv xx-${dataName}_BAM2RDS_ByChunk_chr* $logs_01/
```


## Step 4: Get Methylation value for each CpG site of the read 

- chunks of RDS objects are converted to bed files
- bed file has methylation data for each CpG site on a reads
- Approx Runtime for a chunk (of 1 million long reads) with 10 CPUs: 8 MIN

```bash
CONDAenv_02='GermReadsClassification_Renv'
script_02='02_MethyperCpG_forReads_fromHiFi_ByChunk.R'
OUTDIR_02="${OUTDIR}/02_MethyperCpG_forReads_ByChunk/"
mkdir $OUTDIR_02
logs_02=$OUTDIR_02/logs
mkdir $logs_02
NCPUS_02=12

parallel  sbatch -J ${dataName}_MethyperCpG_forReads_ByChunk_chr{} --time 15:00:00 --cpus-per-task $NCPUS_02 --mem=80GB --account $HPC_accountName  \
  ${codeDir}/runOn_masterCluster.sh $CONDAenv_02 \
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
grep -E ">> completed" $(ls xx-${dataName}_MethyperCpG_forReads_ByChunk_chr*.*.out) | sort -V | wc -l # should be 24

# save logs
mv xx-${dataName}_MethyperCpG_forReads_ByChunk_chr* $logs_02/

```

## Step 5: calculate methylated and unmethylated component prior

- Generate methylation frequency table from bam file.
- Use methylation frequency table to get  perameters of the prior i.e., methylated and unmethylated component
- Submits parallel jobs per chromosome.

```bash
CONDAenv_03a='GermReadsClassification_NGSenv'
script_03a="03a_Generate_HiFiMethyFreqTableF.sh"
MethyType='CpG'
OUTDIR_03=${OUTDIR}/03_${MethyType}_priorData/
mkdir $OUTDIR_03
logs3=$OUTDIR_03/logs
mkdir $logs3
NCPUS_03=24

sbatch -J xx-${dataName}_${MethyType}_freqTable --mem 140g --cpus-per-task=$NCPUS_03  --time=04:00:00 --account $HPC_accountName \
  ${codeDir}/runOn_masterCluster.sh $CONDAenv_03a \
  ${codeDir}/${script_03a}  \
  -B $bam  \
  -D $dataName  \
  -O ${OUTDIR_03}/${dataName}_HiFiMethyFreqTable.tsv \
  -N $NCPUS_03

```

- estimate prior parameters from methylated table #

```bash
CONDAenv_03b='GermReadsClassification_Renv'
script_03b="03b_GA_ParametersEstimation.R"
OUTDFILE_03="${OUTDIR_03}/Parameters.tsv"
NCPUS_03=1

sbatch -J ${dataName}_calPrior_parameters --mem 50GB --cpus-per-task=$NCPUS_03 --account $HPC_accountName --time 1:00:00 \
  ${codeDir}/runOn_masterCluster.sh $CONDAenv_03b \
  ${codeDir}/${script_03b} \
  -I ${OUTDIR_03}/${dataName}_HiFiMethyFreqTable.tsv \
  -D $dataName \
  -O $OUTDFILE_03 \
  -N $NCPUS_03

# save logs
mv xx-${dataName}_${MethyType}_freqTable* $logs3/
mv xx-${dataName}_calPrior_parameters* $logs3/
    
  
```

## Step 6: Calculate  read's log-Likelihood in each Cell Types

- Computes read likelihood in different cell types.
- Runs likelihood calculations for each chromosome in parallel.

```bash
CONDAenv_04='GermReadsClassification_Renv'
script_04='04_Reads_likelihood_InSpermCellTypes_ByChunk.R'
OUTDIR_04="${OUTDIR}/04_HiFiReadsLL_InCellTypes_ByChunk/"
OUTFILE_04=${OUTDIR_04}/${dataName}_${GenomeAssembly}_allchr_LL_of_readInCellType.tsv # final output file witl LL of all chromsomes.
mkdir $OUTDIR_04
logs_04=$OUTDIR_04/logs
mkdir $logs_04
NCPUS_04=10

parallel sbatch -J ${dataName}.ReadsLL_ByChunk_chr{} --time 08:00:00 --cpus-per-task $NCPUS_04 --mem=150GB --account $HPC_accountName \
  ${codeDir}/runOn_masterCluster.sh $CONDAenv_04 \
  ${codeDir}/${script_04} \
  -G $GenomeAssembly \
  -I ${OUTDIR_02}/'MethyperCpG_forReads_'chr{} \
  -C chr{} \
  -D $dataName \
  -T $trainData \
  -P $OUTDFILE_03 \
  -O ${OUTDIR_04}/'Reads_likelihood_'chr{} \
  -N $NCPUS_04 ::: $(seq 1 22) X Y
  
# varify run#
grep -E "error|Error" $(ls xx-${dataName}.ReadsLL_ByChunk_chr*.*.err)
grep -E ">> completed" $(ls xx-${dataName}.ReadsLL_ByChunk_chr*.*.out) | sort -V
grep -E ">> completed" $(ls xx-${dataName}.ReadsLL_ByChunk_chr*.*.out) | sort -V | wc -l # should be 24
cat xx-${dataName}.ReadsLL_ByChunk_chr*.out | grep -E ">>germ frac:"

# save logs
mv xx-${dataName}.ReadsLL_ByChunk_chr*.*.* $logs_04/


```

## Step 7: Merge all chromosomes LL datasets 

- Merges all chromosome-wise likelihood files into a single output

```bash
LL_ChrFiles=($(ls -1 ${OUTDIR_04}/${dataName}_${GenomeAssembly}_chr*_LL_of_readInCellType.tsv | sort -V)) # chromosome LL files 
if [ ${#LL_ChrFiles[@]} == 24 ]; then
  echo -e ">> \xE2\x9C\x85 RUN is successful on all chromosome"
  { head -n 1 ${LL_ChrFiles[0]}; tail -n +2 -q ${LL_ChrFiles[*]}; } > $OUTFILE_04
else
  echo -e ">> \xE2\x9D\x8C RUN is not successful on all chromosome"
fi

head $OUTFILE_04

```

## Step 8: Calculate Germline Fraction

- Find germline fraction in the sample

```bash
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

```

## Step 9: Cleanup intermediate files

```bash
rm 00_BAMChunks/${dataName}_hg38_bamChunks_chr*/sorted.chunk*.bam
rm 01_BAM_RDS_ByChunk/BAMRDS_chr*/sorted.chunk*.rds
rm -r 04_HiFiReadsLL_InCellTypes_ByChunk/Reads_likelihood_chr*
rm 02_MethyperCpG_forReads_ByChunk/MethyperCpG_forReads_chr*/muensterCpGovl_chunk*.RDS
rm 02_MethyperCpG_forReads_ByChunk/MethyperCpG_forReads_chr*/ReadsMethyPerCpG_chunk*.rds
rm -r 04_HiFiReadsLL_InCellTypes_ByChunk/Reads_likelihood_chr*

```

## Training dataset format

- We used BWA-METH pipeline to get per-base CpG methylation metric from EM-Seq / BS-Seq fastq files used for training our model. https://github.com/vinodsinghjnu/BWA_METH_pipeline
- Format the output in desired format of "test_train.bed.gz"

---
This Markdown file contains all the required Bash scripts for running the HiFi Reads Classification Pipeline. Make sure to execute them in order within an SLURM-managed HPC environment.
