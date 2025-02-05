library("optparse");library(tictoc); library(parallel); library(plyranges)


option_list = list(
  make_option(c("-G", "--GenomeAssembly"), type="character", default=NULL,
              help="Human Genome assembly name", metavar="character"),
  make_option(c("-I", "--INPUTFILEDIR"), type="character", default=NULL,
              help="Path of input reads bed File", metavar="character"),
  make_option(c("-O", "--OUTDIR"), type="character", default=NULL,
              help="Output dir", metavar="character"),
  make_option(c("-C", "--CHROMOSOME"), type="character", default=NULL,
              help="name of the chromosome in NCBI format", metavar="character"),
  make_option(c("-D", "--DATANAME"), type="character", default=NULL,
              help="name of the dataset", metavar="character"),
  make_option(c("-T", "--TRAINDATA"), type="character", default=NULL,
              help="Traning data set (rds file)", metavar="character"),
  make_option(c("-P", "--PARAMETERS"), type="character", default=NULL,
              help="Parameters TSV file, Cols1: alpha, col2: beta, col3:weight, for methylated comp ", metavar="character"),
  make_option(c("-N", "--NCORES"), type="integer", default=NULL,
              help="number of CPU cores to use", metavar="number")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$INPUTFILEDIR)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

## INPUT parameters ##
GenomeAssembly=opt$GenomeAssembly
INPUTFILEDIR=opt$INPUTFILEDIR
OUTDIR=opt$OUTDIR
CHR=opt$CHROMOSOME
DATANAME=opt$DATANAME
TRAINDATA=opt$TRAINDATA
PARAMETERS=opt$PARAMETERS
NCORES=opt$NCORES

cat(">>===INPUT parameters===\n\n")
print(paste0('GenomeAssembly: ', GenomeAssembly))
print(paste0('INPUTFILEDIR: (chunk rds files with methylated CpG sites on each read)', INPUTFILEDIR))
print(paste0('OUTDIR: Output dir: ', OUTDIR))
print(paste0('CHR: Chromosome: ', CHR))
print(paste0('DATANAME: Dataname: ', DATANAME))
print(paste0('TRAINDATA: Training data file: ', TRAINDATA))
print(paste0('PARAMETERS: file: ', PARAMETERS))
print(paste0('NCORES: Number of cores: ', NCORES))


args <- commandArgs()
print(paste(c('command:', args), collapse=' '))
cat("<<======\n\n")


# Check if OUTPUT directory exists
if (dir.exists(OUTDIR)) {
  # If the directory exists, delete it and its contents
  unlink(OUTDIR, recursive = TRUE)
}
# Create a new directory
dir.create(OUTDIR)

## Read beta shape parameters for methylated component ##

#betaDist.shapePara=as.data.frame(readxl::read_excel('/home/vinodsingh/mutationalscanning/Workspaces/vinod/CpG_binomial_trainData/HiFidata_MethyBetaShapeParameters_AllSamples.xlsx'))
#rownames(betaDist.shapePara) <- betaDist.shapePara[[1]]
#betaDist.shapePara <- betaDist.shapePara[ , -1]
betaDist.shapePara=read.table(PARAMETERS)
print(betaDist.shapePara)
shapePara=betaDist.shapePara[DATANAME,]


### FUNCTIONS ##

#library(fastmatch)
`%fin%` <- function(x, table) {
  #stopifnot(require(fastmatch))
  fastmatch::fmatch(x, table, nomatch = 0L) > 0L
}



## Likelihood function #
LL_forEachCellType=function(r_data,x, f_muenster_cpg_data_chunk)
{
  cpgMethy_locs=unlist(strsplit(r_data$methyloc, ",")) # cpg locations on the read
  cpgMethy_vals=as.numeric(unlist(strsplit(r_data$methyVal, ",")))/256
  names(cpgMethy_vals)=cpgMethy_locs



  readsCpGLocs_in_genomeCpGdata=cpgMethy_locs[cpgMethy_locs  %fin% names(f_muenster_cpg_data_chunk) ]
  #print(readsCpGLocs_in_genomeCpGdata)
  if(length(readsCpGLocs_in_genomeCpGdata)==0){
    ll_allCpG_ofRead=c(r_data$readname, r_data$chr, rep(0, length(names(mcols(f_muenster_cpg_data_chunk)))), 0 , 0)
    return(ll_allCpG_ofRead)
  }

  rm(cpgMethy_locs)

  read_CpGMethyCounts_fromTrain=f_muenster_cpg_data_chunk[readsCpGLocs_in_genomeCpGdata]
  read_CpGMethyCounts_fromTrain$MethyProb_fromHiFi=cpgMethy_vals[names(read_CpGMethyCounts_fromTrain)]
  rm(cpgMethy_vals)

  CpGs_onRead=length(read_CpGMethyCounts_fromTrain)
  badCpGs_onRead=NA
  

  DataNames=colnames(mcols(read_CpGMethyCounts_fromTrain))[c(1,2,3,4,5,6,7)]
  ll_allCpG_ofRead=c()
  ll_allCpG_ofRead[1]=r_data$readname
  ll_allCpG_ofRead[2]=r_data$chr
  p=read_CpGMethyCounts_fromTrain$MethyProb_fromHiFi # methy P for the read
  p[which(p==0)]= .Machine$double.eps
  rm(r_data)
  
  # likelihood on reads Methyvalues in each celltypes
  #for(d in c(1,2,3,4)){
  for(d in seq(1,length(mcols(read_CpGMethyCounts_fromTrain))-1)){

    MethyCountData=mcols(read_CpGMethyCounts_fromTrain)[,d]

    mC_cnts=sapply(MethyCountData, function(x) x[2])
    umC_cnts=sapply(MethyCountData, function(x) x[3])
    rm(MethyCountData)
    r_cnts=mC_cnts+umC_cnts

    # finalised after agreement ##
    ll_perCpG=log((dbeta(p,shapePara[,'ms1'],shapePara[,'ms2']) * ((mC_cnts + shapePara[,'mWeight']) / (r_cnts +1)))+ (dbeta(p,shapePara[,'ms2'],shapePara[,'ms1']) * ((umC_cnts + (1-shapePara[,'mWeight'])) / (r_cnts +1))))
    ##((1/v) * (1-p)^(v-1) * (mC_cnts +1) / (r_cnts +2)) + ((1/v) * (p)^(v-1) * (umC_cnts +1) / (r_cnts +2))


    ll_allCpG_ofRead[d+2]=sum(ll_perCpG)

  }

  ll_allCpG_ofRead[d+3]=CpGs_onRead
  ll_allCpG_ofRead[d+4]=badCpGs_onRead
  ll_allCpG_ofRead
}

## Read chunk data and do lilelihood calculations for each read in given cell types ##
calculate_chunk=function(idx){
  chunk_data=readRDS(f[idx])
  print(f[idx])
  #print(head(chunk_data))
  c_size=nrow(chunk_data)
  muensterPath <- gsub("ReadsMethyPerCpG", "muensterCpGovl", f[idx])
  muensterPath <- sub("\\.rds$", ".RDS", muensterPath)
  print(muensterPath)
  muenster_cpg_data_chunk=readRDS(muensterPath)
  #muenster_cpg_data_chunk=readRDS(paste0(INPUTFILEDIR,'/',sub(x=basename(f[idx]),'ReadsMethyPerCpG_chunk', 'muensterCpGovl_chunk')))

  # Perform likelihood calculations
  chunk.LL_of_readInCellType_list=lapply(seq(1, c_size), function(x) LL_forEachCellType(chunk_data[x,],x, muenster_cpg_data_chunk))
  rm(chunk_data)
  rm(muenster_cpg_data_chunk)
  chunk.LL_of_readInCellType_list_df=do.call('rbind', chunk.LL_of_readInCellType_list)
  #print(head(chunk.read_perCpG_methy_df))
  rm(chunk.LL_of_readInCellType_list)

  saveRDS(chunk.LL_of_readInCellType_list_df, file = paste0(OUTDIR, '/' ,"ReadsLL_chunk", idx, "_out.RDS"))
  d=nrow(chunk.LL_of_readInCellType_list_df)
  rm(chunk.LL_of_readInCellType_list_df)
  gc()
  print(paste0('>Memory Used (GB) : ',signif(pryr::mem_used()/(10^9),3) ))

  if(d == c_size){
    #return("success!")
    print(paste0('success: Chunk ',idx))
  }else{
    stop("error in this process!")    
  }

  return(NULL)

}


#### MAIN FUNCTION ###

### READ TRAINING DATA ##
#TRAINDATA='/home/vinodsingh/mutationalscanning/Workspaces/vinod/CpG_binomial_trainData/'
#TRAINDATA='/home/vinodsingh/mutationalscanning/Workspaces/vinod/CpG_binomial_trainData/sorted_Cov.g5_Binomial_parameters_allCellTypes.rds'
print(CHR)
muenster_cpg_data=readRDS(TRAINDATA) %>% plyranges::filter(seqnames == CHR)
print(paste0('>TRAINDATA:',TRAINDATA))


## list INPUT chunk files  ##
f <- gtools::mixedsort(Sys.glob(paste0(INPUTFILEDIR,"/ReadsMethyPerCpG_chunk*.rds")))
print(f)
reads_data=lapply(f, function(x) readRDS(x))


## Split "muenster_cpg_data" data  within ranges of each chunk, to reduce RAM usage while parallel processing #
max_locs=sapply(seq(1,length(f)), function(x) max(c(reads_data[[x]]$start,reads_data[[x]]$end)))
min_locs=sapply(seq(1,length(f)), function(x) min(c(reads_data[[x]]$start,reads_data[[x]]$end)))

w2=mclapply(seq(1,length(f)), function(x) saveRDS(muenster_cpg_data %>% filter(seqnames == CHR, start >= min_locs[x], end <= max_locs[x]), file= paste0(INPUTFILEDIR,'/',"muensterCpGovl_chunk", x, ".RDS") ),  mc.cores = NCORES, mc.preschedule = FALSE)
print(w2)
mcolsNames_muenster_cpg_data=names(mcols(muenster_cpg_data))
rm(muenster_cpg_data)

## Store reads with warnings in a seperate variable named warnings_data ## 
## Get reads with CpG ovlerlap less than 90% of ref CpGs, as numeric values ##
all_reads_data=as.data.frame(do.call(rbind, lapply(f, function(x) readRDS(x))))
print(head(all_reads_data))
badOvls=which(all_reads_data$warning!='-')
no_of_Reads=nrow(all_reads_data)
badOvls_data=sapply(all_reads_data$warning[badOvls], function(x) as.numeric(strsplit(x, split=' ')[[1]][3]))
names(badOvls_data)=NULL
warnings_data=all_reads_data$warning
warnings_data[badOvls]=badOvls_data
rm(badOvls)
rm(badOvls_data)
rm(all_reads_data)
##
rm(reads_data)
gc()
##



# Parallel LL caluculation on Chunk ##
ret_dim=mclapply(seq(1,length(f)), function(x) calculate_chunk(x), mc.cores = NCORES,  mc.preschedule = FALSE, mc.cleanup = TRUE )
print(ret_dim)
rm(ret_dim)



# read chromosome chunks liklihood data and merge, Write data of each chromosome in tsv file##
f2 <- gtools::mixedsort(Sys.glob(paste0(OUTDIR,"/ReadsLL_chunk*_out.RDS")))
LL_of_readInCellType_list_df_merge=do.call('rbind', lapply(f2,  readRDS))
LL_of_readInCellType_list_df_merge=as.data.frame(LL_of_readInCellType_list_df_merge)
LL_of_readInCellType_list_df_merge=cbind(LL_of_readInCellType_list_df_merge, warnings_data)
colnames(LL_of_readInCellType_list_df_merge)=c('readsID','chr', mcolsNames_muenster_cpg_data, 'CpGs_onRead','badCpGs_onRead', 'refCpGOvl_perc' )

# convert columns to numeric #
for(i in seq(3,ncol(LL_of_readInCellType_list_df_merge)-2))
{
  LL_of_readInCellType_list_df_merge[,i]=as.numeric(LL_of_readInCellType_list_df_merge[,i])
}

if(nrow(LL_of_readInCellType_list_df_merge)!=no_of_Reads){
  print(paste0('> Warning: Reads < Input reads: ', nrow(LL_of_readInCellType_list_df_merge), ' of ', no_of_Reads ))
  save.image(paste0(dirname(OUTDIR), 'wksp_LL_of_readInCellType.RDATA'))
}


badCpGs_ind=which(LL_of_readInCellType_list_df_merge$CpGs_onRead==0)
if(length(badCpGs_ind)>0){
  LL_of_readInCellType_list_df_merge=LL_of_readInCellType_list_df_merge[-badCpGs_ind,]
}

outFile=paste0(dirname(OUTDIR),'/',DATANAME,'_',GenomeAssembly,'_',CHR,'_LL_of_readInCellType.rds')
print(paste0('>> saved to: ',outFile))
saveRDS(LL_of_readInCellType_list_df_merge, outFile)

outFile2=paste0(dirname(OUTDIR),'/',DATANAME,'_',GenomeAssembly,'_',CHR,'_LL_of_readInCellType.tsv')
write.table(LL_of_readInCellType_list_df_merge, file = outFile2,sep = '\t', quote = F, append = F, row.names = F)


temp_df=LL_of_readInCellType_list_df_merge[,seq(3,9)]
mx=sapply(seq(1,nrow(temp_df)), function(x) which.max(temp_df[x,]))
maxind_tables=table(mx)
print(paste0('>>germ frac: ', sum(maxind_tables[c(1,2,3,4,5)])/sum(maxind_tables)))


print(">> completed")


#r=readRDS(paste0(OUTDIR,'/ob2v_hg38_chr16_LL_of_readInCellType.rds'))
