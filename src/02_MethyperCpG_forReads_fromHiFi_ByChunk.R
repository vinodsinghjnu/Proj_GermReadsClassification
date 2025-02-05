library("optparse")

option_list = list(
  make_option(c("-G", "--GenomeAssembly"), type="character", default=NULL,
              help="genome assembly name", metavar="character"),
  make_option(c("-I", "--INPUTDIR"), type="character", default=NULL,
              help="BAMRDS chunks dir", metavar="character"),
  make_option(c("-C", "--CHROMOSOME"), type="character", default=NULL,
              help="name of the chromosome in NCBI format", metavar="character"),
  make_option(c("-D", "--DATANAME"), type="character", default=NULL,
              help="name of the chromosome in NCBI format", metavar="character"),
  make_option(c("-O", "--OUTDIR"), type="character", default=NULL,
              help="OUT data directory", metavar="character"),
  make_option(c("-N", "--NCORES"), type="integer", default=NULL,
              help="number of CPU cores to use", metavar="number")
);

library('Rsamtools', quietly = T);library(GenomicAlignments,quietly = T);library(rtracklayer,quietly = T); library(logr,quietly = T)
library(BSgenome.Hsapiens.UCSC.hg38,quietly = T); library(Biostrings,quietly = T); library(parallel,quietly = T); library(BSgenome,quietly = T)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt)<1){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


cat(">>======\n\n")

# description ##

GenomeAssembly=opt$GenomeAssembly
INPUTDIR=opt$INPUTDIR
OUTDIR=opt$OUTDIR
CHR=opt$CHROMOSOME
DATANAME=opt$DATANAME
NCORES=opt$NCORES

print(paste0('GenomeAssembly: ', GenomeAssembly))
print(paste0('BAMRDS: INPUT BAMRDS chunk dir: ', INPUTDIR))
print(paste0('CHROMOSOME: ', CHR))
print(paste0('OUTDIR: Output dir: ', OUTDIR))
print(paste0('DATANAME: Name of the sample: ', DATANAME))
print(paste0('NCORES: Number of cores: ', NCORES))

args <- commandArgs()
print(paste(c('command:', args), collapse=' '))
cat("<<======\n\n")


# Check if the directory exists
if (dir.exists(OUTDIR)) {
  # If the directory exists, delete it and its contents
  unlink(OUTDIR, recursive = TRUE)
}
# Create a new directory
dir.create(OUTDIR)


if(GenomeAssembly=='hg38'){
  BS_genome=BSgenome.Hsapiens.UCSC.hg38::Hsapiens
} else if(GenomeAssembly=='hg19'){
  BS_genome=BSgenome.Hsapiens.UCSC.hg19::Hsapiens
}


## fast ## # operations refrence: https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
convert_query_to_reference <- function(query_vector, cigar_string) {
  # Split cigar into numbers and operations using regex
  cigar_elements <- unlist(strsplit(cigar_string, "(?<=\\D)", perl = TRUE))

  # Preallocate for reference vector to avoid repeated memory allocation
  reference_vector <- vector("numeric", length(query_vector))
  ref_index <- 1
  query_index <- 1

  for (element in cigar_elements) {
    number <- as.numeric(gsub("\\D", "", element))
    operation <- gsub("\\d", "", element)

    if (operation %in% c("M", "=")) {
      # Match: Copy query sequence to ref sequence ( comsumes both query and ref)
      reference_vector[ref_index:(ref_index + number - 1)] <- query_vector[query_index:(query_index + number - 1)]
      ref_index <- ref_index + number
      query_index <- query_index + number
    } else if (operation == "D") {
      # Deletion: Insert gaps by NA in the reference sequence
      ref_index <- ref_index + number
    } else if (operation == "I") {
      # Insertion: Skip these positions in the query sequence
      query_index <- query_index + number
    } else if (operation == "X") {
      # Mismatch: Insert mismatched elements ( comsumes both query and ref)
      reference_vector[ref_index:(ref_index + number - 1)] <- query_vector[query_index:(query_index + number - 1)]
      ref_index <- ref_index + number
      query_index <- query_index + number
    } else if (operation == "S") {
      # Soft clip: Skip clipped bases in the query sequence, no change in ref sequence
      query_index <- query_index + number
    }
  }

  # Trim excess NA values
  return(reference_vector[1:(ref_index - 1)])
}


getMethylated_GpGs_data=function(aln_read)
{
  #x=1
  #print(paste0('> ',x,'/',l))

  #aln_read=chr_MethyReads[x]
  names(aln_read)=aln_read@elementMetadata@listData$qname
  #aln_read@elementMetadata@listData$ML[[1]]
  #aln_read@elementMetadata@listData$MM_rle[[1]]
  #aln_read@elementMetadata@listData$CnCG_ind_rle.list[[1]]$c.ind_rle
  #aln_read@elementMetadata@listData$CnCG_ind_rle.list[[1]]$cg.ind_rle
  st=start(aln_read)
  en=end(aln_read)
  strand=as.vector(aln_read@strand)
  warn='-'

  if(strand=='+'){
    read_seq=aln_read@elementMetadata@listData$seq
    #read_CGind=start(matchPattern('CG',as.vector(read_seq)))
    #read_Cind=start(matchPattern('C',as.vector(read_seq)))

  }else if(strand=='-'){
    read_seq=reverseComplement(aln_read@elementMetadata@listData$seq)
    #read_CGind=end(matchPattern('CG',as.vector(read_seq))) # location on ref strand
    #read_Cind=start(matchPattern('G',as.vector(read_seq))) # location on ref strand
  }

  #ref_CGind=start(matchPattern('CG',BS_genome$chr16[st:en]))
  #ref_Cind=start(matchPattern('C',BS_genome$chr16[st:en]))


  #### Refcords for MM tag ##
  c_rle=aln_read@elementMetadata@listData$CnCG_ind_rle.list[[1]]$c.ind_rle
  #length(intersect(read_Cind,which(c_rle==1)))/ length(which(c_rle==1)) ## all reads C ovl read  Cs
  MM_rle=aln_read@elementMetadata@listData$MM_rle[[1]]
  MethyC_loc=which(c_rle==1)[which(MM_rle==1)]

  if(strand=='-'){
    MethyC_loc=rev(width(read_seq)-MethyC_loc+1)
  }
  #intersect(read_Cind,MethyC_loc)
  #intersect(read_CGind,MethyC_loc)

  methy_c_rle=rep(0,length(c_rle)); methy_c_rle[MethyC_loc]=1; methy_c_rle=Rle(methy_c_rle)
  methy_cg=paste0(methy_c_rle,collapse='')
  MethyCpGs_string_ref=sequenceLayer(BStringSet(methy_cg), cigar(aln_read))
  MethyCpGs_locs_ref=gregexpr(MethyCpGs_string_ref,pattern = '1')[[1]]
  if(strand=='+'){
    MethyCpG_locs_refCord=MethyCpGs_locs_ref+st-1
  }else if(strand=='-'){
    MethyCpG_locs_refCord=MethyCpGs_locs_ref+st-2
  }

  p=plyranges::as_granges(data.frame(seqnames=seqnames(aln_read),start=MethyCpG_locs_refCord, end=MethyCpG_locs_refCord+1), strand=strand(aln_read))
  p$seq=getSeq(BS_genome,p)
  t=table(p$seq)

  if(is.element('CG', names(t))){
    if( t['CG']/sum(t) < .90)
    {
      #print(">> methylated CpGs dont Map to CpGs sites on refrence" )
      #print('Dinucs Freq at Methylated sites, may be most of the sites are overlapping with SNPs ')
      #print(paste0('>Read: ',x,' ',names(aln_read)) )
      #print(t)
      #warning(paste0(">Read: CG% ",signif(t['CG']/sum(t),2),' ',names(aln_read), '>> methylated CpGs dont Map to CpGs sites on refrence, may be overlapping to SNPs'))
      warn=paste0(">Read: CG% ",signif(t['CG']/sum(t),2),' ',names(aln_read), '>> methylated CpGs dont Map to CpGs sites on refrence, may be overlapping to SNPs')

      #print(x)
      #print(aln_read)
      #stop()
    }
  }else{
    warn=paste0('> No CGs in the read')
    d=data.frame(readname=names(aln_read), chr=as.vector(seqnames(aln_read)) ,start=st, end=en, strand=strand, methyloc='-',methyVal='-',  warning=warn )
    return(d)
  }

  #### methy values for Refcords ##

  ML=aln_read@elementMetadata@listData$ML[[1]]
  ML_ind_alongseq=Rle(rep(0,width(read_seq)))
  if(strand=='+'){
    ML_ind_alongseq[MethyC_loc]=ML
  }else if(strand=='-'){
    ML_ind_alongseq[rev(MethyC_loc)]=ML
  }

  methyVals_ref=convert_query_to_reference(as.vector(ML_ind_alongseq), cigar(aln_read))
  methyVals_ref=methyVals_ref[MethyCpGs_locs_ref]
  names(methyVals_ref)=paste0(seqnames(aln_read),'_',MethyCpG_locs_refCord)
  df=data.frame(readname=names(aln_read), chr=as.vector(seqnames(aln_read)) ,start=st, end=en, strand=strand, methyloc=paste(names(methyVals_ref),collapse=','),methyVal=paste(methyVals_ref,collapse=','),  warning=warn )

}


#f <- gtools::mixedsort(Sys.glob(paste0(INPUTDIR,paste0('/*_',CHR,'_chunk*.rds'))))
f <- gtools::mixedsort(Sys.glob(paste0(INPUTDIR,paste0('/sorted.chunk*.rds'))))
print(f)

nreads=0
for(idx in 1:length(f))
{
  print(paste0("---- chunk: ", idx,' of ',length(f), "------"))
  chunk_data=readRDS(f[idx])
  nreads=nreads+length(chunk_data)
  print(f[idx])
  print(head(chunk_data))
  chunk.read_perCpG_methy_list=mclapply(seq(1, length(chunk_data)), function(x) getMethylated_GpGs_data(chunk_data[x]), mc.cores = getOption("mc.cores", NCORES))
  rm(chunk_data)
  chunk.read_perCpG_methy_df=do.call('rbind', chunk.read_perCpG_methy_list)
  emptyReads_ind=which(chunk.read_perCpG_methy_df$methyloc=='-')
  if(length(emptyReads_ind)>0){
    chunk.read_perCpG_methy_df = chunk.read_perCpG_methy_df[-emptyReads_ind,]
  }

  #print(head(chunk.read_perCpG_methy_df))
  rm(chunk.read_perCpG_methy_list)
  head(chunk.read_perCpG_methy_df)
  #saveRDS(chunk.read_perCpG_methy_df, file = paste0(rds_Dir, '/' ,"chunk_out", idx, ".RDS"))

  #saveRDS(chunk.read_perCpG_methy_df, file = sub(x=f[idx],'_chunk','_ReadsMethyPerCpG_chunk' ))
  saveRDS(chunk.read_perCpG_methy_df, file = paste0(OUTDIR,'/', sub(x=basename(f[idx]),'sorted.chunk','ReadsMethyPerCpG_chunk' )))

  rm(chunk.read_perCpG_methy_df)
  gc()
  print(paste0("====Finished chunk: ", idx, "======="))
}


outFile=paste0(OUTDIR,'/',DATANAME,'_',GenomeAssembly,'_',CHR,'_ReadsMethyPerCpG.bed')

f2 <- gtools::mixedsort(Sys.glob(paste0(OUTDIR,"/ReadsMethyPerCpG_chunk*.rds")))
print(f2)
#f2=stringr::str_sort(f2, numeric = TRUE)
read_perCpG_methy_df_merge=do.call('rbind', lapply(f2,  readRDS))
colnames(read_perCpG_methy_df_merge)=c('readname','chr', 'start', 'end', 'strand', 'methyloc', 'methyVal', 'warning' )
write.table(read_perCpG_methy_df_merge, file = outFile,sep = '\t', quote = F, append = F, row.names = F)

if(nrow(read_perCpG_methy_df_merge)!=nreads){
  print(paste0('>> Warning: Reads < Input reads: ', nrow(read_perCpG_methy_df_merge), ' of ', nreads ))
}

if(nrow(read_perCpG_methy_df_merge)==nreads)
{
  print('>> Successful: Same reads in INPUT and OUTPUT')
}

print(">> completed")
