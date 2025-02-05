library("optparse")

option_list = list(
  make_option(c("-I", "--BAMDIR"), type="character", default=NULL,
              help="BAM files directory of a chromosome", metavar="character"),
  make_option(c("-N", "--NCPUS"), type="integer", default=5,
              help="CPU", metavar="CPU required: integer value"),
  make_option(c("-O", "--OUTDIR"), type="character", default=NULL,
              help="OUT data directory", metavar="character")
);
## Libraries required##
library('Rsamtools', quietly = T);library(GenomicAlignments,quietly = T);library(rtracklayer,quietly = T); library(logr,quietly = T)
library('parallel',quietly = T)

## INPUTS ##
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt)<1){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


cat(">>======\n\n")

# description ##

BAMDIR=opt$BAMDIR
OUTDIR=opt$OUTDIR
NCPUS=opt$NCPUS


print(paste0('BAMDIR: INPUT BAM dir: ', BAMDIR))
print(paste0('NCPUS: ', NCPUS))
print(paste0('OUTDIR: Output dir: ', OUTDIR))

args <- commandArgs()
print(paste(c('command:', args), collapse=' '))
cat("<<======\n\n")

##

# Del if OUTDIR exists and create a new
if (dir.exists(OUTDIR)) {
  unlink(OUTDIR, recursive = TRUE)
}
dir.create(OUTDIR)


### FUNCTIONS ###
convert2rle=function(v){  #
  #v=c(5,12,0)
  v_m=as.vector((sapply(v, function(x) c(x,1))))
  rle_v=Rle(rep(c(0,1), floor(length(v_m)/2)) , v_m)
  #which(rle_v==1)
  rle_v[which(rle_v != 0)]=1
  return(rle_v)

}

## Convert pattern index to rle object ##
convert_patLoc2rle=function(x, bam.gr_qwidth_f, C_loc.list_f, CG_loc.list_f  ){
  c.ind=rep(0,bam.gr_qwidth_f[x])
  cg.ind=c.ind

  c.ind[C_loc.list_f[[x]]] = 1
  cg.ind[CG_loc.list_f[[x]]] = 1

  return(CnCG_ind.rle=list(c.ind_rle=Rle(c.ind), cg.ind_rle=Rle(cg.ind)))
}



operation_on_chunk=function(chunk_file){

    print(paste0("===Chunk",chunk_file," SAM data=="))

    print('==BAM stats==')
    test=idxstatsBam(chunk_file)
    #print(test) #  
    print(test[which(test[,'mapped']>0),])
    print("==BAM file seqinfo==")
    quickBamFlagSummary(chunk_file)

    bam.gr_f <- readGAlignments(chunk_file)
    aln.methy_f <- scanBam(chunk_file, param = param)

    aln.methy_first <- aln.methy_f[[1]]
    names(aln.methy_first)
    rm(aln.methy_f)
    gc()
    #str(aln.methy_first[["tag"]])
    #str(aln.methy_first[["tag"]][["ML"]])

    # aln.methy_first$tag$MM[1]
    # aln.methy_first$tag$ML[[1]]
    # aln.methy_first$tag$ML[2][[1]]

    MM_tag_list=lapply(seq(1, length(aln.methy_first$tag$ML)), function(x) as.numeric(unlist(strsplit(aln.methy_first$tag$MM[x], "\\,|\\;| "))[-1] ) )
    MM_tag_name=sapply(seq(1, length(aln.methy_first$tag$ML)), function(x) unlist(strsplit(aln.methy_first$tag$MM[x], "\\,|\\;| "))[1])

    MM_tag_name_types1=sapply(MM_tag_name, function(x) length(grep( 'C+', x)))# new addition 18 Apr 2024
    if(any(MM_tag_name_types1>1)){print(MM_tag_name_types1);stop("> MM_tag_name_types1: This is an error message: Multiple C+ base modificationfound in any read")} # new addition 18 Apr 2024
    MM_tag_name_types2=sapply(MM_tag_name, function(x) length(grep( 'G-', x))) # new addition 18 Apr 2024
    if(any(MM_tag_name_types2>1)){print(MM_tag_name_types2);stop("> MM_tag_name_types2: This is an error message: Multiple  C+ base modification found or they are at opposite strand")} # new addition 18 Apr 2024


    ML_tag_list=aln.methy_first[["tag"]][["ML"]]

    NMtag_perReads=sapply(MM_tag_list, length)
    badInd1=which(NMtag_perReads==0)
    #unlist(strsplit(aln.methy_first$tag$MM[1], "\\,|\\;| "))[1]
    badInd=which(is.na(MM_tag_name))
    badInd=c(badInd,badInd1)
    if(length(badInd)>0){
      #if( any(names(table(MM_tag_name[-badInd]))!='C+m') && any(names(table(MM_tag_name[-badInd]))!='C+m?')){stop("This is an error message: unidentified base modification found")}
      # table(names(table(MM_tag_name[-badInd])) %in% c('C+m','C+m?'))['FALSE']>0
      # any(names(table(names(table(MM_tag_name[-badInd])) %in% c('C+m','C+m?')))=='FALSE')
      if(any(names(table(names(table(MM_tag_name[-badInd])) %in% c('C+m','C+m?')))=='FALSE')){stop("This is an error message: unidentified base modification found")}
      if( !any( any(names(table(MM_tag_name[-badInd]))=='C+m') || any(names(table(MM_tag_name[-badInd]))=='C+m?'))){
        print(names(table(MM_tag_name[-badInd])))
        stop("This is an error message: Not available for this base modification")
      }
      MM_tag_list.filt=MM_tag_list[-badInd]
      MM_tag_list_rle.filt=lapply(seq(1,length(MM_tag_list.filt)), function(x) convert2rle(MM_tag_list.filt[[x]]))

      ML_tag_list.filt=ML_tag_list[-badInd]
    }else{
      MM_tag_list.filt=MM_tag_list
      MM_tag_list_rle.filt=lapply(seq(1,length(MM_tag_list.filt)), function(x) convert2rle(MM_tag_list.filt[[x]]))
      ML_tag_list.filt=ML_tag_list
    }

    print("==Print few tags==")
    print(head(MM_tag_name))
    print(head(MM_tag_list.filt,2))
    print(head(MM_tag_list_rle.filt[[1]]))

    # x=1
    #which(as.vector(MM_tag_list_rle.filt[[x]])==1)
    #ML_tag_list.filt[[x]]
    # length(ML_tag_list.filt[[x]])== length(which(as.vector(MM_tag_list_rle.filt[[x]])==1))

    #aln.methy_first$tag$MM[bandInd]

    ##

    rm(MM_tag_list)
    rm(MM_tag_name)
    rm(ML_tag_list)
    gc()

    ### gr of bam # passed in function
    print('>>>>>')
    print(head(bam.gr_f))

    bam.gr_f@elementMetadata@listData$seq = aln.methy_first$seq
    bam.gr_f@elementMetadata@listData$qual = aln.methy_first$qual
    bam.gr_f@elementMetadata@listData$qname = aln.methy_first$qname
    bam.gr_f@elementMetadata@listData$flag = aln.methy_first$flag
    bam.gr_f@elementMetadata@listData$mapq = aln.methy_first$mapq
    bam.gr_f@elementMetadata@listData$mapq = aln.methy_first$mapq

    rm(aln.methy_first)
    gc()

    # table(strand(bam.gr_f[which(strand(bam.gr_f)=='+')]))
    # table(strand(bam.gr_f[which(strand(bam.gr_f)=='-')]))
    bam.gr_f@elementMetadata@listData$seq[which(strand(bam.gr_f)=='-')]= reverseComplement(bam.gr_f@elementMetadata@listData$seq[which(strand(bam.gr_f)=='-')])

    # find C and GC locations on the read and convert to rle objects for data compression #
    C_loc.ranges=vmatchPattern('C',bam.gr_f@elementMetadata@listData$seq)
    C_loc.list=lapply(C_loc.ranges, function(x) x@start)
    rm(C_loc.ranges)

    CG_loc.ranges=vmatchPattern('CG',bam.gr_f@elementMetadata@listData$seq)
    CG_loc.list=lapply(CG_loc.ranges, function(x) x@start)
    rm(CG_loc.ranges)

    print(bam.gr_f)

    bam.gr_qwidth=qwidth(bam.gr_f)
    l=length(C_loc.list)

    CnCG_ind_rle.list=lapply(seq(1,l), function(x) convert_patLoc2rle(x, bam.gr_qwidth, C_loc.list, CG_loc.list ))

    rm(C_loc.list)
    rm(CG_loc.list)

    gc()

    print("==varify C, GC and MM base overlaps==")

    #### varify C, GC and MM base overlaps ##
    #####
    #x=sample( seq(1,(min(badInd)-1)), size = 1)
    if(length(badInd)>0){
      x=seq(1, length(bam.gr_f))[-badInd][10] # any random read after removing bad index
    }else {
      x=seq(1, length(bam.gr_f))[10]
    }
    print(bam.gr_f[x])

    #sum(as.vector(MM_tag_list_rle.filt[[x]]))
    #which(as.vector(MM_tag_list_rle.filt[[x]])==1)
    #which(as.vector(CnCG_ind_rle.list[[x]]$c.ind_rle)==1)
    #aln.methy_first$seq[x]
    #aln.methy_first$tag$MM[x]
    #aln.methy_first$tag$ML[x]

    se=(bam.gr_f@elementMetadata@listData$seq[x])
    #if(as.vector(strand(bam.gr_f[x]))=='-') {se=reverseComplement(se) }
    p=vmatchPattern('C',se)[[1]]@start[which(as.vector(MM_tag_list_rle.filt[[x]])==1)] ## all modifications at C of the read
    MMs.in.Cs= (!any(is.na(p)) & length(which(as.vector(MM_tag_list_rle.filt[[x]])==1))==length(p)) #  modified C and MMtags counts are equal
    #print(" modified C and MMtags counts are equal")
    print(paste0(">> modified C and MMtags counts are equal ? : ",MMs.in.Cs))

    cpg.loc=vmatchPattern('CG',se)[[1]]@start  # CpG locations on the read
    #intersect(cpg.loc, p)

    #print("all MM overlap with CpGs")
    MMs.in.CpGs= (length(intersect(cpg.loc, p))== length(which(as.vector(MM_tag_list_rle.filt[[x]])==1))) # all MM overlap with CpGs
    mod_cpgs=(length(p)/length(cpg.loc))*100
    print(paste0(">> all MM overlap with CpGs: ", MMs.in.CpGs, ' ovl is= ', mod_cpgs, ' %' ))
    #print("CpGs >=  modifications")
    CpGs.GrEq.MMs= (length(cpg.loc)>= length(which(as.vector(MM_tag_list_rle.filt[[x]])==1)))
    print(paste0(">> CpGs >=  modifications : ",CpGs.GrEq.MMs))


    #####

    #CnCG_ind_rle.list[[1]]

    bam.gr_f@elementMetadata@listData$CnCG_ind_rle.list = CnCG_ind_rle.list
    if(length(badInd)>0){
      bam.gr2=bam.gr_f[-badInd]
      rm(badInd)
    }else{
      bam.gr2=bam.gr_f
    }

    rm(CnCG_ind_rle.list); rm(bam.gr_f)
    gc()

    bam.gr2@elementMetadata@listData$MM_rle = MM_tag_list_rle.filt
    bam.gr2@elementMetadata@listData$ML = ML_tag_list.filt

    rm(MM_tag_list_rle.filt)
    rm(ML_tag_list.filt)

    #print(bam.gr2)

    #return(bam.gr2)

    #saveRDS(bam.gr2, paste0(OUTDIR,tools::file_path_sans_ext(basename(BAMFILE)),'_up.rds') )
    #rdsFile=paste0(OUTDIR,'/',tools::file_path_sans_ext(basename(BAMFILE)),'_chunk',chunkInd_f,'.rds')
    rdsFile=paste0(OUTDIR,'/',tools::file_path_sans_ext(basename(chunk_file)),'.rds')
    saveRDS(bam.gr2, rdsFile )

    #con <- gzfile(rdsFile, "wb", compression = 1)
    #saveRDS(bam.gr2, con)
    #close(con)

    print(paste0("===Chunk ",chunk_file," Completed=="))
    print("")  # Adds a blank line

    return(TRUE)

}

### Main Code ##

f <- gtools::mixedsort(Sys.glob(paste0(BAMDIR,paste0('/sorted.chunk*.bam'))))
print('== CHUNKS FILES ==')
print(f)
param <- ScanBamParam(what=c("qname","flag","rname", "strand", "pos", "qwidth", "mapq", "cigar","mrnm","mpos","isize","seq","qual"),tag=c("ML","MM"))

#runOut_success=mclapply(seq(1,length(f)), function(x) operation_on_chunk(f[x]), mc.cores=NCPUS, mc.preschedule=TRUE, mc.cleanup=TRUE )
runOut_success=lapply(seq(1,length(f)), function(x) operation_on_chunk(f[x]) )

## Varify RUN is succces or failure ##
false_indices <- which(!unlist(runOut_success))
if (length(false_indices) > 0) {
  print(f[false_indices])
  stop(paste("Error: FALSE element(s) found at index/indices:", paste(false_indices, collapse = ", ")))
} else {
  print(">> All RUNS are successful")
}
print("===Task Completed==")

######





