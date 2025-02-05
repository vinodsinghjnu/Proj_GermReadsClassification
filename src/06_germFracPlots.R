


# List files
files <- Sys.glob("/Users/au734493/mount/genomedk_mutScanning//*_chrBAMs_hg38//*_BinoClassificationResults//05_GermFrac/germFrac2.txt")
print(files)
data=lapply(files, function(x) read.table(x, header = T))

# header
samples=gsub('_chrBAMs_hg38','',stringr::str_split_i(files, '/', 7))




#read.table(files[1], header = T)

frac=as.data.frame(do.call(cbind, sapply(seq(1,length(data)), function(x) data[[x]]['germFrac'])))
reads_remained=as.data.frame(do.call(cbind, sapply(seq(1,length(data)), function(x) data[[x]]['readsRemained'])))
colnames(frac)=samples; rownames(frac)=seq(1,26)
colnames(reads_remained)=samples; rownames(reads_remained)=seq(1,26)
reads_remained=reads_remained/100
frac=frac/100


## PLOT ##

dataNames=c('aa_sperm', 'ob006','kipenzi','lenny','marlock','martin', 'mongo','sethi','tarzan')
library(ggplot2);library(patchwork);library(plotly); library('gridExtra')
# kipenzi: G_T, lenny: PTM_T , marlock:C_T, martin: C_T, mongo: C_T ,sethi: , tarzan GIB_T: 

plotFun=function(f_data, yaxisLabel, ind){
  plotData=f_data # head(frac)
  plotData$LL_diff=seq(0,nrow(plotData)-1)
  plotData_boxdf=reshape2::melt(plotData, 'LL_diff')
  plotData_boxdf$LL_diff=as.character(plotData_boxdf$LL_diff)
  plotData_boxdf$variable=factor(x=plotData_boxdf$variable, levels = dataNames  )
  #plotData_boxdf$shape=factor(c(rep(21,6), rep(18,3), rep(13,2), rep(8,2), 25 ), levels=c(21, 18, 13, 8, 25))
  plotData_boxdf$shape=plyr::mapvalues(plotData_boxdf$variable, from = dataNames, to = c(21, 18, 13, 10, rep(25,3), 11, 12 ))
  
  addline_format <- function(x,...){
    gsub('\\s','\n',x)
  }
  
  #modified_xticks=addline_format(c('Ob2v HumanTestis', 'Ob003h HumanTestis', 'Da1 HumanSperm', 'La4 HumanSperm','Stephan ChimpTestis',
  #                  'Carl ChimpTestis', 'Samson GorillaTestis', 'Mbewe GorillaTestis' ))
  
  #modified_xticks=c('HumanTestis: Ob2v', 'HumanTestis: Ob003h', 'HumanTestis: Ob006', 'HumanTestis: Ob007', 'HumanTestis: Ob008', 'HumanTestis: Ob009', 'HumanSperm: Da1 ', 'HumanSperm: La4','ChimpTestis: Stephan',
  #                  'ChimpTestis: Carl', 'GorillaTestis: Samson', 'GorillaTestis: Mbewe' )
  
  
  plotData_boxdf$LL_diff=factor( x=plotData_boxdf$LL_diff ,  levels = sort(unique(as.numeric(plotData_boxdf$LL_diff))))
  p=ggplot(data = plotData_boxdf, aes(x = LL_diff, y = value, group = variable)) +
    ylab(yaxisLabel) +
    geom_line(aes(color = variable), show.legend = (ind != 2) ) +  # Hide legend if ind is 2
    geom_point(aes(shape = shape, color = variable), show.legend = (ind != 1)) +     # Hide legend if ind is 1
    scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
    guides(color = ifelse(ind == 2, "none", "legend"),        # Hide color legend if ind is 2
           shape = ifelse(ind == 1, "none", "legend")) +
    theme_minimal()  # Optional: Add a minimal theme for better aesthetics
  
  #library(plotly)
  #p2=ggplotly(p)
  p
  
}




