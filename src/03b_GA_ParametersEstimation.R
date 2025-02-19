library("optparse")

option_list = list(
  make_option(c("-I", "--INPUTFILE"), type="character", default=NULL,
              help="Methylation Frequency Table", metavar="character"),
  make_option(c("-D", "--DATANAME"), type="character", default=NULL,
              help="name of the chromosome in NCBI format", metavar="character"),
  make_option(c("-O", "--OUTFILE"), type="character", default=NULL,
              help="OUT data directory", metavar="character"),
  make_option(c("-N", "--NCPUS"), type="character", default=NULL,
              help="Number of CPUS for GA", metavar="character")
);

## Libraries ##
library(EnvStats);library(tibble);library(patchwork); library(GA); library(ggplot2); library(parallel)



opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt)<1){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


cat(">>======\n\n")

# description ##

INPUTFILE=opt$INPUTFILE
DATANAME=opt$DATANAME
OUTFILE=opt$OUTFILE
NCPUS=opt$NCPUS


print(paste0('INPUTFILE: INPUT Methylation table: ', INPUTFILE))
print(paste0('DATANAME: Name of the sample: ', DATANAME))
print(paste0('OUTFILE: GA parameters file: ', OUTFILE))
print(paste0('NCPUS: Number of CPUS: ', NCPUS))


args <- commandArgs()
print(paste(c('command:', args), collapse=' '))
cat("<<======\n\n")

#INPUTFILE='/home/vinodsingh/mutationalscanning/Workspaces/vinod/ob2v_chrBAMs_hg38/BinoClassificationResults/priorData_ob2v_MathylationFreq_CpG/ob2v_HiFiMethyFreqTable.tsv'
#INPUTFILE='/Users/au734493/mount/genomedk_mutScanning/ob2v_chrBAMs_hg38/BinoClassificationResults/priorData_ob2v_MathylationFreq_CpG/ob2v_HiFiMethyFreqTable.tsv'
#DATANAME='ob2v'


allDatasets_HiFiMethyFreqTable_df=read.table(INPUTFILE,header = T,  sep = '\t')
rownames(allDatasets_HiFiMethyFreqTable_df) <- allDatasets_HiFiMethyFreqTable_df[[1]]
print(head(allDatasets_HiFiMethyFreqTable_df))


l=sample(x=seq(1,256,1)/257,size = 1000000, replace = T, prob = as.vector(allDatasets_HiFiMethyFreqTable_df[,'Freq'])[seq(1,256)]/sum(allDatasets_HiFiMethyFreqTable_df[,'Freq']))
dist_l=table(l)/sum(l)
#dist_l=dist_l[sort(as.numeric(names(dist_l)), decreasing = F)]
#names(dist_l)=seq(0,255)
para.fit=ebeta(l, method = "mle")
l_beta=dbeta(x=seq(1,256,1)/257, para.fit$parameters['shape1'],para.fit$parameters['shape2']) # fit beta distribution on the methylation data
x_ticks=seq(1,256,1)/257

# GA objective function#
fit_fun <- function(s1,s2,lamda) {
  beta1 <- dbeta(x_ticks, shape1 = s2, shape2 = s1) # unmethy
  beta2 <- dbeta(x_ticks, shape1 = s1, shape2 = s2) # methy
  err=sum((l_beta-((1-lamda)*beta1+(lamda*beta2)))^2) # sum of least square errors
  err
}


GA <- ga(type = "real-valued",
         fitness = function(x) -fit_fun(s1=x[1], s2=x[2], lamda=x[3]),
         lower = c(1,.001,0), upper = c(10, 1, 1), # methylated parameters
         #fitness = function(x) -fit_fun2(s1=x[1],  lamda=x[2]), # methylated parameters
         #lower = c(1,0), upper = c(50, 1),
         #monitor = if(interactive()) gaMonitor else FALSE,
         popSize = 50, pmutation=0.03, maxiter = 10000, run=500,
         keepBest = TRUE, seed = 123, parallel = FALSE)
summary(GA)
#plot(GA)

my.m <- dbeta(x_ticks, shape1 =GA@solution[1], shape2 = GA@solution[2])
my.um <- dbeta(x_ticks, shape1 = GA@solution[2], shape2 = GA@solution[1])

#MIXTURE of methylated and unmethylated componets
my.mix_um_m <- ((1-GA@solution[3])*my.um) + (GA@solution[3]*my.m)
my.er=sum((l_beta-my.mix_um_m)^2) # LSE

## Write Parameter file #
my_paraf=cbind(data.frame(GA@solution[1], GA@solution[2], GA@solution[3]),my.er)
row.names(my_paraf)=DATANAME
colnames(my_paraf)=c('ms1','ms2','mWeight','err')
write.table(my_paraf, file = OUTFILE, quote = F, sep = '\t', row.names = TRUE)
saveRDS(my_paraf, file = gsub(".tsv$", ".rds", OUTFILE))
#my_paraf2=read.table(file = 'test.tsv', sep = '\t')


# plot data #
my.plot_df=tibble(x=x_ticks,obs=l_beta, um=my.um, m=my.m, est=my.mix_um_m )
my.plot_df.melt=reshape2::melt(my.plot_df, measure.vars = c('obs', 'um', 'm','est') )
my.plot_df.melt$variable=factor(my.plot_df.melt$variable, levels = c('um','m','obs','est' ))

p1=ggplot(my.plot_df.melt, aes(x = x, y = value, color =variable)) +
  geom_line()+
  scale_color_manual(values=c('blue','red','black','grey'))+
  ggtitle(paste0(DATANAME,' GA:mW=',signif(GA@solution[3],3), ' , ms1:', signif(GA@solution[1],3), ' , ms2:', signif(GA@solution[2],3), ' , er:', signif(my.er,3)  ) )+
  theme_minimal()

p1$data$variable <- factor(p1$data$variable, levels = c('um','m',"obs", "est"), labels = c('Unmethylated','Methylated', "Observed", "Estimated"))
p1=p1+ xlab('Methylation value') + ylab('Probablity Density Function \n (beta distribution)')+ theme(legend.position = c(0.5, 0.6))+theme(legend.title=element_blank())


ggsave(gsub(".tsv$", "_plot.png", OUTFILE),p1, width=6, height=4.5, units='in')
