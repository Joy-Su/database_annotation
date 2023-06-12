library(optparse)
library(ggplot2)
library(reshape2)
option_list<-list(
  make_option(c('-i','--input'),
              type = "character", default = NULL,
              help = '.KEGG.classification.xls'),
  make_option(c('-o','--outdir'),
              type = 'character',
              help = 'outdir path')
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
KObarplot <- function(input, outdir){
  dat <- read.table(input,header=T,sep="\t",check.names=F)
  dat$Classification_level2 <- factor(dat$Classification_level2, levels=dat$Classification_level2)
  p <- ggplot(dat, aes(x=Classification_level2,y=Gene_number,fill=Classification_level1)) + ylab("Counts of Genes") +xlab("") + labs(title="KEGG Classification")
  p <- p + geom_bar(stat="identity",position="dodge",width=0.7) +geom_text(aes(label=Gene_number),vjust=0.4,hjust=-0.2,size=2.5) +coord_flip() 
  p <- p + theme(panel.border=element_rect(fill=NA,colour="black"))
  p <- p + theme(panel.background = element_rect(fill="transparent",colour=NA),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.background = element_rect(fill="transparent",colour=NA) )
  p <- p + theme(legend.title=element_blank())+ theme(legend.key.height=unit(0.5,"cm"),legend.text=element_text(size=8))
  p <- p + theme(axis.text.y=element_text(size=8,color="black")) + theme(axis.text.x=element_text(hjust=1, size=8,color="black"))
  ggsave(paste0(outdir,sub('.xls','.pdf',input)),width=12,height=6,plot=p)
  ggsave(paste0(outdir,sub('.xls','.png',input)),type="cairo-png",width=12,height=6,plot=p)
}

                   
argv <- commandArgs(T)
input <- opt$input
outdir <- opt$outdir
KObarplot(input, outdir)

