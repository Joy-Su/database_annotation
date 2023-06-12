library(ggplot2)
library(optparse)
option_list<-list(
  make_option(c('-i','--input'),
              type = "character", default = NULL,
              help = 'input file'),
  make_option(c('-o','--outdir'),
              type = 'character', default = NULL,
              help = 'output directory'),
  make_option(c('-n','--name'),
              type = 'character', default = 'Unigene',
              help = 'prefix')
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

count <- read.table(file=opt$input, header=T, quote="", sep="\t", check.names=F)
count <-count[,1:3]
p=ggplot(data=count, aes(x=Class_Definition, y=Genes_Count,fill=Class_Definition))+
  geom_bar(stat="identity", width = 0.5, position = position_dodge(0.7))+
  guides(fill=FALSE)+
  theme(panel.background=element_rect(fill='transparent',color ="gray"),  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,color = "black",size=9))+
  theme(panel.grid =element_blank())+
  geom_text(mapping = aes(label = count$Genes_Count),size=3.5,vjust=-0.2)+
  xlab("Class definition") + ylab("Genes number")
ggsave(paste0(opt$outdir,'/',opt$name,".CAZy.class.pdf"),width=8,height=7,plot=p)
ggsave(paste0(opt$outdir,'/',opt$name,".CAZy.class.png"),type="cairo-png",width=8,height=7,plot=p)
