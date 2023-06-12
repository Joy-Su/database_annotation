### ARGS: [enrichment result] [output dir] [title]

#==========================================================
# the input enrichment result like :
# COG	number	COG_description
# A	45	A: RNA processing and modification
#==========================================================

args<-commandArgs(TRUE)
input=args[1]
dir=args[2]

library(ggplot2)
library(reshape2)
df<-read.table(input,head=T,sep="\t",fill=TRUE)

# bar plot
#labx = paste(title," Classification",sep="")
p <- ggplot(data=df, aes(x=eggNOG, y=number, fill=eggNOG_description)) +
geom_bar(stat="identity",width=0.7,position=position_dodge()) +
xlab("eggNOG class") + ylab("Number of genes") + 
labs(title="eggNOG Function Classification")

p<-p+theme(panel.border=element_rect(fill=NA,colour="black"))
p <- p + theme(panel.background = element_rect(fill = "transparent",colour =NA),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.background = element_rect(fill = "transparent",colour =NA))

ggsave(filename=paste(dir,".barplot.pdf",sep=""), plot=p, height=7, width=15)
ggsave(filename=paste(dir,".barplot.png",sep=""), type="cairo-png", plot=p, height=7, width=15)
