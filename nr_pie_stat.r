library(ggplot2)
library(reshape2)
library(optparse)
option_list<-list(
  make_option(c('-i','--input'),
              type = "character", default = NULL,
              help = 'input file'),
  make_option(c('-m','--limit_max'),
              type = 'integer', default = 10,
              help = 'limit max')
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

nr_pie_stat<-function(input, limit_max=10){
  par(mar=c(2,2,2,2))
  
  df<-read.table(input,header=T,sep="\t",check.names=F)
  name = as.vector(df$Species_Name)
  pct=round(df$Gene_Number/sum(df$Gene_Number) * 100, 2)
  myLabel = paste(name, "(", pct, "%)", sep = '')
  p=ggplot(df, aes(x = "", y = Gene_Number, fill = Species_Name)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    theme_bw() +
    theme(axis.ticks = element_blank()) +
    theme(axis.text.x = element_blank()) +
    labs(x = "", y = "", title = paste0("Top ",limit_max," species distribution")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank(), legend.position = "right") + 
    scale_fill_discrete(breaks = df$Species_Name, labels = myLabel) +
    theme(panel.grid=element_blank()) + 
    theme(panel.border=element_blank())
  ggsave(sub('.xls','.pdf',basename(input)),width=8,height=7,plot=p)
  ggsave(sub('.xls','.png',basename(input)),type="cairo-png",width=8,height=7,plot=p)
  
}

argv <- commandArgs(T)
input <- opt$input
limit_max<-opt$limit_max
nr_pie_stat(input, limit_max)
