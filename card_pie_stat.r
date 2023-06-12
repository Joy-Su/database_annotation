library(ggplot2)
library(optparse)
library(reshape2)
option_list<-list(
  make_option(c('-i','--input'),
              type = "character", default = NULL,
              help = 'input file')
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

par(mar=c(2,2,2,2))
df<-read.table(opt$input,header=T,sep="\t",check.names=F)
name = as.vector(df$AROID)
pct=round(df$Num/sum(df$Num) * 100, 2)
myLabel = paste(name, "(", pct, "%)", sep = "")
p=ggplot(df, aes(x = "", y = Num, fill = AROID)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "", title = "Functional Categories") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank(), legend.position = "right") + 
  scale_fill_discrete(breaks = df$AROID, labels = myLabel) +
  theme(panel.grid=element_blank()) + 
  theme(panel.border=element_blank())
ggsave(paste0(sub('\\.[^.]+$','',opt$input),".pdf"),width=8,height=7,plot=p)
ggsave(paste0(sub('\\.[^.]+$','',opt$input),".png"),type="cairo-png",width=8,height=7,plot=p)