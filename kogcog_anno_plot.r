#!/usr/local/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript bin/kog_anno_plot.r in.stat prefix")
	print("1) in.stat: the stat file for KOG anno")
	print("2) prefix.png/pdf: the filename for output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 ) {
	print(args)
	usage()
	stop("the length of args != 2")
}

# load library
library(ggplot2)
library(grid)

# get args
args <-commandArgs(TRUE)

# reading data
data <- read.delim(args[1], header=TRUE, sep="\t")

# plot
df <- data.frame(Frequency=data[,3], group=data[,1])
labels <- paste(data[,1], data[,2], sep=": ")
p <- ggplot(data=df, aes(x=group, y=Frequency)) + geom_bar(aes(fill=group), stat="identity")
p <- p + scale_fill_discrete(name="", breaks=sort(data[,1]), labels=sort(labels))
p <- p + theme(legend.key.size=unit(0.5, "cm")) 
p <- p + labs(x="Function Class", y="Gene number", title="KOG Function Classification")
p <- p + theme(panel.border=element_rect(fill=NA,colour="black"))
p <- p + theme(panel.background = element_rect(fill = "transparent",colour =NA),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),plot.background = element_rect(fill = "transparent",colour =NA))

# output plot
png(filename=paste(args[2],".png",sep=""), height = 2500, width = 5000, res = 300, units = "px")
print(p)
dev.off()

pdf(paste(args[2],".pdf",sep=""), height = 7, width = 15)
print(p)
dev.off()
