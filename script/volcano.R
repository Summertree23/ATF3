args<-commandArgs(T)
input<-args[1]
output<-args[2]
library(ggplot2)
mydata<-read.table(input,header=T,sep="\t")
mydata$log2FoldChange = -(mydata$log2FoldChange)
library(dplyr)
mydata<-mydata%>%mutate(threshold = ifelse(log2FoldChange>= 1 & padj <0.05 ,"NDST", ifelse(log2FoldChange<=-1 & padj<0.05 , "control", "non-significant")))
mydata$log<- -log10(mydata$padj)
pdf(output)
g<-ggplot(mydata, aes(x=log2FoldChange,y=log))+geom_point(aes(colour = threshold), size=1) +scale_colour_manual(values = c("control"="#5ac0d0", "NDST"="#f2938e","non-significant"= "black"))
g+geom_vline(xintercept=c(-1,1),linetype="dashed",color="black")+geom_hline(yintercept=1.3, linetype="dashed", color = "black")+theme_classic()
print (g)
dev.off()
