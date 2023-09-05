suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
library("Hmisc")
library(reshape2)
library(ggplot2)
library( "RColorBrewer" )
library(pheatmap)

args<-commandArgs(T)
input<-args[1]
output<-args[2]

expr<-read.table(input,header=T,sep="\t")
res2<-rcorr(as.matrix(expr),type = "pearson")
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

# Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <-cormat[hc$order, hc$order]
}

cormat <- res2$r
min(cormat)

breaks<-seq(0.5,1,by=0.025)
colour<-colorRampPalette(brewer.pal(7,"RdYlBu"))(21)
pdf(output)
pheatmap(cormat,fontsize_row=4,display_numbers =round(cormat,2),show_colnames= F,treeheight_col = 0,show_rownames=T,breaks=breaks,color=rev(colour),border_color = "NA")
dev.off()