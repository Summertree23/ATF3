library("pheatmap")
suppressMessages(library( "gplots" ))
library( "RColorBrewer" )
args<-commandArgs(T)
input<-args[1]
input2<-args[2]
output<-args[3]

data<-readRDS("vsd_matrix.rds")
comparison<-read.table("comparison.xls",header=T,sep="\t")
#comparison$new<-paste(comparison$comparison_time,comparison$Sample,sep="_")
group<-NULL
for( i in unlist(strsplit(input2,"--"))){
    a<-grep(i,comparison$comparison)
    group<-c(group,a)
}
extract<-comparison[group,]$sample

peak<-read.table(input,header=T,sep="\t")
peak2 <- peak[which( peak$padj < 0.05), ]
peak2 <- peak2[which(peak2$log2FoldChange > 1 | peak2$log2FoldChange < -1),]
target<-data[as.character(rownames(peak2)),extract]

color = colorRampPalette(brewer.pal(3,"RdBu"))(21)
breaks<-seq(-1.5,1.5,by=0.15)
pdf(output)
#out<-pheatmap(target2, show_rownames=F, cluster_rows=T,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),annotation_col=mergerd_new3,clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks,annotation_colors = list(subtype=c(CL="#F8766D",MS="#00BA38",PN="#619CFF"),cluster=c("1"="#7570B3","2"="#D95F02","3"="#E6AB02")))
#cluster=c(1="#fc8d59",2="#ffffbf",3="#99d594")
pheatmap(target, show_rownames=F, cluster_rows=F,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)
pheatmap(target, show_rownames=F, cluster_rows=T,cluster_cols = F,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)
pheatmap(target, show_rownames=F, cluster_rows=T,cluster_cols = T,clustering_distance_cols="euclidean",scale="row",show_colnames= T,color = rev(color),clustering_distance_rows="euclidean", cex=1, clustering_method="complete", border_color=FALSE,fontsize_row=0.2,fontsize_col=4,breaks=breaks)
dev.off()


