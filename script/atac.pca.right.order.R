args<-commandArgs(T)
input<-args[1]
input2<-args[2]

library(factoextra)
library(matrixStats)
data<-readRDS(input)
#data3<-data[a,]
data4<-scale(t(data),center = TRUE, scale = F)

samples<-read.table(input2,header=T,sep="\t")
sample_2<-data.frame(duplicates=samples[,2],comparison=samples[,3])
rownames(sample_2)<-samples[,1]
#sample_2<-unique(sample_2)
w<-NULL
for ( i in colnames(data)){
    a<-which(rownames(sample_2) == i)
    w<-c(w,a)
}
sample_3<-sample_2[w,]

res.pca <- prcomp(data4, scale = F)
pdf("atac.pca.pdf")
fviz_pca_ind(res.pca,
             col.ind = sample_3$comparison, # Color by the quality of representation
#             geom = "point",
#             col.ind.sup = c("#8dd3c7"),
#             fill.ind= c("#8dd3c7", "#fb8072", "#fdb462"),
#             palette = c("#8dd3c7", "#fb8072", "#fdb462"),
             repel = TRUE     # Avoid text overlapping
)+theme_classic()+theme(text=element_text(size=20,,family="Helvetica")) + theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

sampleDists <- dist( t( data ) )
sampleDistMatrix <- as.matrix( sampleDists )
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
par(oma=c(10,4,4,2))
pdf("distance.samplematrix.atac.pdf")
heatmap.2( sampleDistMatrix, trace="none", col=colours,margins=c(14,14))
dev.off()
