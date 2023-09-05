args<-commandArgs(T)
input<-args[1]
input2<-args[2]
input3<-args[3]

comparisons<-unlist(strsplit(input3,"--"))

data<-readRDS(input)
samples<-read.table(input2,header=T,sep="\t")
sample_2<-data.frame(duplicates=samples[,2],comparison=samples[,3])
rownames(sample_2)<-samples[,1]
#sample_2<-unique(sample_2)
w<-NULL
for ( i in colnames(data)){
    a<-which(rownames(sample_2) == i)
    w<-c(w,a)
}
sample_4<-sample_2[w,]
sample_3<-data.frame(group=sample_4[,2])
rownames(sample_3)<-rownames(sample_4)

library(DESeq2)
atacDDS <- DESeqDataSetFromMatrix(round(data,digits=0), sample_3, ~ group)
#atacDDS <- collapseReplicates ( atacDDS2, groupby = atacDDS2$collapse)
#atacDDS <- DESeqDataSetFromMatrix(round(data2,digits=0), merged_new2, ~ cluster )
#atacDDS_deseq <- DESeq(atacDDS,fitType="local")
#atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = 11,parallel=T)
atacDDS_deseq <- DESeq(atacDDS,minReplicatesForReplace = Inf,parallel=T)
vsd <- vst(atacDDS_deseq, blind=FALSE)
matrix<-assay(vsd,blind=FALSE)
saveRDS(matrix,file="vsd_matrix.rds")

#atac_Rlog <- rlog(atacDDS_deseq,fitType = "local")
result1 <- results(atacDDS_deseq, contrast=c('group',comparisons[1],comparisons[2]))

#filter 20 
#resSig <- result1[which(result1$padj < 0.05 & result1$pvalue < 0.05), ]
resSig <- result1[which( result1$padj < 0.05), ]
resSig2 <- resSig[which(resSig$log2FoldChange > 1),]
resSig2 <- resSig2[which(resSig2$baseMean > 10 ),]
file_name<-paste(comparisons[1],comparisons[2],"high.xls",sep="_")
file_all<-paste(comparisons[1],comparisons[2],"all.xls",sep="_")
write.table(resSig2,file=file_name,quote=F,sep="\t")
write.table(result1,file=file_all,quote=F,sep="\t")
