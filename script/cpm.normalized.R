#Rscript cpm.normalized.R peak.filter.counts.xls comparison.xls

args<-commandArgs(T)
    
input1<-args[1]
input2<-args[2]

data<-read.table(input1,header=T,sep="\t")
samples<-read.table(input2,header=T,sep="\t")
sample_2<-data.frame(duplicates=samples[,2],comparison=samples[,3])
rownames(sample_2)<-samples[,1]
w<-NULL
for ( i in colnames(data)){
    a<-which(rownames(sample_2) == i)
    w<-c(w,a)
}
sample_3<-sample_2[w,]
s<-rownames(sample_2)[w]
sample_4<-data.frame(sample_3)
rownames(sample_4)<-s

suppressMessages(library(edgeR))
print ("begin")
b<-NULL
#for(i in 1:dim(data2)[1]){
#    a<-ifelse(data2[i,] > 5,1,0)
#    b<-rbind(b,a)
#}
data3<-(data > 10) 
for(i in 1:ncol(data3)){
    if(is.logical(data3[, i]) == TRUE) data3[, i] <- as.numeric(data3[, i])
    }

print ("end")
cell_104_data<-data[which(rowSums(data3) > 2),]

dat_edgeR_norm  <- cpm(cell_104_data,log = F,prior.count = 5)
exprs_with_time<-log(x=(dat_edgeR_norm + 1),base=2)


saveRDS(exprs_with_time,file="cpm_normalized.counts.rds")
write.table(exprs_with_time,file="cpm_normalized.counts.xls",quote=F,sep="\t")


