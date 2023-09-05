#Rscript filter.counts.R peak.new.xls

args<-commandArgs(T)

input<-args[1]

data2<-read.table(input,header=T,sep="\t")
data<-data2[,2:dim(data2)[[2]]]
rownames(data)<-data2[,1]
data3<-(data > 10) 
for(i in 1:ncol(data3)){
    if(is.logical(data3[, i]) == TRUE) data3[, i] <- as.numeric(data3[, i])
    }

print ("end")
cell_104_data<-data[which(rowSums(data3) > 2),]

saveRDS(cell_104_data,file="peak.filter.counts.rds")
write.table(cell_104_data,file="peak.filter.counts.xls",quote=F,sep="\t")


