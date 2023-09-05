args<-commandArgs(TRUE)

suppressWarnings(library(ggplot2))
suppressWarnings(library(RColorBrewer))

input <- args[1]
output <- args[2]

input_file <- input
output_file<- output
data<-read.table(input_file,header=FALSE)
data2<-data.frame(region=data[,1],num=rep(1,dim(data)[1]))
df2<-aggregate(num~region,data2,sum)
df3<-data.frame(region=df2$region,cluster=rep(input,length(dim(data)[1])),num=df2$num/sum(df2$num)*100)
write.table(df3,file=output_file,quote=FALSE,row.names = F,sep="\t")
