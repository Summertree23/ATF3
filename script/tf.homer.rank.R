library(dplyr)
library(ggplot2)
library(ggrepel)

args<-commandArgs(TRUE)
input<-args[1]
output<-args[2]
cluster1<-read.table(input,header=T,sep="\t")
#cluster1_1<-data.frame(Motif_name=cluster1$Motif_name,Log_P.value=-cluster1$Log_P.value,percent=cluster1$of.Target.Sequences.with.Motif)
#cluster1_1<-data.frame(Motif_name=cluster1$Motif.Name,Log_P.value= -cluster1$Log.P.value,percent=cluster1$X..of.Target.Sequences.with.Motif)
cluster1_1<-data.frame(Motif_name=cluster1$Motif.Name,Log_P.value= -cluster1$Log.P.value,fdr=cluster1$q.value..Benjamini.,percent=cluster1$of.Target.Sequences.with.Motif,rank=seq(1,length(cluster1$Motif.Name),1))
#sig<-ifelse(cluster1_1$percent>31.70,cluster1_1$Log_P.value > 2,"Not sig")
#Label<-ifelse(cluster1_1$percent>29  & cluster1_1$fdr < 0.05,"Top10 Enrichment","Rest")
Label<-ifelse(cluster1_1$rank < 11,"Top10 significant","Rest")
pdf(output)
#ggplot(cluster1_1, aes(x = percent, y = Log_P.value)) +geom_point()+geom_text(aes(label=ifelse(percent>31.70,as.character(Motif_name),'')),hjust=0,vjust=0)
#ggplot(cluster1_1, aes(x = percent, y = Log_P.value)) +geom_point(color = dplyr::case_when(cluster1_1$percent > 31.70 ~ "#1b9e77",cluster1_1$percent < 31.70 ~ "#d95f02"))+geom_text(aes(label=ifelse(percent>31.70,as.character(Motif_name),'')),hjust=0,vjust=0)
#ggplot(cluster1_1, aes(x = percent, y = Log_P.value)) +geom_point(aes(col=sig))+scale_color_manual(values=c("#1b9e77", "#d95f02"))+geom_text_repel(data=filter(cluster1_1, percent>31.70  & fdr < 0.05), aes(label=as.character(Motif_name)))+ theme_classic()
#ggplot(cluster1_1, aes(x = rank, y = Log_P.value)) +geom_point(aes(col=Label))+scale_color_manual(values=c("#000000", "#d94602"))+geom_text_repel(data=filter(cluster1_1, rank < 11), aes(label=as.character(Motif_name)),colour="#d94602")+ylab("-log Pvalue")+xlab("Enrichment percentage")+geom_hline(yintercept=4.5, linetype="dashed", color = "black")+geom_vline(xintercept=29,linetype="dashed",color="black")+theme_classic()
ggplot(cluster1_1, aes(x = rank, y = Log_P.value)) +geom_point(aes(col=Label))+scale_color_manual(values=c("#000000", "#d94602"))+geom_text_repel(data=filter(cluster1_1, rank < 11), aes(label=as.character(Motif_name)),colour="#d94602")+ylab("-log Pvalue")+xlab("rank")+theme_classic()
dev.off()

