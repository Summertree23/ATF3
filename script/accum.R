
library(ggplot2)
df<-read.table("accum.txt",header=T)
pdf("accum.pdf")
ggplot(df, aes(x = cluster, y = num))+scale_fill_brewer(palette="RdBu")+geom_col(aes(fill = region ), width = 0.7)+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(df, aes(x = cluster, y = num))+scale_fill_brewer(palette="Blues")+geom_col(aes(fill = region ), width = 0.7)+theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
