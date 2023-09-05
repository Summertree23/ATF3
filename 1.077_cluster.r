load("filtered_seurat.doublet.filteredh5.RData")

suppressPackageStartupMessages({
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(scran)
    library(pheatmap)
    library(RCurl)
    library(AnnotationHub)
    library(EnsDb.Hsapiens.v86)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(dplyr)
    library(AUCell)
    library(GSEABase)
    library(future)
})

plan("multiprocess", workers = 30)
plan()

###cell cycle
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv")
cell_cycle_genes <- read.csv(text = cc_file)
ah <- AnnotationHub()
ahDb <- query(ah,
              pattern = c("Homo sapiens", "EnsDb"),
              ignore.case = TRUE)
id <- ahDb %>%
        mcols() %>%
        rownames() %>%
        tail(n = 1)
edb <- ah[[id]]
annotations <- genes(edb,
                     return.type = "data.frame")
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "S") %>%
        pull("gene_name")

# Acquire the G2M phase genes
g2m_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "G2/M") %>%
        pull("gene_name")

write.table(s_genes,file="s_gene.xls",quote=F,row.names=F)

write.table(g2m_genes,file="g2m_genes.xls",quote=F,row.names=F)

getwd()



#filtered_seurat<-b

quantile(filtered_seurat$nCount_RNA)
length(filtered_seurat$nCount_RNA)

#gene expressiom
par(mar = c(4, 8, 2, 1))
C <- filtered_seurat@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

ribo.genes<- c(rownames(filtered_seurat@assays$RNA@data)[grep("^RP[1:9]",rownames(filtered_seurat@assays$RNA@data))],rownames(filtered_seurat@assays$RNA@data)[grep("^RP[L,S]",rownames(filtered_seurat@assays$RNA@data))])

removed.genes<-c(ribo.genes)



#filtered_seurat<-b
filtered_seurat<-NormalizeData(filtered_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_seurat<-FindVariableFeatures(filtered_seurat,selection.method="vst",nfeatures=10000,verbose=T)
filtered_seurat <- CellCycleScoring(filtered_seurat, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
filtered_seurat<-RunPCA(filtered_seurat,features = c(s_genes, g2m_genes))
PCAPlot(object = filtered_seurat)
filtered_seurat$CC.Difference <- filtered_seurat$S.Score - filtered_seurat$G2M.Score
filtered_seurat<-ScaleData(filtered_seurat,vars.to.regress=c("CC.Difference"))
filtered_seurat<-RunPCA(filtered_seurat,features = VariableFeatures(filtered_seurat),npcs=20,nfeatures.print=10)

ElbowPlot(filtered_seurat,reduction="pca",ndims=20)
filtered_seurat<-RunUMAP(filtered_seurat,dims=1:20)
filtered_seurat<-RunTSNE(filtered_seurat,dims=1:20)
filtered_seurat<-RunUMAP(filtered_seurat,reduction="pca",dims=1:20,n.components=2,n.neighbors=30,n.epochs=200,min.dist=0.3,leanring.rate=1,spread=1)


#build gene-expression rankings for each cell 
filtered_seurat.list.core.matrix<-as.matrix(filtered_seurat@assays$RNA@data)
cells_rankings <- AUCell_buildRankings(filtered_seurat.list.core.matrix, nCores=1, plotStats=FALSE)
#calculate enrichment for the gene signatures
geneSets <- GeneSet(c(s_genes,g2m_genes), setName="geneSet1")
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)


#Determine the cells with the given gene signatures or active gene sets
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
cells_assignment$geneSet1$aucThr$thresholds
cells_assignment$geneSet1$aucThr$selected

#ploting
geneSetName<-rownames(cells_AUC)[grep("geneSet1", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.046)
abline(v=0.046)

pdf("cell.cycling.threshold.pdf")
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.046)
abline(v=0.046)
dev.off()

AUC<-data.frame(t(as.data.frame(cells_AUC@assays@data$AUC)))
colnames(AUC)<-"AUC"
#add AU#add AUC
filtered_seurat<-AddMetaData(filtered_seurat,metadata=AUC)
filtered_seurat$S_G2M_AUC = AUC

FeaturePlot(filtered_seurat, features = 'S_G2M_AUC',pt.size = 2.5,reduction="umap")
FeaturePlot(filtered_seurat,reduction="umap",features="S_G2M_AUC",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))


pdf("umap.cell.cycling.pdf")
FeaturePlot(filtered_seurat, features = 'S_G2M_AUC',pt.size = 2.5,reduction="umap")
dev.off()

#filtered_seurat$cell_cycle = as.numeric(filtered_seurat$S_G2M_AUC > 0.04)

#table(filtered_seurat$cell_cycle,filtered_seurat$`RNA_snn_res.0.8`)

library(ggplot2)
cyclying<-data.frame(a=c("cluster0","cluster1","cluster2"),b=c(0.21*100,0,0.15*100))
ggplot(data=cyclying, aes(x=a, y=b))+
  geom_bar(stat="identity", fill="steelblue")+theme_classic()

pdf("cell.cycling.pdf")
ggplot(data=cyclying, aes(x=a, y=b))+
  geom_bar(stat="identity", fill="steelblue")+theme_classic()
dev.off()


DimPlot(filtered_seurat,reduction="pca",group.by="orig.ident")+NoAxes()+ggtitle("PCA raw_data")
DimPlot(filtered_seurat,reduction="umap",group.by="orig.ident")+NoAxes()+ggtitle("UMAP raw_data")
DimPlot(filtered_seurat,reduction="tsne",group.by="orig.ident")+NoAxes()+ggtitle("TSNE raw_data")

filtered_seurat<-FindNeighbors(filtered_seurat,dims=1:20,k.param=60,prune.SNN=1/20)
names(filtered_seurat@graphs)
filtered_seurat



for( res in c(0.1,0.2,0.3,0.4,0.6,0.8,1)){
    filtered_seurat<-FindClusters(filtered_seurat,graph.name='RNA_snn',resolution=res,algorithm=1)
    #filtered_seurat<-FindClusters(filtered_seurat,resolution=res)
}

library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col<-sample(col_vector, n)

DimPlot(filtered_seurat, reduction = "umap", group.by = "RNA_snn_res.0.1",cols = col) + ggtitle("louvain_0.1")
DimPlot(filtered_seurat,reduction="umap",group.by="RNA_snn_res.0.2",cols = col) + ggtitle("louvain_0.2")
DimPlot(filtered_seurat,reduction="umap",group.by="RNA_snn_res.0.3",cols = col) + ggtitle("louvain_0.3")
DimPlot(filtered_seurat,reduction="umap",group.by="RNA_snn_res.0.4",cols = col) + ggtitle("louvain_0.4")
DimPlot(filtered_seurat,reduction="tsne",group.by="RNA_snn_res.0.4",cols = col) + ggtitle("louvain_0.4")
DimPlot(filtered_seurat, reduction = "umap", group.by = "RNA_snn_res.0.6",cols = col) + ggtitle("louvain_0.6")
DimPlot(filtered_seurat,reduction="umap",group.by="RNA_snn_res.0.8",cols = col) + ggtitle("louvain_0.8")
DimPlot(filtered_seurat, reduction = "umap", group.by = "RNA_snn_res.1",cols = col) + ggtitle("louvain_1")

save(filtered_seurat,file='glio_RNA.edge.cluster.Rdata')

load("glio_RNA.edge.cluster.Rdata")

library(RColorBrewer)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col<-sample(col_vector, n)

DimPlot(filtered_seurat,reduction="pca",group.by="orig.ident")+NoAxes()+ggtitle("PCA raw_data")
DimPlot(filtered_seurat, reduction = "pca", group.by = "RNA_snn_res.0.8",cols = col) + ggtitle("louvain_1")

DimPlot(filtered_seurat, reduction = "umap", group.by = "RNA_snn_res.0.8",cols = col) + ggtitle("louvain_1")
DimPlot(filtered_seurat, reduction = "tsne", group.by = "RNA_snn_res.0.8",cols = col) + ggtitle("louvain_1")

pdf("edge.cluster.pdf")
DimPlot(filtered_seurat, reduction = "umap", group.by = "RNA_snn_res.0.3",cols = col) + ggtitle("louvain_0.2")
dev.off()


write.table(filtered_seurat$`RNA_snn_res.0.4`,file="cluster.edge.txt",quote=F,sep="\t")



#length(filtered_seurat$RNA_snn_res.1.5)
print ("1")
table(filtered_seurat$RNA_snn_res.1)
print ("0.2")
table(filtered_seurat$RNA_snn_res.0.2)
print ("0.3")
table(filtered_seurat$RNA_snn_res.0.3)
print ("0.4")
table(filtered_seurat$RNA_snn_res.0.4)
print ("0.6")
table(filtered_seurat$RNA_snn_res.0.6)
print ("0.8")
table(filtered_seurat$RNA_snn_res.0.8)

#compute differential expression 

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(enrichR)
    library(rafalib)
})

filtered_seurat<-FindClusters(filtered_seurat,resolution=0.8,algorithm=1)
markers_genes <- FindAllMarkers(filtered_seurat, log2FC.threshold = 0.2, test.use = "wilcox",
    min.pct = 0.1, min.diff.pct = 0.1, only.pos = TRUE ,max.cells.per.ident = 500,
    assay = "RNA")

#markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj) -> top25


write.table(markers_genes,file="differential_genes.difference.xls",quote=F,sep="\t")


VlnPlot(filtered_seurat, features = c("DDIT4","TRIB3"),group.by = "RNA_snn_res.0.8",cols=c("#e41a1c", "#377eb8","#A9A9A9"))

FeaturePlot(filtered_seurat,reduction="umap",features="DDIT4",order=T)+NoLegend()+NoAxes()+NoGrid()
FeaturePlot(filtered_seurat,reduction="umap",features="TRIB3",order=T)+NoLegend()+NoAxes()+NoGrid()


markers_genes %>% group_by(cluster) %>% top_n(-10,p_val_adj) -> top5
#filtered_seurat<-ScaleData(filtered_seurat,features=as.character(unique(top5$gene)),assay="RNA")
DoHeatmap(filtered_seurat, features = as.character(unique(top5$gene)), group.by = "RNA_snn_res.0.8",assay = "RNA")+theme(text = element_text(size = 12))

markers_genes<-markers_genes[which(markers_genes$p_val_adj < 0.01),]
markers_genes %>% group_by(cluster) %>% top_n(-100,p_val_adj) -> top
#filtered_seurat<-ScaleData(filtered_seurat,features=as.character(unique(top5$gene)),assay="RNA")
DoHeatmap(filtered_seurat, features = as.character(unique(top$gene)), group.by = "RNA_snn_res.0.8",assay = "RNA")+theme(text = element_text(size = 12))


cluster00<-markers_genes[grep("0",markers_genes$cluster),]
cluster11<-markers_genes[grep("1",markers_genes$cluster),]
markers<-rbind(cluster00,cluster11)

markers %>% group_by(cluster) %>% top_n(-50,p_val_adj) -> top5
mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
#filtered_seurat<-ScaleData(filtered_seurat,features=as.character(unique(top5$gene)),assay="RNA")
DoHeatmap(filtered_seurat, features = as.character(unique(top5$gene)), group.by = "RNA_snn_res.0.8",assay = "RNA")+scale_fill_gradientn(colours = rev(mapal))+theme(text = element_text(size = 5))

pdf("heatmap.diff.gene.pdf")
DoHeatmap(filtered_seurat, features = as.character(unique(top5$gene)), group.by = "RNA_snn_res.0.8",assay = "RNA")+scale_fill_gradientn(colours = rev(mapal))+theme(text = element_text(size = 5))
dev.off()


markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj) -> top25


#mypar(2,5,mar=c(4,6,3,1))
for (i in unique(top25$cluster)){
    barplot(sort(setNames(top25$avg_log2FC,top25$gene)[top25$cluster==i],F),
    horiz=T,las=1,main=paste0(i,"vs.rest"),border="white",yaxs = "i")
    abline(v=c(0,0.25),lty=c(1,2))
}


VlnPlot(filtered_seurat, features = c("CDK1","CDK4","CDK6","ARID1A","MKI67","CDKN1B","CDKN1C","CDKN1A","ATF3"),group.by = "RNA_snn_res.0.8",cols=c("#e41a1c", "#377eb8","#A9A9A9"))

pdf("addtional.boxplot.pdf")
VlnPlot(filtered_seurat, features = c("ARID1A","CENPF","SOX9","SOX4","ATF4"),group.by = "RNA_snn_res.0.8",cols=c("#e41a1c", "#377eb8","#A9A9A9"))
dev.off()



pdf("additional.umap.pdf")
FeaturePlot(filtered_seurat,reduction="umap",features="ARID1A",order=T,pt.size = 2.5,cols = c("grey", "red"))
FeaturePlot(filtered_seurat,reduction="umap",features="CENPF",order=T,pt.size = 2.5,cols = c("grey", "red"))
FeaturePlot(filtered_seurat,reduction="umap",features="SOX9",order=T,pt.size = 2.5,cols = c("grey", "red"))
FeaturePlot(filtered_seurat,reduction="umap",features="SOX4",order=T,pt.size = 2.5,cols = c("grey", "red"))
FeaturePlot(filtered_seurat,reduction="umap",features="ATF4",order=T,pt.size = 2.5,cols = c("grey", "red"))

dev.off()

plot<-VlnPlot(filtered_seurat, features = c("CDK1","CDK4","CDK6","ARID1A","MKI67","CDKN1B","CDKN1C","CDKN1A","ATF3"),group.by = "RNA_snn_res.0.8",cols=c("#e41a1c", "#377eb8","#A9A9A9"))

pdf("plot.marker.pdf")
plot
dev.off()

FeaturePlotSingle<- function(obj, feature, metadata_column, ...){
  all_cells<- colnames(obj)
  groups<- levels(obj@meta.data[, metadata_column])
    }
p_list<- FeaturePlotSingle(filtered_seurat, feature= "CLCN3", metadata_column = "orig.ident", pt.size = 1, order =TRUE)

FeaturePlot(filtered_seurat,reduction="umap",features="NFIL3",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
FeaturePlot(filtered_seurat,reduction="umap",features="PNRC1",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

#+ scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))


pdf("overlaped.gene.pdf")
FeaturePlot(filtered_seurat,reduction="umap",features="NFIL3",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
FeaturePlot(filtered_seurat,reduction="umap",features="PNRC1",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
VlnPlot(filtered_seurat, features = c("ATF4"),group.by = "RNA_snn_res.1",cols=c("#e41a1c", "#377eb8","#A9A9A9"))
dev.off()

pdf("overlaped.gene2.pdf")
FeaturePlot(filtered_seurat,reduction="umap",features="ATF4",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
FeaturePlot(filtered_seurat,reduction="umap",features="ATF3",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
VlnPlot(filtered_seurat, features = c("ATF4"),group.by = "RNA_snn_res.1",cols=c("#e41a1c", "#377eb8","#A9A9A9"))

dev.off()

getwd()

VlnPlot(filtered_seurat, features = c("NFIL3","PNRC1"),group.by = "RNA_snn_res.0.8",cols=c("#e41a1c", "#377eb8","#A9A9A9"))

pdf("overlaped.vlnplot.gene.pdf")
VlnPlot(filtered_seurat, features = c("NFIL3","PNRC1"),group.by = "RNA_snn_res.0.8",cols=c("#e41a1c", "#377eb8","#A9A9A9"))
dev.off()

pdf("atf3.gene.plot.pdf")
FeaturePlot(filtered_seurat,reduction="umap",features="ATF3",order=T,pt.size = 2.5,cols = c("grey", "red"))+theme_classic()
dev.off()

pdf("gene.pdf")
FeaturePlot(filtered_seurat,reduction="umap",features="CLCN3",order=T,pt.size = 2.5,cols = c("grey", "red"))
FeaturePlot(filtered_seurat,reduction="umap",features="PLOD1",order=T,pt.size = 2.5,cols = c("grey", "red"))
FeaturePlot(filtered_seurat,reduction="umap",features="MXD1",order=T,pt.size = 2.5,cols = c("grey", "red"))
FeaturePlot(filtered_seurat,reduction="umap",features="NR4A1",order=T,pt.size = 2.5,cols = c("grey", "red"))
dev.off()



#pathway enrichment 




library(enrichR)


gene_rank <- setNames( markers_genes$avg_log2FC, casefold(rownames(markers_genes),upper=T) )


library(msigdbr)


msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
msigdbgmt <- as.data.frame(msigdbgmt)

unique(msigdbgmt$gs_subcat)


msigdbgmt_subset <- msigdbgmt[msigdbgmt$gs_subcat == "CP:KEGG",]
gmt <- lapply( unique(msigdbgmt_subset$gs_name),function(x){msigdbgmt_subset [msigdbgmt_subset$gs_name == x ,"gene_symbol"]} )
names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name,"_",msigdbgmt_subset$gs_exact_source))

library(presto)

pbmc.genes <- wilcoxauc(filtered_seurat, 'RNA_snn_res.0.8')

dim(pbmc.genes)

dplyr::count(pbmc.genes, group)


library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
#library(tibble)



m_df<- msigdbr(species = "Homo sapiens", category = "H")#我们使用H

head(m_df)
dim(m_df)



fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets$HALLMARK_UNFOLDED_PROTEIN_RESPONSE
write.table(fgsea_sets$HALLMARK_UNFOLDED_PROTEIN_RESPONSE,file="unfolded_protein.xls",quote=F,sep="\t",row.names=F)
fgsea_sets$HALLMARK_CHOLESTEROL_HOMEOSTASIS
write.table(fgsea_sets$HALLMARK_CHOLESTEROL_HOMEOSTASIS,file="cholesterol.xls",quote=F,sep="\t",,row.names=F)
fgsea_sets$HALLMARK_HYPOXIA
write.table(fgsea_sets$HALLMARK_HYPOXIA,file="hypoxia.xls",quote=F,sep="\t",,row.names=F)
fgsea_sets$HALLMARK_TNFA_SIGNALING_VIA_NFKB
write.table(fgsea_sets$HALLMARK_TNFA_SIGNALING_VIA_NFKB,file="nfkb.xls",quote=F,sep="\t",row.names=F)


pbmc.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)  

cluster0.genes<- pbmc.genes %>%
  dplyr::filter(group == "0") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

#library(tibble)
#ranks<- deframe(cluster0.genes)

#head(ranks)

ranks<-as.vector(unlist(cluster0.genes$auc))
names(ranks)<-cluster0.genes$feature
head(ranks)

length(ranks)

saveRDS(ranks,file="ranks.cluster0.rds")

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)


fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) 

table.cluster0<-as.matrix(fgseaResTidy %>% filter(padj < 0.05))
write.table(table.cluster0,file="cluster0.pathway.significant.xls",quote=F,row.names=F,sep="\t")

# 显示top20信号通路
ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 10)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + theme_classic()
 # theme_minimal() ####以7.5进行绘图填色

pdf("cluster0.hallmark.pathway.pdf")
ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 10)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +theme_classic()
dev.off()

plotEnrichment(fgsea_sets[["HALLMARK_FATTY_ACID_METABOLISM"]],
               ranks) + labs(title="HALLMARK_FATTY_ACID_METABOLISM cluster0")

plotEnrichment(fgsea_sets[["HALLMARK_UNFOLDED_PROTEIN_RESPONSE"]],
               ranks) + labs(title="HALLMARK_UNFOLDED_PROTEIN_RESPONSE cluster0")

plotEnrichment(fgsea_sets[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]],
               ranks) + labs(title="HALLMARK_CHOLESTEROL_HOMEOSTASIS cluster0")

plotEnrichment(fgsea_sets[["HALLMARK_HYPOXIA"]],
               ranks) + labs(title="HALLMARK_HYPOXIA cluster0")

plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB cluster1")

plotEnrichment(fgsea_sets[["HALLMARK_HYPOXIA"]],
               ranks) + labs(title="HALLMARK_HYPOXIA cluster1")

plotEnrichment(fgsea_sets[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]],
               ranks) + labs(title="HALLMARK_CHOLESTEROL_HOMEOSTASIS cluster1")

data1<-read.table("cluster0.11.txt",header=F)
data2<-read.table("cluster0.22.txt",header=F)
data3<-read.table("cluster0.33.txt",header=F)

x<-list(OXIDATIVE=data1[,1],MYC=data2[,1],PROTEIN=data3[,1])


library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
  )

pdf("venn.plot.cluster0.pdf")
ggvenn(
  x, 
  fill_color = c("#b3e2cd", "#fdcdac", "#cbd5e8"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()



cluster1.genes<- pbmc.genes %>%
  dplyr::filter(group == "1") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)



ranks<-as.vector(unlist(cluster1.genes$auc))
names(ranks)<-cluster1.genes$feature
head(ranks)

saveRDS(ranks,file="ranks.cluster1.rds")

fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)


u<-fgseaResTidy %>% dplyr::select( -ES, -nMoreExtreme)  %>%
  arrange(desc(NES))

e<-as.matrix(as.data.frame(u))

write.table(e,file="077.pathway.output.xls",quote=F,sep="\t",row.names=F)



fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaResTidy %>%
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  arrange(padj) %>% head(40) 

table.cluster1<-as.matrix(fgseaResTidy %>% filter(padj < 0.05))
write.table(table.cluster1,file="cluster1.pathway.significant.xls",quote=F,row.names=F,sep="\t")

table.cluster1

# 显示top20信号通路
ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 10)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +theme_classic()
 # theme_minimal() ####以7.5进行绘图填色

pdf("cluster1.hallmark.pathway.pdf")
ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 5), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 10)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") +theme_classic()
dev.off()

plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB group1")

#library(enrichplot)
#gseaplot2(fgseaResTidy,geneSetID = 1, title = fgseaResTidy$pathway[1])

plotEnrichment(fgsea_sets[["HALLMARK_HYPOXIA"]],
               ranks) + labs(title="HALLMARK_HYPOXIA cluster1")

plotEnrichment(fgsea_sets[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]],
               ranks) + labs(title="HALLMARK_CHOLESTEROL_HOMEOSTASIS group1")

plotEnrichment(fgsea_sets[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]],
               ranks) + labs(title="HALLMARK_CHOLESTEROL_HOMEOSTASIS cluster0")

plotEnrichment(fgsea_sets[["HALLMARK_FATTY_ACID_METABOLISM"]],
               ranks) + labs(title="HALLMARK_FATTY_ACID_METABOLISM cluster1")

plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
               ranks) + labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB cluster1")

plotEnrichment(fgsea_sets[["HALLMARK_HYPOXIA"]],
               ranks) + labs(title="HALLMARK_HYPOXIA cluster1")

plotEnrichment(fgsea_sets[["HALLMARK_CHOLESTEROL_HOMEOSTASIS"]],
               ranks) + labs(title="HALLMARK_CHOLESTEROL_HOMEOSTASIS cluster1")





fgseaRes[which(fgseaRes$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),]

fgseaRes[which(fgseaRes$pathway == "HALLMARK_FATTY_ACID_METABOLISM"),]

fgseaRes[which(fgseaRes$pathway == "HALLMARK_CHOLESTEROL_HOMEOSTASIS"),]


as.character(fgseaRes[which(fgseaRes$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),][,8])

as.character(fgseaRes[which(fgseaRes$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),][,8])

write.table(as.character(fgseaRes[which(fgseaRes$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),][,8]),file="table1.xls")

as.character(fgseaRes[which(fgseaRes$pathway == "HALLMARK_HYPOXIA"),][,8])

as.character(fgseaRes[which(fgseaRes$pathway == "HALLMARK_CHOLESTEROL_HOMEOSTASIS"),][,8])

as.character(fgseaRes[which(fgseaRes$pathway == "HALLMARK_MTORC1_SIGNALING"),][,8])

as.character(fgseaRes[which(fgseaRes$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),][,8])

getwd()



data1<-read.table("c11.txt",header=F)
data2<-read.table("c22.txt",header=F)
data3<-read.table("c33.txt",header=F)
data4<-read.table("c44.txt",header=F)

x<-list(nfkb=data1[,1],hypoxia=data2[,1],cholesterol=data3[,1],unfold=data4[,1])


library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )

pdf("venn.plot.pdf")
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()



data1<-read.table("c11.txt",header=F)
data2<-read.table("c22.txt",header=F)
data3<-read.table("c33.txt",header=F)

x<-list(nfkb=data1[,1],hypoxia=data2[,1],cholesterol=data3[,1])


library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#fbb4ae", "#b3cde3", "#ccebc5"),
  stroke_size = 0.5, set_name_size = 4
  )

pdf("venn.plot.cluster1.pdf")
ggvenn(
  x, 
  fill_color = c("#fbb4ae", "#b3cde3", "#ccebc5"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()




c0_string<-read.table("/disk1/xilu/collaborate/xiaonan/single_RNA/077/cluster0.mp.top15.txt",sep="\t",header=T)
c1_string<-read.table("/disk1/xilu/collaborate/xiaonan/single_RNA/077/cluster1.mp.top15.txt",sep="\t",header=T)


c0_string

c0_string$group=rep("c1",6)
c1_string$group=rep("c2",6)

c0_string$term.description<-factor(c0_string[,2])
c1_string$term.description<-factor(c1_string[,2])

all_pathway<-rbind(c1_string,c0_string)
all_pathway$group<-factor(all_pathway$group,levels=c("c1","c2"))

all_pathway

library(ggplot2)

bubble_plot<-ggplot(all_pathway,aes(x=group,y=term.description,size=strength,color=false.discovery.rate	))+theme_minimal()+geom_point(alpha=1)+theme(text = element_text(size=5),axis.text.y = element_text(hjust=1))+scale_color_gradient(low="red",high="red")
bubble_plot

pdf("bubbleplot.pdf")
bubble_plot
dev.off()

getwd()
