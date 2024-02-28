

#-----Intergrated transcriptome,Biotag and Clicktag
library(Seurat)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)


#load data
smk <- read_tsv("tag16_umi_tag.tsv",col_names = FALSE)
WTA <- Read10X("CINTERseq_CT26hHER2")
tag.smp <- read_tsv("tag5-7.txt", col_names = T)
projectname <- c("tag16")
biotin <- read_tsv("tagbiotin_umi_tag.tsv")


biotin <- biotin[,-3]
biotin <- t(biotin)
colnames(biotin) <- biotin[1,]
biotin <- biotin[-1,]
biotin <- t(biotin)
rownames(biotin) <- "biotin"

biotin.scal <- as.numeric(biotin)
biotin.scal <- log10(biotin.scal+10)%>%t()
colnames(biotin.scal) <- colnames(biotin)
rownames(biotin.scal) <- "biotin"
biotin.scal <- rbind(colnames(biotin.scal), biotin.scal)

WTA <- as.matrix(WTA)
WTA <- rbind(colnames(WTA), WTA)
he <- inner_join( as.data.frame(t(biotin.scal)),as.data.frame(t(WTA)), by="V1")%>%t()
colnames(he) <- he[1,]
he <- he[-1,]
he <- he[,order(he[1,],decreasing = T)]
he <- he[-1,]
mix <- CreateSeuratObject(he, project = projectname)

biotin.scal <- biotin.scal[-1,]
biotin.scal <- t(biotin.scal)
rownames(biotin.scal) <- "biotin"
mix[['biotin']] <- CreateAssayObject(counts = biotin.scal)

tag <- strsplit(tag.smp$tag, split = "_")%>%unlist()%>%table()%>%names()
mu <- c("Multiplet","Multiplet","Undetermined","Undetermined")%>%
  matrix(nrow = 2,ncol = 2,byrow = T)
colnames(mu) <- colnames(tag.smp)
tag.smp <- rbind(mu,tag.smp)


smk[1,1] <- "cell"
cell_ln <- smk[1,]
colnames(smk) <- cell_ln
smk <- smk[-1,]
smk <- t(smk)
cell_ln <- smk[1,]
colnames(smk) <- cell_ln
smk <- smk[-1,]

smk.used <- matrix(ncol = ncol(smk) )
for (i in 1:length(tag)) {
  smk.used <- rbind(smk.used, smk[tag[i],])
}
smk.used <- smk.used[-1,]
rownames(smk.used) <- tag


joint.bcs <- intersect(colnames(he), colnames(smk))


smk.tag <- smk[nrow(smk), joint.bcs]
smk.data <- smk.used[, joint.bcs]



mix[['HTO']] <- CreateAssayObject(counts = smk.data)

mix@meta.data$tag <- smk.tag

for (i in 1:nrow(tag.smp)) {
  mix@meta.data[which(mix@meta.data$tag == tag.smp$tag[i]), "sample"] <- tag.smp$sample[i]
}


DefaultAssay(mix) <- "HTO"
mix <- NormalizeData(mix, normalization.method = 'CLR')
#RidgePlot(mix, assay = "HTO", features = "SMK7")
mix <- ScaleData(mix, features = rownames(mix), verbose = FALSE)
mix <- FindVariableFeatures(mix, selection.method = "vst", nfeatures = length(tag))
mix <- RunPCA(mix, features = rownames(mix), approx = FALSE)


#SMK tsne
mix <- RunTSNE(mix, dims = 1:5, perplexity = 100, check_duplicates = FALSE)
DimPlot(mix, label = T, group.by = "sample", reduction = 'tsne')+
  ggsave("SMK tag tsne.jpg", width = 15, height = 10)+
  ggsave("SMK tag tsne.eps", width = 15, height = 10)
mix <- RunUMAP(mix, dims = 1:length(tag))
DimPlot(mix, reduction = "umap", label = T,group.by = "sample")+
  ggsave("SMK tag umap.jpg", width = 15, height = 10)+
  ggsave("SMK tag umap.pdf", width = 15, height = 10)
#SMK heatmap
DoHeatmap(mix, features = tag, group.by = "sample", assay = 'HTO')+
  ggsave('smk heatmap.jpg', width = 15, height = 10)+
  ggsave('smk heatmap.pdf', width = 15, height = 10)


haha <- mix
Idents(haha) <- haha@meta.data$sample
shun_xu <- tag.smp$sample
Idents(haha) <- factor(Idents(haha), levels = shun_xu)
VlnPlot(haha, 'biotin', pt.size = 0 )

VlnPlot(haha, 'biotin', pt.size = 0 )+geom_boxplot(width=.2, col='black', fill="white")+NoLegend()

#Remove Multiplet和undetermined
haha <- subset(haha, sample=='Multiplet'|sample=='Undetermined',invert=T)
VlnPlot(haha,features = "biotin",pt.size=0)+geom_boxplot(width=.2, col='black', fill="white")+geom_hline(aes(yintercept=2.5),colour="red",linetype="dashed",size=1)

# Remove label2(label1 and label2 are biological repeats)
haha1<- subset(haha,subset=sample=="label2",invert=T)
VlnPlot(haha1,features = "biotin",pt.size=0)+geom_boxplot(width=.2, col='black', fill="white")+geom_hline(aes(yintercept=2.5),colour="red",linetype="dashed",size=1)
ggsave(filename = "biotin expression.eps")


#Save rds
DefaultAssay(haha1) <- "RNA"
saveRDS(object = haha1,file = "biotin-sample.rds")

DefaultAssay(haha) <- "RNA"
saveRDS(object = haha,file = "biotin-sample-total.rds")
rm(list=ls())
gc()


#---------Clustering

library(Seurat)
library(tidyverse)
library(corrplot)
library(dplyr)
library(grid)
library(cowplot)
PROX1 <- readRDS("biotin-sample-total.rds")

PROX1[["percent.mt"]] <- PercentageFeatureSet(PROX1, pattern = "^mt-") 

VlnPlot(PROX1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0,ncol = 3)
ggsave(filename = "QC1.png",width = 8, height = 5 ) 

PROXA2 <- subset(PROX1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt<5)

set.seed(1234)
PROXA2 <- NormalizeData(PROXA2, normalization.method = "LogNormalize")
all.genes <- rownames(PROXA2)
PROXA2 <- ScaleData(PROXA2, features = all.genes)

PROXA2 <- FindVariableFeatures(PROXA2, selection.method = "vst", nfeatures = 2000)
PROXA2 <- RunPCA(PROXA2, features = VariableFeatures(object = PROXA2), verbose = FALSE)

ElbowPlot(PROXA2)
PROXA2 <- FindNeighbors(PROXA2, dims = 1:15)
PROXA3 <- FindClusters(PROXA2, resolution = 1.2)

#Run tsNE
PROXA3 <- RunTSNE(PROXA3, dims = 1:15)
DimPlot(PROXA3, reduction = "tsne", label = T )
ggsave(filename = "all_tsne.png",width = 8, height = 5 )

#Run umap
set.seed(1234)
PROXA3 <- RunUMAP(PROXA3, dims = 1:15)
DimPlot(PROXA3, reduction = 'umap', label = T)
ggsave(filename = "all_umap.png", width = 8, height = 5)

PROXA3.markers <- FindAllMarkers(PROXA3, 
                                 only.pos = TRUE,
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)
write.csv(PROXA3.markers, "cluster_markers.csv")


top10 <- PROXA3.markers %>% group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(PROXA3, features = top10$gene) + NoLegend()
ggsave(filename = "top10 markers.png")
write.csv(top10, "top10 gene.csv")
rm(all_tsne,all_umap, tsn_ump, PROXA3.markers)



#Annotation
cluster <- c("NK cells", "Monocytes","Monocytes","Macrophages","T cells","NK cells","T cells","Neutrophils","cDCs","pDCs","Neutrophils","NK cells","MSCs","NK cells","Macrophages","T cells","Mast cells","NK cells","cDCs","T cells","unknown","T cells","unknown")
names(cluster) <- levels(PROXA3)
PROXA3 <- RenameIdents(PROXA3, cluster)
PROXA3[["cell.types"]]<- Idents(object = PROXA3)
PROXA3<-subset(PROXA3,idents = 'unknown', invert = TRUE) 
DimPlot(PROXA3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "rename umap.png")
DimPlot(PROXA3, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "sample") + NoLegend()
saveRDS(PROXA3, 'rename umap all.rds')

PROXA3 <- readRDS("rename umap all.rds")

#demultiplexed umqp plot
iso_label1 <- subset(PROXA3,subset=sample=="label2",invert = TRUE)
DimPlot(iso_label1, reduction = "umap", pt.size = 0.5, group.by="sample")
ggsave(filename = "iso_label1 umap.png")
ggsave(filename = "iso_label1 umap.eps")      

DimPlot(iso_label1, reduction = "umap",  pt.size = 0.5, split.by = "sample")
ggsave(filename = "iso_label1 split.png")
ggsave(filename = "iso_label1 split.eps") 


#Biotag signal of anti-HER2-PPa group (label1)
label1 <- subset(PROXA3,subset=sample=="label1")
saveRDS(label1, 'label1.rds')


FeaturePlot(label1, "biotin", cols = c("#7fa4b8","#e5ecc5", "#c03932"), 
                         max.cutoff = 3, 
                         pt.size = 0.2,
                         min.cutoff = 1.5)

#Split violin plot
VlnPlot(iso_label1, "biotin", group.by = "cell.types", split.by = "sample", split.plot = T)+geom_boxplot(width=.2, col='black', outlier.shape = NA)+NoLegend()




