
#-----Intergrated transcriptome,Biotag and Clicktag-------#
library(Seurat)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)

#load data
smk <- read_tsv("sampletag_umi_tag.tsv",col_names = FALSE)
WTA <- Read10X("0329WTAplus")
tag.smp <- read_tsv("tag0-4.txt", col_names = T)

projectname <- c("dLN")
biotin <- read_tsv("biotintag_umi_tag.tsv")


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
VlnPlot(haha,features = "biotin",pt.size=0)+geom_boxplot(width=.2, col='black', fill="white")+geom_hline(aes(yintercept=2.0),colour="red",linetype="dashed",size=1)
ggsave('Biotin_sample.jpg', width = 15, height = 10)
ggsave('Biotin_sample.eps', width = 15, height = 10)

#Save RDS
DefaultAssay(haha) <- "RNA"
saveRDS(object = haha,file = "biotin-sample-dLN.rds")

rm(list=ls())
gc()


#------------------Clustering--------------#

library(Seurat)
library(tidyverse)
library(ggplot2)
library(corrplot)
library(dplyr)
library(grid)
library(cowplot)
PROX1 <- readRDS("biotin-sample-dLN.rds")

PROX1[["percent.mt"]] <- PercentageFeatureSet(PROX1, pattern = "^mt-") 

VlnPlot(PROX1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.5,ncol = 3)
ggsave(filename = "QC.png",width = 8, height = 5 ) 


PROXA2 <- NormalizeData(PROX1, normalization.method = "LogNormalize")
all.genes <- rownames(PROXA2)
PROXA2 <- ScaleData(PROXA2, features = all.genes)


PROXA2 <- FindVariableFeatures(PROXA2, selection.method = "vst", nfeatures = 2000)
PROXA2 <- RunPCA(PROXA2, features = VariableFeatures(object = PROXA2), verbose = FALSE)

ElbowPlot(PROXA2)
PROXA2 <- FindNeighbors(PROXA2, dims = 1:20)
PROXA3 <- FindClusters(PROXA2, resolution = 1.2)

#Run tsNE
set.seed(1234)
PROXA3 <- RunTSNE(PROXA3, dims = 1:20)
DimPlot(PROXA3, reduction = "tsne", label = T )
ggsave(filename = "all_tsne.png",width = 8, height = 5 )

#Run umap
PROXA3 <- RunUMAP(PROXA3, dims = 1:20)
DimPlot(PROXA3, reduction = 'umap', label = T)
ggsave(filename = "all_umap.png", width = 8, height = 5)

PROXA3.markers <- FindAllMarkers(PROXA3, 
                                 only.pos = TRUE,
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)
top10 <- PROXA3.markers %>% group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "top10 gene.csv")

#Annotation
celltype <- c("B cells", "B cells","B cells","T cells","T cells","T cells","B cells","T cells","T cells","T cells","T cells","B cells","Proliferating cells","MP","Plasma cells","B cells","B cells","NK")
names(celltype) <- levels(PROXA3)
PROXA3 <- RenameIdents(PROXA3, celltype)
PROXA3[["celltype"]]<- Idents(object = PROXA3)

DimPlot(PROXA3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "rename umap.png")
DimPlot(PROXA3, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "sample") + NoLegend()
saveRDS(PROXA3, 'rename umap all1.rds')




# Remove ICB and label2 (label1 and label2,ICB1 and ICB2 are biological repeats)

Iso_label1_ICB1_dln <- subset(PROXA3, sample=="isotype"|sample=="label1"|sample=="ICB1")
DimPlot(Iso_label1_ICB1_dln, reduction = "umap", pt.size = 0.5,label = T)+NoLegend()
ggsave("Iso_label1_ICB1 umap.eps")
ggsave("Iso_label1_ICB1 umap.png")

VlnPlot(Iso_label1_ICB1_dln, features = "biotin",pt.size=0, group.by = "celltype",split.by = "sample")

#Marker genes for cell clusters

feature <- c( "Cd79a","Cd79b","Ms4a1","Bank1","Mki67","Ccnb2","Pclaf","Top2a","Fcer1g","Fscn1",
              "Ifitm3","Lyz2", "Cd68","Trac","Trbc1","Cd3d","Cd3e", "Cd3g","Klra4","Klrd1","Klre1",
              "Ncr1",   "Nkg7",   "Iglv1",  "Xbp1",   "Jchain", "Iglc1","Ighg1")
DotPlot(Iso_label1_ICB1_dln, features = feature)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

ggsave("Marker genes for cell clusters.eps")
ggsave("Marker genes for cell clusters.png")

saveRDS(Iso_label1_ICB1_dln,"Iso_label1_ICB1_dln.rds")


#Violin plot showing the biotin labeling signal of B cell and T cell clusters

B.T.Iso_label1_ICB1 <- subset(Iso_label1_ICB1_dln, celltype=="B cells"|celltype=="T cells")
VlnPlot(B.T.Iso_label1_ICB1, features = "biotin",pt.size=0, group.by = "celltype",split.by = "sample")





