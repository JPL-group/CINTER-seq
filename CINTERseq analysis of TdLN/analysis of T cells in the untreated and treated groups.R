#------------------------------Analysis of T cells----------------------# 
library(Seurat)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)

PROXA3 <- readRDS('rename umap all1.rds')


#--------T cell subclusters clustering

Tcell <- subset(PROXA3,idents = 'T cells')

set.seed(2333)
Tcell <- NormalizeData(Tcell, normalization.method = "LogNormalize")
all.genes <- rownames(Tcell)
Tcell <- ScaleData(Tcell, features = all.genes)


Tcell <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 2000)
Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell), verbose = FALSE)

ElbowPlot(Tcell)
Tcell <- FindNeighbors(Tcell, dims = 1:15)
Tcell1 <- FindClusters(Tcell, resolution = 1.6)

#Run tsNE
Tcell1 <- RunTSNE(Tcell1, dims = 1:15)
DimPlot(Tcell1, reduction = "tsne", label = T )


#Run umap
Tcell1 <- RunUMAP(Tcell1, dims = 1:15)
DimPlot(Tcell1, reduction = 'umap', label = T)


#remove cluster 6, 9 and 12
FeaturePlot(Tcell1, features = c("Cd79a","Cd4","Cd8a","Cd3e"), label = T)
Tcell2 <- subset(Tcell1, idents = c('6','9','12'), invert=T)

#Re-clustering
haha <- Tcell2
haha= NormalizeData(haha, normalization.method = "LogNormalize")
all.genes = rownames(haha)
haha = ScaleData(haha, features = all.genes)

haha = FindVariableFeatures(haha, selection.method = "vst", nfeatures = 2000)
haha = RunPCA(haha, features = VariableFeatures(object = haha), verbose = FALSE)

ElbowPlot(haha)
haha = FindNeighbors(haha, dims = 1:15)
haha = FindClusters(haha, resolution = 1)

#tsNE降维
haha <- RunTSNE(haha, dims = 1:15)
DimPlot(haha, reduction = "tsne", label = T )
ncol(haha)

#umap降维
haha <- RunUMAP(haha, dims = 1:15)
DimPlot(haha, reduction = 'umap', label = T)

n=length(unique(haha@meta.data$seurat_clusters))
celltype=data.frame(ClusterID=0:(n-1),
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,5,13),2]='cluster 0'
celltype[celltype$ClusterID %in% c(1,3,10),2]='cluster 1'
celltype[celltype$ClusterID %in% c(2),2]='cluster 2'
celltype[celltype$ClusterID %in% c(6,9,12),2]='cluster 3'
celltype[celltype$ClusterID %in% c(7,8),2]='cluster 4'
celltype[celltype$ClusterID %in% c(4,11),2]='cluster 5'
haha@meta.data$cluster <- "NA"
for(i in 1:nrow(celltype)){
  haha@meta.data[which(haha@meta.data$seurat_clusters == celltype$ClusterID[i]),'cluster'] <- celltype$celltype[i]}
table(haha@meta.data$cluster)
DimPlot(haha, group.by = "cluster", label = T)

saveRDS(haha,"Tcell.reclustering.rds")


#umap plot of untreated group (label1) and aPD-1 group (ICB1) 
Tcell <- readRDS("Tcell.reclustering.rds")
Idents(Tcell) <- "cluster"
# Remove ICB and label2 (label1 and label2,ICB1 and ICB2 are biological repeats)
T_iso_lable1_ICB1 <- subset(Tcell,sample=="isotype"|sample=="label1"|sample=="ICB1")

label1_ICB1 <- subset(T_iso_lable1_ICB1,sample=="isotype",invert=T)
DimPlot(label1_ICB1, pt.size=0.5, label=T, split.by="sample")
DimPlot(label1_ICB1, pt.size=0.5, label=T)
#Biotag signal projected over the UMAP plot
FeaturePlot(label1_ICB1, "biotin", cols = c("#7fa4b8","#e5ecc5", "#c03932"),   max.cutoff =2.2, min.cutoff = 1)

#Top genes expressed in T cell clusters
PROXA3 <- T_iso_lable1_ICB1
Idents(PROXA3 ) <- "cluster"
PROXA3.markers <- FindAllMarkers(PROXA3, 
                                 only.pos = TRUE,
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)

top5 <- PROXA3.markers %>% group_by(cluster) %>%  top_n(n = 5, wt = avg_log2FC)

marker <- c("Tnfrsf4","Ctla4","Icos","Maf","S100a6","Cd8a","Klf2","Fos","Nkg7","Dusp2","Sltm","Rn7sk","Slfn1",  
             "Ly6a","Bst2","Zbp1","Ly6c2","Samd9l","Lef1","Cited2","Dusp10","Igfbp4","Nr4a1","Cd4","Ccl5",
            "Klrc1" ,"Gzmm","Klrk1","Il2rb")
DotPlot(PROXA3, features = marker)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

ggsave("dotplot T marker.png")
ggsave("dotplot T marker.eps")


#Proportion of cells
table(label1_ICB1@meta.data[["sample"]])
yy <- table(Idents(label1_ICB1),label1_ICB1@meta.data[["sample"]])

Cellratio <- prop.table(yy, margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1)) 
ggplot(Cellratio) +
  geom_bar(aes(x = Freq, y = Var2, fill =Var1),stat= 'identity',width=0.7,size=0.5,colour = '#222222')+
  theme_classic() +
  labs(x= 'Ratio',y= 'Sample')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color='black',size=0.5,linetype='solid'))
ggsave(filename = "T_prop-1.table.png")
ggsave(filename = "T_prop-1.table.eps")



#-------------------TCR clone analysis (CDR3 a chaim)

library(scRepertoire)
ori_TCR <- read_xlsx("achain.dLN.xlsx")
com_TCR <- combineTCR(ori_TCR,samples = "TCR", ID = "TCR",cells ="T-AB")
lengthContig(com_TCR, cloneCall="nt", chain = "both")
#ggsave("TCR/lengthContig.pdf", plot = lengthContig, width = 12, height = 10)
#ggsave("TCR/lengthContig.png", plot = lengthContig, width = 12, height = 10)
clonalHomeostasis(com_TCR, cloneCall = "aa",cloneTypes = c(rare = 0.001,  Medium = 0.01,expanded = 1))
#ggsave("TCR/clonalHomeostasis.pdf", plot = umap_TBNK, width = 12, height = 10)
#ggsave("TCR/clonalHomeostasis.png", plot = umap_TBNK, width = 12, height = 10)
clonalProportion(com_TCR, cloneCall = "aa") 

com_TCR$TCR_TCR$barcode<- str_remove_all(string = com_TCR$TCR_TCR$barcode, pattern = "TCR_TCR_")


colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

T_iso_lable1_ICB1 <- readRDS("T_iso_lable1_ICB1.rds")

##dLN.label1 (untreated) TCR clone

dLN.label1 <- subset(T_iso_lable1_ICB1, sample=="label1")

table(dLN.label1@meta.data$sample)


dLN.label1_TCR <- combineExpression(com_TCR, 
                                    dLN.label1, 
                                    cloneCall="aa",proportion = F,
                                    cloneTypes = c(None=0,Single = 1, 
                                                   Small = 3, Medium = 5, 
                                                   expanded = 30))
slot(dLN.label1_TCR, "meta.data")$cloneType <- factor(slot(dLN.label1_TCR, "meta.data")$cloneType, 
                                                      levels = c("expanded (5 < X <= 30)",  
                                                                 "Medium (3 < X <= 5)", 
                                                                 "Small (1 < X <= 3)", 
                                                                 "Single (0 < X <= 1)",
                                                                 NA))

colorblind_vector <- colorRampPalette(c( "#7AC5FF", "#0348A6"))
DimPlot(dLN.label1_TCR, group.by = "cloneType",pt.size = 2, split.by = "sample") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")

ggsave("dLN.label1_TCR_achain.clonetype.png")
ggsave("dLN.label1_TCR_achain.clonetype.eps")

##dLN.ICB1 (aPD-1) TCR clone
dLN.ICB1 <- subset(T_iso_lable1_ICB1, sample=="ICB1")

table(dLN.ICB1@meta.data$sample)


dLN.ICB1_TCR <- combineExpression(com_TCR, 
                                  dLN.ICB1, 
                                  cloneCall="aa",proportion = F,
                                  cloneTypes = c(None=0,Single = 1, 
                                                 Small = 3, Medium = 5, 
                                                 expanded = 30))
slot(dLN.ICB1_TCR, "meta.data")$cloneType <- factor(slot(dLN.ICB1_TCR, "meta.data")$cloneType, 
                                                    levels = c("expanded (5 < X <= 30)",  
                                                               "Medium (3 < X <= 5)", 
                                                               "Small (1 < X <= 3)", 
                                                               "Single (0 < X <= 1)",
                                                               NA))

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
DimPlot(dLN.ICB1_TCR, group.by = "cloneType",pt.size = 2, split.by = "sample") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")

ggsave("dLN.ICB1_TCR_achain.clonetype.png")
ggsave("dLN.ICB1_TCR_achain.clonetype.eps")


saveRDS(dLN.ICB1_TCR,"dLN.ICB1_TCR.achain.rds")

dLN.ICB1_TCR <- readRDS("dLN.ICB1_TCR.achain.rds")


#median value of biotin signal
dLN.ICB1_TCR@assays$biotin@counts@x %>% median()


#Non-TCR cells removed
dLN.ICB1_TCR.have <- subset(dLN.ICB1_TCR,subset=cloneType != "NA")
DimPlot(dLN.ICB1_TCR.have, group.by = "cloneType",pt.size = 1) +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")

ggsave("dLN.ICB1_TCR_have.achain.clonetype.png")
ggsave("dLN.ICB1_TCR_have.achain.clonetype.eps")

saveRDS(dLN.ICB1_TCR.have,"dLN.ICB1_TCR.have.achain.rds")

dLN.ICB1_TCR.have <- readRDS("dLN.ICB1_TCR.have.achain.rds")


# TCR clone distribution over Biotag signal
haha <- dLN.ICB1_TCR.have
haha@meta.data$TCR_AA <- haha@meta.data$CTaa
haha@meta.data[which(haha@meta.data$Frequency == '2'),"TCR_AA"] <- 'double clone' 
haha@meta.data[which(haha@meta.data$Frequency == '1'),"TCR_AA"] <- 'single clone' 
VlnPlot(haha, "biotin", group.by = "TCR_AA")
table(haha@meta.data$TCR_AA)

Idents(haha) <- haha@meta.data$TCR_AA
shunxu <- c("CAADSNYQLIW_NA", "CAASASSGSWQLIF_NA","CAASANSGTYQRF_NA","CAMREGRTGNYKYVF_NA",
            "CVVGDRGSALGRLHF_NA", "double clone", "single clone")
shunxu <- rev(shunxu)
Idents(haha) <- factor(Idents(haha), levels = shunxu)
VlnPlot(haha, "biotin", pt.size = 0.5)+
  coord_flip()+
  geom_hline(aes(yintercept=1.60206),
             colour="red",linetype="dashed",size=1)+
  NoLegend()

ggsave("dLN.TCR-achain distribution.png")
ggsave("dLN.TCR-achain distribution.eps")







#-----------------Correlation analysis of ICB.CD4(aPD-1)


ICB1 <- subset(T_iso_lable1_ICB1, sample=="ICB1")
label1 <- subset(T_iso_lable1_ICB1, sample=="label1")


ICB1.CD4 <- subset(ICB1, cluster=="cluster 1"|cluster=="cluster 2"|cluster=="cluster 4")


srt_corr <- ICB1.CD4

srt_corr@assays$RNA@scale.data <- rbind(as.matrix(srt_corr@assays$biotin@data),
                                        as.matrix(srt_corr@assays$RNA@scale.data))

srt_corr.scale <- srt_corr@assays$RNA@scale.data
View(srt_corr.scale) 
srt_corr_corr <- correlatePairs(srt_corr.scale, pairings=list(1,2:55472))

corr <- srt_corr_corr %>%as.data.frame()
genelist <- corr$gene2 %>% grep(pattern="^Rp", 
                                value = T, 
                                invert = T) %>% grep(pattern="mt", 
                                                     value = T,
                                                     invert = T)
corr <- subset(corr, corr$gene2 %in% genelist)
corr <- subset(corr, corr$rho != "NaN")

write.csv(corr, "ICB1.CD4_biotin.correlation.csv")

ICB1.CD4.corr <- corr
ICB1.CD4.corr <- ICB1.CD4.corr[order(ICB1.CD4.corr$rho),]


corr.point <-  function(corr, label.gene, title="ICB1.CD4 biotin corr"){
  library(ggplot2)
  library(ggrepel)
  ggplot(corr, aes(x = 1:nrow(corr), y= rho))+
    geom_point()+
    geom_text_repel(aes(label= ifelse(gene2 %in% label.gene, gene2, "")),
                    color="red", min.segment.length = 0.1,
                    nudge_x = 100, nudge_y = 0, max.overlaps = 3000
    )+
    labs(x="genes", y = "Correlation value", title = title)+
    geom_point(aes(x = ifelse(gene2 %in% label.gene, 1:nrow(corr), NA), 
                   y = ifelse(gene2 %in% label.gene, rho, NA)),
               color = "red", size = 2)}


#highlighted genes
Pos10 <- filter(ICB1.CD4.corr, ICB1.CD4.corr$rho >= 0.1789)
Pos10$gene2
Neg10 <- filter(ICB1.CD4.corr, ICB1.CD4.corr$rho <= -0.0977)
Neg10$gene2

genelist1 <- c("Il2rb","Lgals1","S100a6","Snx2","Maf","Itgb1","Sri","Mapre2","Nkg7","Rhd",
               "Vps37b","Gimap6","Txnip","Klf2","Urb1","Lef1","Limk2","Pnrc1","Eef1a1","Tpt1")

corr.point(corr = ICB1.CD4.corr, genelist1)

ggsave("ICB1.CD4.corr.eps")
ggsave("ICB1.CD4.corr.png")



#-----------------Correlation analysis of ICB.CD8(aPD-1)

ICB1.CD8 <- subset(ICB1, cluster=="cluster 0"|cluster=="cluster 3"|cluster=="cluster 5")


srt_corr <- ICB1.CD8

srt_corr@assays$RNA@scale.data <- rbind(as.matrix(srt_corr@assays$biotin@data),
                                        as.matrix(srt_corr@assays$RNA@scale.data))

srt_corr.scale <- srt_corr@assays$RNA@scale.data
View(srt_corr.scale) 
srt_corr_corr <- correlatePairs(srt_corr.scale, pairings=list(1,2:55472))

corr <- srt_corr_corr %>%as.data.frame()
genelist <- corr$gene2 %>% grep(pattern="^Rp", 
                                value = T, 
                                invert = T) %>% grep(pattern="mt", 
                                                     value = T,
                                                     invert = T)
corr <- subset(corr, corr$gene2 %in% genelist)
corr <- subset(corr, corr$rho != "NaN")

write.csv(corr, "ICB1.CD8_biotin.correlation.csv")

ICB1.CD8.corr <- corr
ICB1.CD8.corr <- ICB1.CD8.corr[order(ICB1.CD8.corr$rho),]


corr.point <-  function(corr, label.gene, title="ICB1.CD8 biotin corr"){
  library(ggplot2)
  library(ggrepel)
  ggplot(corr, aes(x = 1:nrow(corr), y= rho))+
    geom_point()+
    geom_text_repel(aes(label= ifelse(gene2 %in% label.gene, gene2, "")),
                    color="red", min.segment.length = 0.1,
                    nudge_x = 100, nudge_y = 0, max.overlaps = 3000
    )+
    labs(x="genes", y = "Correlation value", title = title)+
    geom_point(aes(x = ifelse(gene2 %in% label.gene, 1:nrow(corr), NA), 
                   y = ifelse(gene2 %in% label.gene, rho, NA)),
               color = "red", size = 2)}


#highlighted genes
Pos10 <- filter(ICB1.cd8, ICB1.cd8$rho >= 0.15476)
Pos10$gene2
Neg10 <- filter(ICB1.cd8, ICB1.cd8$rho <= -0.0806)
Neg10$gene2

genelist1 <- c("Cd52","Rhd","Gzmm","Nudt16l1", "Samd3","Prr13","Tmsb4x","Il2rb","Ly6c2","Ccl5",
               "Klf2","Gm42418","Gm48432","Malat1","Calr","Tjap1","Jund","Armcx3","Rgcc","Igf1r")

corr.point(corr = ICB1.CD8.corr, genelist1)

ggsave("ICB1.CD8.corr.eps")
ggsave("ICB1.CD8.corr.png")




#-----------------Correlation analysis of label1.CD4(untreated)

label1.CD4 <- subset(label1, cluster=="cluster 1"|cluster=="cluster 2"|cluster=="cluster 4")


srt_corr <- label1.CD4

srt_corr@assays$RNA@scale.data <- rbind(as.matrix(srt_corr@assays$biotin@data),
                                        as.matrix(srt_corr@assays$RNA@scale.data))

srt_corr.scale <- srt_corr@assays$RNA@scale.data
View(srt_corr.scale) 
srt_corr_corr <- correlatePairs(srt_corr.scale, pairings=list(1,2:55472))

corr <- srt_corr_corr %>%as.data.frame()
genelist <- corr$gene2 %>% grep(pattern="^Rp", 
                                value = T, 
                                invert = T) %>% grep(pattern="mt", 
                                                     value = T,
                                                     invert = T)
corr <- subset(corr, corr$gene2 %in% genelist)
corr <- subset(corr, corr$rho != "NaN")

write.csv(corr, "label1.CD4_biotin.correlation.csv")

label1.CD4.corr <- corr
label1.CD4.corr <- label1.CD4.corr[order(label1.CD4.corr$rho),]


corr.point <-  function(corr, label.gene, title="label1.CD4 biotin corr"){
  library(ggplot2)
  library(ggrepel)
  ggplot(corr, aes(x = 1:nrow(corr), y= rho))+
    geom_point()+
    geom_text_repel(aes(label= ifelse(gene2 %in% label.gene, gene2, "")),
                    color="red", min.segment.length = 0.1,
                    nudge_x = 100, nudge_y = 0, max.overlaps = 3000
    )+
    labs(x="genes", y = "Correlation value", title = title)+
    geom_point(aes(x = ifelse(gene2 %in% label.gene, 1:nrow(corr), NA), 
                   y = ifelse(gene2 %in% label.gene, rho, NA)),
               color = "red", size = 2)}


#highlighted genes
Pos10 <- filter(label1.CD4.corr, label1.CD4.corr$rho >= 0.1856)
Pos10$gene2
Neg10 <- filter(label1.CD4.corr, label1.CD4.corr$rho <= -0.12765)
Neg10$gene2

genelist1 <- c("Surf1", "Maf","Itgae","Tnfrsf9","Ndfip1","S100a6","Tnfrsf18", "Srgn","Tnfrsf4","Icos",
               "Eef1a1","Cox7a2l", "Tpt1","Cxcr4", "Siah1a","Igfbp4",  "Txnip","Evl","Nacc2","Nde1")
corr.point(corr = label1.CD4.corr, genelist1)

ggsave("label1.CD4.corr.eps")
ggsave("label1.CD4.corr.png")


#-----------------Correlation analysis of label1.CD8(untreated)

label1.CD8 <- subset(label1, cluster=="cluster 0"|cluster=="cluster 3"|cluster=="cluster 5")


srt_corr <- label1.CD8

srt_corr@assays$RNA@scale.data <- rbind(as.matrix(srt_corr@assays$biotin@data),
                                        as.matrix(srt_corr@assays$RNA@scale.data))

srt_corr.scale <- srt_corr@assays$RNA@scale.data
View(srt_corr.scale) 
srt_corr_corr <- correlatePairs(srt_corr.scale, pairings=list(1,2:55472))

corr <- srt_corr_corr %>%as.data.frame()
genelist <- corr$gene2 %>% grep(pattern="^Rp", 
                                value = T, 
                                invert = T) %>% grep(pattern="mt", 
                                                     value = T,
                                                     invert = T)
corr <- subset(corr, corr$gene2 %in% genelist)
corr <- subset(corr, corr$rho != "NaN")

write.csv(corr, "label1.CD8_biotin.correlation.csv")

label1.CD8.corr <- corr
label1.CD8.corr <- label1.CD8.corr[order(label1.CD8.corr$rho),]


corr.point <-  function(corr, label.gene, title="label1.CD8 biotin corr"){
  library(ggplot2)
  library(ggrepel)
  ggplot(corr, aes(x = 1:nrow(corr), y= rho))+
    geom_point()+
    geom_text_repel(aes(label= ifelse(gene2 %in% label.gene, gene2, "")),
                    color="red", min.segment.length = 0.1,
                    nudge_x = 100, nudge_y = 0, max.overlaps = 3000
    )+
    labs(x="genes", y = "Correlation value", title = title)+
    geom_point(aes(x = ifelse(gene2 %in% label.gene, 1:nrow(corr), NA), 
                   y = ifelse(gene2 %in% label.gene, rho, NA)),
               color = "red", size = 2)}


#highlighted genes
Pos10 <- filter(label1.cd8, label1.cd8$rho >= 0.1714)
Pos10$gene2
Neg10 <- filter(label1.cd8, label1.cd8$rho <= -0.1168)
Neg10$gene2

genelist1 <- c("Grk2","Hopx","Dcaf13","Sh2d3c","Ccl5","Tmem160","Rap1a","Klra13-ps","Golim4","Ndufs3",
               "Eef1a1","Dapl1","Dtx1","Mta2","Tpt1","Tcf4","Cfap410","Zfp568", "Shkbp1","H2-T23" )

corr.point(corr = label1.CD8.corr, genelist1)
ggsave("label1.CD8.corr.eps")
ggsave("label1.CD8.corr.png")
