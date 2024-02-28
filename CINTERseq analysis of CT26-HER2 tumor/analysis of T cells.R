

#-----------Analysis of T cells-----------# 
library(Seurat)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)

PROXA3 <- readRDS("rename umap all.rds")

#T cell subclusters annotation

Tcell <- subset(PROXA3,idents = 'T cells')
set.seed(1234)
Tcell <- FindNeighbors(Tcell, reduction = "pca", dims = 1:15) 
Tcell <- FindClusters(Tcell,resolution = 0.6, algorithm = 1) 
Tcell <- RunTSNE(object=Tcell,dims.use=1:15, do.fast=TRUE,check_duplicates = FALSE)
Tcell <- RunUMAP(Tcell, reduction = "pca", dims = 1:15)
DimPlot(Tcell,label = TRUE)  

T_cluster <- c("Naive/Memory CD4","Cytotoxic effector CD8","Cytotoxic effector CD8","Activation/dysfunction CD8","Treg","γδT","Effector memory CD4","Activation/dysfunction CD8")
names(T_cluster) <- levels(Tcell)
Tcell <- RenameIdents(Tcell, T_cluster)
Tcell[["T_cluster"]]<- Idents(object = Tcell)

DimPlot(Tcell, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(filename = "T cell clusters/Tcell rename.png",width =12, height = 8)
ggsave(filename = "T cell clusters/Tcell rename.eps")

saveRDS(Tcell, 'T cell clusters/T cell rename_total.rds')

label1_T <- subset(Tcell,subset=sample=="label1")

saveRDS(label1_T, 'T cell clusters/label1_T.rds')

#Marker genes expressed for T cell clusters

Tcell.marker <- c("Mki67","Ifng","Pdcd1","Ctla4","Foxp3","Tcf7","Il7r","Tnfrsf9","Gzma","Gzmb","Nkg7","Lag3","Trdc","Ccl5")   

DotPlot(label1_T, features = Tcell.marker) 
ggsave(filename = "T_cluster functional markers.eps")
ggsave(filename = "T_cluster functional markers.png")


# Top genes expressed in T cell clusters
PROXA3.markers <- FindAllMarkers(label1_T, 
                                 only.pos = TRUE,
                                 min.pct = 0.25, 
                                 logfc.threshold = 0.25)
write.csv(PROXA3.markers, "T_cluster_markers.csv")

top10 <- PROXA3.markers %>% group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC)

library(viridis)
DoHeatmap(label1_T, features = top10$gene)+scale_fill_viridis(discrete=FALSE)

ggsave(filename = "T_cluster_top10 markers.png")
ggsave(filename = "T_cluster_top10 markers.eps")


#-----Cellular landscapes over Biotag signal

library(RColorBrewer)
display.brewer.all()

#biotin low:L3, biotin medium:L2, biotin high:L1
srt_corr <- label1_T
srt_corr@meta.data$level.biotin <- "NA"                     
srt_corr@meta.data$level.biotin[1:233] <-  "L3"                    
srt_corr@meta.data$level.biotin[234:466] <-  "L2"
srt_corr@meta.data$level.biotin[467:701] <-  "L1"

VlnPlot(srt_corr,features = "biotin",pt.size=0,group.by = "level.biotin")
ggsave(filename = "Tcell biolevel.png")
ggsave(filename = "Tcell biolevel.eps")

Idents(srt_corr) <- 'T_cluster'
table(srt_corr@meta.data[["level.biotin"]])
yy <- table(Idents(srt_corr),srt_corr@meta.data[["level.biotin"]])

col <- brewer.pal(9, "Set3")
Cellratio <- prop.table(yy, margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1)) 
ggplot(Cellratio) +
  geom_bar(aes(x = Freq, y = Var2, fill =Var1),stat= 'identity',width=0.7,size=0.5,colour = '#222222')+
  theme_classic() +
  labs(x= 'Ratio',y= 'Sample')+
  coord_flip()+
  scale_fill_manual(values = col)+
  theme(panel.border = element_rect(fill=NA,color='black',size=0.5,linetype='solid'))
ggsave(filename = "Tcell.table.png")
ggsave(filename = "Tcell.table.eps")


#UMAP plot showing T cell landscapes over Biotag signal
DimPlot(srt_corr, reduction = "umap", pt.size =2, split.by = "level.biotin",cols=col)
ggsave(filename = "Tcell_split.png")
ggsave(filename = "Tcell_split.eps")


#----Correlation analysis of total T cells

srt_corr <- label1_T 

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

write.csv(corr, "total T_biotin correlation.csv")


t.corr <- corr
t.corr <- t.corr[order(t.corr$rho),]


corr.point <-  function(corr, label.gene, title="T cell biotin corr"){
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
genelist1 <- c("Lag3","Tnfrsf9","Havcr2",
               "Tcf7","Malat1","Txnip","Klf2")

corr.point(corr = t.corr, genelist1)

ggsave("T cell biotin corr.eps")



#------Correlation analysis of γδT

T_γ <- subset(label1_T, T_cluster=="γδT")


srt_corr <- T_γ 

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

write.csv(corr, "γδT_biotin correlation.csv")


γδT.corr <- corr

γδT.corr <- γδT.corr[order(γδT.corr$rho),]


corr.point <-  function(corr, label.gene, title="γδT.corr biotin corr"){
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
Pos10 <- filter(γδT.corr, γδT.corr$rho >= 0.3412)
Pos10$gene2
Neg10 <- filter(γδT.corr, γδT.corr$rho <= -0.2850)
Neg10$gene2


genelist1 <- c("Stat3","Chmp1a","Crnkl1","Elmo1","Pik3ap1", "H2az2","S100a9","Tnfrsf9", "Edf1", "Cpa3","Psma5",
               "Whamm","Grk6","Hnrnpa3","Fbxo25", "Ubn1", "Rbfa", "Pigx","Naa16","Rnf115","Klf4" )

corr.point(corr = γδT.corr  , genelist1)

ggsave("γδT.corr_biotin corr.eps")



#-------Correlation analysis of Treg

Treg <- subset(label1_T, T_cluster=="Treg")

srt_corr <- Treg 

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

write.csv(corr, "Treg_biotin correlation.csv")

Treg.corr <- corr

Treg.corr <- Treg.corr[order(Treg.corr$rho),]


corr.point <-  function(corr, label.gene, title="Treg biotin corr"){
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
Pos10 <- filter(Treg.corr, Treg.corr$rho >= 0.3365)
Pos10$gene2
Neg10 <- filter(Treg.corr, Treg.corr$rho <= -0.311)
Neg10$gene2


genelist1 <- c("Tnfrsf4","Lbh" ,"Ccdc88c", "Lag3","Sell","Pnpla8","Gmfg","Ly6a","Satb1","Glg1",
               "Lta","Tmsb4x","Lrch1","Ms4a4b", "Emp3", "Gm12216","Ucp2", "Ybx1")

corr.point(corr = Treg.corr, genelist1)

ggsave("Treg_biotin corr.eps", width = 6, height = 5.5)



#------Correlation analysis of total CD4

CD4 <- subset(label1_T, T_cluster=="Naive/Memory CD4"|T_cluster=="Effector memory CD4")
ncol(CD4)

srt_corr <- CD4 

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

write.csv(corr, "CD4_biotin correlation.csv")

CD4.corr <- corr
CD4.corr <- CD4.corr[order(CD4.corr$rho),]


corr.point <-  function(corr, label.gene, title="CD4biotin corr"){
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
genelist1 <- c("Lag3","Lgals1","Nrgn","Lgals7","Cd200","Nkg7","Runx2","Cd3g",
               "Klf2","Malat1","Tcf7","Txnip","Itga4","Klf3")

corr.point(corr = CD4.corr, genelist1)

ggsave("CD4_biotin corr.eps")




#--------Correlation analysis of total CD8

CD8 <- subset(label1_T, T_cluster=="Cytotoxic effector CD8"|T_cluster=="Activation/dysfunction CD8")
ncol(CD8)

srt_corr <- CD8

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

write.csv(corr, "CD8_biotin correlation.csv")

CD8.corr <- corr
CD8.corr <- CD8.corr[order(CD8.corr$rho),]


corr.point <-  function(corr, label.gene, title="CD8 biotin corr"){
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
genelist1 <- c("Tnfrsf9","Eif3l","Il2ra","Zfp36l1","Havcr2","Lag3","Stat5b","Runx2",
               "Ctla4","Stat3","Tcf7","Ctsw")

corr.point(corr = CD8.corr, genelist1)

ggsave("CD8_biotin corr.eps")





#--------Correlation analysis of Naive/Memory CD4

CD4_naive.memory <- subset(label1_T, T_cluster=="Naive/Memory CD4")
ncol(CD4_naive.memory)

srt_corr <- CD4_naive.memory 

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

write.csv(corr, "CD4_naive.memory_biotin correlation.csv")


CD4_nai.mem.corr <- corr

CD4_nai.mem.corr <- CD4_nai.mem.corr[order(CD4_nai.mem.corr$rho),]


corr.point <-  function(corr, label.gene, title="CD4_nai.mem biotin corr"){
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
Pos10 <- filter(CD4_nai.mem.corr, CD4_nai.mem.corr$rho >= 0.2121)
Pos10$gene2
Neg10 <- filter(CD4_nai.mem.corr, CD4_nai.mem.corr$rho <= -0.1889)
Neg10$gene2

#highlighted genes
genelist1 <- c("Lgals1","Ifitm1","Ifitm3","Bhlhe40","Rabgap1l","Cd52","Jdp2","Lag3", "Ahcyl1","Runx2","Tnfrsf19",
               "Bzw2","Tex30","Zrsr2","Rab3gap1", "Ctdsp2","Smim10l1", "Ubc","Itm2b", "Crbn", "H3f3b","Itga4","Txnip","Klf3","Tcf7")

corr.point(corr = CD4_nai.mem.corr, genelist1)

ggsave("CD4_nai.mem_biotin corr.eps")








#--------Correlation analysis of effector memory like CD4

CD4_effector.mem <- subset(label1_T, T_cluster=="Effector memory CD4")
ncol(CD4_effector.mem)

srt_corr <- CD4_effector.mem

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

write.csv(corr, "CD4_effector.mem_biotin correlation.csv")


CD4_effector.mem.corr <- corr


CD4_effector.mem.corr <- CD4_effector.mem.corr[order(CD4_effector.mem.corr$rho),]


corr.point <-  function(corr, label.gene, title="CD4_effector.mem biotin corr"){
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
Pos10 <- filter(CD4_Exh.corr, CD4_Exh.corr$rho >= 0.39475)
Pos10$gene2
Neg10 <- filter(CD4_Exh.corr, CD4_Exh.corr$rho <= -0.3074)
Neg10$gene2


genelist1 <- c("Skp1a","Iscu","Mark2","Cox5b","Smim19", "Nfatc1", "Sec23b","Dctn3","Yaf2","Zyx" ,
               "Pitpna","Dpy19l1","Tbc1d1","Znrf2","Xpr1","Zeb1","Ring1","Med13l","Nek9","Fuca2")

corr.point(corr = CD4_effector.mem.corr, genelist1)

ggsave("CD4_effector.mem.corr.eps")


#-----------------------TCR clone analysis (CDR3 a chaim)

library(scRepertoire)
ori_TCR <- read_xlsx("a.chain.xlsx")
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

label1_T_TCR <- combineExpression(com_TCR, 
                                  label1_T, 
                                  cloneCall="aa",proportion = F,
                                  cloneTypes = c(None=0,Single = 1, 
                                                 Small = 3, Medium = 5, 
                                                 expanded = 10))
slot(label1_T_TCR, "meta.data")$cloneType <- factor(slot(label1_T_TCR, "meta.data")$cloneType, 
                                                    levels = c("expanded (5 < X <= 10)",  
                                                               "Medium (3 < X <= 5)", 
                                                               "Small (1 < X <= 3)", 
                                                               "Single (0 < X <= 1)",
                                                               NA))
DimPlot(label1_T_TCR, group.by = "cloneType",pt.size = 2) +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")
ggsave("TCR_achain.clonetype.png")
ggsave("TCR_achain.clonetype.eps")

ncol(label1_T_TCR)
#Median value of biotin in total T cells
label1_T_TCR@assays$biotin@counts@x %>% median()


#Non-TCR cells removed
label1_T_TCR1 <- subset(label1_T_TCR,subset=cloneType != "NA")
DimPlot(label1_T_TCR1, group.by = "cloneType",pt.size = 1) +
  scale_color_manual(values = colorblind_vector(5), na.value="grey")


# TCR clone distribution over Biotag signal
haha <- label1_T_TCR1
haha@meta.data$TCR_AA <- haha@meta.data$CTaa
haha@meta.data[which(haha@meta.data$Frequency == '2'),"TCR_AA"] <- 'double clone' 
haha@meta.data[which(haha@meta.data$Frequency == '1'),"TCR_AA"] <- 'single clone' 
VlnPlot(haha, "biotin", group.by = "TCR_AA")
Idents(haha) <- haha@meta.data$TCR_AA
shunxu <- c("CAQSRGGRALIF_NA", "CALSDLEGNNRIFF_NA","CAASDSNNRIFF_NA","CAMRGYTGANTGKLTF_NA",
            "CALGHQGGRALIF_NA","CAARVWKTASLGKLQF_NA","CALGGYGGSGNKLIF_NA",
            "CAMRNRGSALGRLHF_NA","CARGYGGSGNKLIF_NA","CAIGGDSNYQLIW_NA",
            "CALGRTGNTRKLIF_NA", "CALSERGSALGRLHF_NA", "CILRDGNEKITF_NA","CSASSGGNYKYVF_NA",
            "double clone", "single clone")
shunxu <- rev(shunxu)
Idents(haha) <- factor(Idents(haha), levels = shunxu)
VlnPlot(haha, "biotin", pt.size = 0.5)+
  coord_flip()+
  geom_hline(aes(yintercept=2.045323),
             colour="red",linetype="dashed",size=1)+
  NoLegend()

ggsave("TCR-achain distribution.png")
ggsave("TCR-achain distribution.eps")








