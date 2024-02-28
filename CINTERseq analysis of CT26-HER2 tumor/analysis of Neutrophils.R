#-----------Analysis of Neutrophils----------#
library(Seurat)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)

label1 <- readRDS("label1.rds")

#clustering
library(RColorBrewer)
col <- brewer.pal(6, "Accent")


label1_Neutrophils <- subset(label1,idents = 'Neutrophils')
set.seed(1234)
label1_Neutrophils <- FindNeighbors(label1_Neutrophils, reduction = "pca", dims = 1:15) 
label1_Neutrophils <- FindClusters(label1_Neutrophils,resolution = 0.6, algorithm = 1) 
label1_Neutrophils <- RunTSNE(object=label1_Neutrophils,dims.use=1:15, do.fast=TRUE,check_duplicates = FALSE)
label1_Neutrophils <- RunUMAP(label1_Neutrophils, reduction = "pca", dims = 1:15)
DimPlot(label1_Neutrophils,label = TRUE, cols = col)  
ggsave(filename = "Neutrophils clusters/label1_Neutrophils clusters.png")
ggsave(filename = "Neutrophils clusters/label1_Neutrophils clusters.eps")

saveRDS(label1_Neutrophils,"label1_Neutrophils.rds")

# Top genes expressed in neutrophils clusters
label1_Neutrophils.markers <- FindAllMarkers(label1_Neutrophils, 
                                             only.pos = TRUE,
                                             min.pct = 0.25, 
                                             logfc.threshold = 0.25)
top10 <- label1_Neutrophils.markers %>% group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "Neutrophils clusters/label1_Neutrophils top10markers.csv")

label1_Neutrophils_top10markers <- read_csv("Neutrophils clusters/label1_Neutrophils top10markers.csv")

marker <- unique(label1_Neutrophils_top10markers$gene)

marker1 <- marker %>% grep(pattern="^Rp", 
                           value = T, 
                           invert = T) %>% grep(pattern="mt", 
                                                value = T,
                                                invert = T)

DotPlot(label1_Neutrophils, features = marker1)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

ggsave("Neutrophils clusters/dotplot Neutrophils marker.png")
ggsave("Neutrophils clusters/dotplot Neutrophils marker.eps")

#--------Cellular landscapes over Biotag signal

#biotin low:L3, biotin medium:L2, biotin high:L1
srt_corr <- label1_Neutrophils
srt_corr@meta.data$level.biotin <- "NA"                     
srt_corr@meta.data$level.biotin[1:183] <-  "L3"                    
srt_corr@meta.data$level.biotin[184:367] <-  "L2"
srt_corr@meta.data$level.biotin[368:551] <-  "L1"


VlnPlot(srt_corr,features = "biotin",pt.size=0,group.by = "level.biotin")
ggsave(filename = "Neutropils biolevel.png")
ggsave(filename = "Neutropils biolevel.eps")

Idents(srt_corr) <- 'seurat_clusters'
table(srt_corr@meta.data[["level.biotin"]])
yy <- table(Idents(srt_corr),srt_corr@meta.data[["level.biotin"]])

col <- brewer.pal(6, "Accent")
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
ggsave(filename = "Neutropils.table.png")
ggsave(filename = "Neutropils.table.eps")


#------------Correlation analysis 

srt_corr <- label1_Neutrophils 
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

write.csv(corr, "Neutrophils_biotin correlation.csv")


neu.corr <- corr
neu.corr <- neu.corr[order(neu.corr$rho),]


corr.point <-  function(corr, label.gene, title="Neutrophils biotin corr"){
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
Pos10 <- filter(neu.corr, neu.corr$rho >= 0.2409)
Pos10$gene2
Neg10 <- filter(neu.corr, neu.corr$rho <= -0.20267)
Neg10$gene2

genelist1 <- c("Fth1","Ccl4","Ccl3","Nfkbia","Il1rn","Cd274","Tnfrsf1b","Cxcl2","Clec4e","Hilpda",
               "Actb","Mmp8","Retnlg","Pglyrp1","Mtus1","Msrb1","Tmsb4x","Mmp9","Taldo1","Lsp1")



corr.point(corr = neu.corr, genelist1)
ggsave("Neutrophils biotin corr.eps")


#GO enrichment analysis of the top 100 genes positively correlated with biotin signal
library(org.Mm.eg.db)

pos.gene <- filter(neu.corr, neu.corr$rho > 0)

pos.gene.top100 <- pos.gene %>% top_n(n = 100, wt = rho)

ego_BP <- enrichGO(gene          = pos.gene.top100$gene2,
                   #universe     = row.names(dge.cell.types),
                   OrgDb         = 'org.Mm.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)

dotplot(ego_BP ,showCategory=10)+ ggtitle("barplot for Biological process")


#Representative genes positively and negatively correlated with biotin signal
p<-subset(label1_Neutrophils, Ccl4>0) %>% FeatureScatter(feature1 = "biotin", 
                                                          feature2 = "Ccl4", 
                                                          span = 1, 
                                                          jitter = F,
                                                          pt.size = 0,
                                                          group.by = "orig.ident",
                                                          cols = "gray")
p+ labs(x = "Biotin")


p<-subset(label1_Neutrophils, Il1rn>0) %>% FeatureScatter(feature1 = "biotin", 
                                                         feature2 = "Il1rn", 
                                                         span = 1, 
                                                         jitter = F,
                                                         pt.size = 0,
                                                         group.by = "orig.ident",
                                                         cols = "gray")
p+ labs(x = "Biotin")


p<-subset(label1_Neutrophils, Mmp8>0) %>% FeatureScatter(feature1 = "biotin", 
                                                          feature2 = "Mmp8", 
                                                          span = 1, 
                                                          jitter = F,
                                                          pt.size = 0,
                                                          group.by = "orig.ident",
                                                          cols = "gray")
p+ labs(x = "Biotin")

p<-subset(label1_Neutrophils, Retnlg>0) %>% FeatureScatter(feature1 = "biotin", 
                                                         feature2 = "Retnlg", 
                                                         span = 1, 
                                                         jitter = F,
                                                         pt.size = 0,
                                                         group.by = "orig.ident",
                                                         cols = "gray")
p+ labs(x = "Biotin")




#Normalized expression of selected genes
FeaturePlot(label1_Neutrophils, features = c("Cstb","Hilpda","Gadd45b","Il1rn","Ccl4","Cd274"),ncol = 3)


