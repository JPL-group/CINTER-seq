#-----------Analysis of NK cells----------#

library(Seurat)
library(readr)
library(dplyr)
library(cowplot)
library(ggplot2)

label1 <- readRDS("label1.rds")

#clustering
library(RColorBrewer)
col <- brewer.pal(5, "Set2")
label1_NK <- subset(label1,idents = 'NK cells')
set.seed(1234)
label1_NK <- FindNeighbors(label1_NK, reduction = "pca", dims = 1:15) 
label1_NK <- FindClusters(label1_NK,resolution = 0.6, algorithm = 1) 
label1_NK <- RunTSNE(object=label1_NK,dims.use=1:15, do.fast=TRUE,check_duplicates = FALSE)
label1_NK <- RunUMAP(label1_NK, reduction = "pca", dims = 1:15)
DimPlot(label1_NK,label = TRUE, cols = col)  
ggsave(filename = "NK cell clusters/label1_NK clusters.png")
ggsave(filename = "NK cell clusters/label1_NK clusters.eps")

saveRDS(label1_NK,'NK cell clusters/label1_NK.rds')


# Top genes expressed in NK clusters
label1_NK.markers <- FindAllMarkers(label1_NK, 
                                    only.pos = TRUE,
                                    min.pct = 0.25, 
                                    logfc.threshold = 0.25)
top10 <- label1_NK.markers %>% group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "NK cell clusters/label1_NK top10markers.csv")

label1_NK_top10markers <- read_csv("NK cell clusters/label1_NK top10markers.csv")
marker <- unique(label1_NK_top10markers$gene)

marker1 <- marker %>% grep(pattern="^Rp", 
                           value = T, 
                           invert = T) %>% grep(pattern="mt", 
                                                value = T,
                                                invert = T)

DotPlot(label1_NK, features = marker1)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

ggsave("NK cell clusters/dotplot NK marker.png")
ggsave("NK cell clusters/dotplot NK marker.eps")


#--------Cellular landscapes over Biotag signal

#biotin low:L3, biotin medium:L2, biotin high:L1
srt_corr <- label1_NK
ncol(srt_corr)

srt_corr@meta.data$level.biotin <- "NA"                     
srt_corr@meta.data$level.biotin[1:286] <-  "L3"                    
srt_corr@meta.data$level.biotin[287:571] <-  "L2"
srt_corr@meta.data$level.biotin[572:857] <-  "L1"


VlnPlot(srt_corr,features = "biotin",pt.size=0,group.by = "level.biotin")
ggsave(filename = "NK cells biolevel.png")
ggsave(filename = "NK cells biolevel.eps")

Idents(srt_corr) <- 'seurat_clusters'
table(srt_corr@meta.data[["level.biotin"]])
yy <- table(Idents(srt_corr),srt_corr@meta.data[["level.biotin"]])

col <- brewer.pal(5, "Set2")

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
ggsave(filename = "NK cells.table.png")
ggsave(filename = "NK cells.table.eps")

#------------Correlation analysis 

srt_corr <- label1_NK 
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

write.csv(corr, "NK_biotin correlation.csv")


NK.corr <- corr
NK.corr <- NK.corr[order(NK.corr$rho),]


corr.point <-  function(corr, label.gene, title="NK cells biotin corr"){
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
Pos10 <- filter(corr, corr$rho >= 0.168)
Pos10$gene2
Neg10 <- filter(corr, corr$rho <= -0.0896)
Neg10$gene2

genelist1 <- c( "Tnfrsf9","Adam8","Spp1","Ndufb6","Icam1","Otud5","Serpina3g","Gzmc","Zbtb32","Ssr4",
                "Ccl5","Gm42418","Malat1","Pald1", "Ptprc","Tspyl4","Jak1")


corr.point(corr = NK.corr, genelist1)
ggsave("NK biotin corr.eps")




#-----------Analysis of Macrophages----------#
#clustering
col <- brewer.pal(12, "Set3")[7:12] 

label1_Macrophages <- subset(label1,idents = 'Macrophages')
set.seed(1234)
label1_Macrophages <- FindNeighbors(label1_Macrophages, reduction = "pca", dims = 1:15)
label1_Macrophages <- FindClusters(label1_Macrophages,resolution = 0.6, algorithm = 1) 
label1_Macrophages <- RunTSNE(object=label1_Macrophages,dims.use=1:15, do.fast=TRUE,check_duplicates = FALSE)
label1_Macrophages <- RunUMAP(label1_Macrophages, reduction = "pca", dims = 1:15)
DimPlot(label1_Macrophages,label = TRUE,cols = col) 

ggsave(filename = "Macrophages clusters/label1_Macrophages clusters.png")
ggsave(filename = "Macrophages clusters/label1_Macrophages clusters.eps")

saveRDS(label1_Macrophages,'Macrophages clusters/label1_Macrophages.rds')

# Top genes expressed in Macrophage clusters
label1_Macrophages.markers <- FindAllMarkers(label1_Macrophages, 
                                             only.pos = TRUE,
                                             min.pct = 0.25, 
                                             logfc.threshold = 0.25)
top10 <- label1_Macrophages.markers %>% group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "Macrophages clusters/label1_Macrophages top10markers.csv")


label1_Macrophages_top10markers <- read_csv("Macrophages clusters/label1_Macrophages top10markers.csv")


marker <- unique(label1_Macrophages_top10markers$gene)

marker1 <- marker %>% grep(pattern="^Rp", 
                           value = T, 
                           invert = T) %>% grep(pattern="mt", 
                                                value = T,
                                                invert = T)

DotPlot(label1_Macrophages, features = marker1)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

ggsave("Macrophages clusters/dotplot Macrophages marker.png")
ggsave("Macrophages clusters/dotplot Macrophages marker.eps")


#--------Cellular landscapes over Biotag signal
#biotin low:L3, biotin medium:L2, biotin high:L1
srt_corr <- label1_Macrophages
ncol(srt_corr)

srt_corr@meta.data$level.biotin <- "NA"                     
srt_corr@meta.data$level.biotin[1:134] <-  "L3"                    
srt_corr@meta.data$level.biotin[135:268] <-  "L2"
srt_corr@meta.data$level.biotin[269:402] <-  "L1"


VlnPlot(srt_corr,features = "biotin",pt.size=0,group.by = "level.biotin")
ggsave(filename = "Macrophages biolevel.png")
ggsave(filename = "Macrophages biolevel.eps")

Idents(srt_corr) <- 'seurat_clusters'
table(srt_corr@meta.data[["level.biotin"]])
yy <- table(Idents(srt_corr),srt_corr@meta.data[["level.biotin"]])

col <- brewer.pal(12, "Set3")[7:12]

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
ggsave(filename = "Macrophages.table.png")
ggsave(filename = "Macrophages.table.eps")



#------------Correlation analysis 

srt_corr <- label1_Macrophages 
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

write.csv(corr, "Macrophages_biotin correlation.csv")


Macro.corr <- corr
Macro.corr <- Macro.corr[order(Macro.corr$rho),]


corr.point <-  function(corr, label.gene, title="Macrophages biotin corr"){
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
Pos10 <- filter(corr, corr$rho >= 0.257)
Pos10$gene2
Neg10 <- filter(corr, corr$rho <= -0.2031)
Neg10$gene2

genelist1 <- c("Fth1","Anpep","Ptges","Hmox1","Ptgs2","Slpi","Ccl3","Mmp14","Nos2","Cd84",
               "Plac8","Ly6e", "Gbp8","Mndal","Tspan13","Ms4a4c","Gbp2b","Calhm6","Ly6c2" )

corr.point(corr = Macro.corr, genelist1)

ggsave("Macrophages biotin corr.eps")
ggsave("Macrophages biotin corr.png")




#--------------Analysis of Monocytes-------------#
#clustering
col <- brewer.pal(6, "Paired")
label1_Monocytes <- subset(label1,idents = 'Monocytes')
set.seed(1234)
label1_Monocytes <- FindNeighbors(label1_Monocytes, reduction = "pca", dims = 1:15)
label1_Monocytes <- FindClusters(label1_Monocytes,resolution = 0.6, algorithm = 1) 
label1_Monocytes <- RunTSNE(object=label1_Monocytes,dims.use=1:15, do.fast=TRUE,check_duplicates = FALSE)
label1_Monocytes <- RunUMAP(label1_Monocytes, reduction = "pca", dims = 1:15)
DimPlot(label1_Monocytes,label = TRUE,cols = col)  
ggsave(filename = "Monocytes clusters/label1_Monocytes clusters.png")
ggsave(filename = "Monocytes clusters/label1_Monocytes clusters.eps")

# Top genes expressed in Monocytes clusters
label1_Monocytes.markers <- FindAllMarkers(label1_Monocytes, 
                                           only.pos = TRUE,
                                           min.pct = 0.25, 
                                           logfc.threshold = 0.25)
top10 <- label1_Monocytes.markers %>% group_by(cluster) %>%  top_n(n = 10, wt = avg_log2FC)

write.csv(top10, "Monocytes clusters/label1_Monocytes top10markers.csv")

label1_Monocytes_top10markers <- read_csv("Monocytes clusters/label1_Monocytes top10markers.csv")

marker <- unique(label1_Monocytes_top10markers$gene)

marker1 <- marker %>% grep(pattern="^Rp", 
                           value = T, 
                           invert = T) %>% grep(pattern="mt", 
                                                value = T,
                                                invert = T)

DotPlot(label1_Monocytes, features = marker1)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))

ggsave("Monocytes clusters/dotplot Monocytes marker.png")
ggsave("Monocytes clusters/dotplot Monocytes marker.eps")



#--------Cellular landscapes over Biotag signal
#biotin low:L3, biotin medium:L2, biotin high:L1
srt_corr <- label1_Monocytes
ncol(srt_corr)

srt_corr@meta.data$level.biotin <- "NA"                     
srt_corr@meta.data$level.biotin[1:284] <-  "L3"                    
srt_corr@meta.data$level.biotin[285:568] <-  "L2"
srt_corr@meta.data$level.biotin[569:854] <-  "L1"


VlnPlot(srt_corr,features = "biotin",pt.size=0,group.by = "level.biotin")
ggsave(filename = "Monocytes biolevel.png")
ggsave(filename = "Monocytes biolevel.eps")

Idents(srt_corr) <- 'seurat_clusters'
table(srt_corr@meta.data[["level.biotin"]])
yy <- table(Idents(srt_corr),srt_corr@meta.data[["level.biotin"]])

col <- brewer.pal(9, "Paired")

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
ggsave(filename = "Monocytes.table.png")
ggsave(filename = "Monocytes.table.eps")


#------------Correlation analysis 

srt_corr <- label1_Monocytes 
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

write.csv(corr, "Monocytes_biotin correlation.csv")


Mono.corr <- corr
Mono.corr <- Mono.corr[order(Mono.corr$rho),]


corr.point <-  function(corr, label.gene, title="Monocytes biotin corr"){
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
Pos10 <- filter(corr, corr$rho >= 0.201)
Pos10$gene2
Neg10 <- filter(corr, corr$rho <= -0.095)
Neg10$gene2

genelist1 <- c("Thbs1","Clec4d","Mcemp1","Marcksl1","Ier3","Nfkbie","Cd14","Ccr1","Cxcl2","Clec4e",
               "Gm42418","Gpx1","Fos","Cd74","Ace","Ptma","Cmss1","Dusp1","Tab3","Vamp1")



corr.point(corr = Mono.corr, genelist1)
ggsave("Monocytes biotin corr.eps")



