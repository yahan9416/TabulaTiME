#Calculating the signature score for each cell type, we could well describe the function of different cell types.
#In this function, all signarute gene are defined in your dataset.


library(MAESTRO)
library(Seurat)
library(ggplot2)
library(ggpubr)
SeuratObj
cytotoxicity_fea=c("PRF1","IFNG","GNLY","NKG7","GZMB","GZMA","GZMH","KLRK1","KLRB1","KLRD1","CTSW","CST7","CCL4","CCL3")
Exhausteed_fea=c('PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'CD244', 'CD160', 'TIGIT')
Regulatory_fea=c("IL2RA","IL1R2","FOXP3","LAYN","TNFRSF9","FANK1","RTKN2","IL1R1","CUL9","IKZF2")
SeuratObj <- AddModuleScore(
  object = SeuratObj,
  features = list(cytotoxicity_fea),
  assay="RNA",
  ctrl = 5,
  name = 'Cytotoxicity'
)

SeuratObj <- AddModuleScore(
  object = SeuratObj,
  features = list(Exhausteed_fea),
  assay="RNA",
  ctrl = 5,
  name = 'Exhausted'
)

SeuratObj <- AddModuleScore(
  object = SeuratObj,
  features = list(cytotoxicity_fea),
  assay="RNA",
  ctrl = 5,
  name = 'Regulatory'
)
#然后将两个feature 当成x,y 轴画出Point 图。
#或者我们可以尝试按照Source 分开来看。

temp_met=SeuratObj@meta.data
data1=aggregate(Exhausted1~curated_anno,temp_met,mean)
data2=aggregate(Cytotoxicity1~curated_anno,temp_met,mean)
data=cbind(data1,data2)
data=data[,c(1,2,4)]
colnames(data)=c("Celltype","Exhausted_Score","Cytotoxic_Score")
data$Celltype=as.vector(unlist(data$Celltype))
data$Celltype[which(data$Celltype == "CD8_GZMK_exh")]="CD8_PDCD1_exh"
data$Celltype[which(data$Celltype == "CD8_GZMK_em")]="CD8_GZMK_mem"
data$Celltype=as.factor(data$Celltype)
library(ggrepel)

setwd("/fs/home/hanya/Project/TME_Immune_difference/TISCH_data_10X_2105/Lineage_seperate_analysis/CD8TNK/Basic_analysis/")
color_man=c("#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C","#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#A5AA99", "#BCBD22", "#B279A2", "#EECA3B", "#17BECF", "#FF9DA6", "#778AAE", "#1B9E77","#A6761D", "#526A83", "#B82E2E", "#80B1D3", "#68855C", "#D95F02","#BEBADA", "#AF6458", "#D9AF6B", "#9C9C5E", "#625377", "#8C785D","#88CCEE", "#E73F74", "#FFFFB3", "#CCEBC5", "#332288", "#A65628","#0096FF", "#F3D4F4", "#FDCDAC", "#548235", "#9271CB", "#917F1D")
library(ggrepel)
p=ggplot(data, aes(x=Cytotoxic_Score, y=Exhausted_Score,fill=Celltype,color=Celltype)) +
  geom_point(alpha=0.7) + theme_bw()+
  geom_text_repel(aes(label = data$Celltype),
                  size = 1.5) +scale_fill_manual(values=color_man)+scale_color_manual(values = color_man)+theme(legend.position="none")
ggsave(file="Celltype_Cyto_Exha_score_Mean1.pdf",p,width = 3.2,height = 3)
