#Rogue score estimate the purity of each cluster or cell types

library(ROGUE)
library("dplyr")
#Calculate the ROGUE value of each putative cluster for each sample.
args <- commandArgs(trailingOnly = TRUE)
Seurat_obj <- args[1]


SeuratObj=readRDS(Seurat_obj)
DefaultAssay(SeuratObj)<-"integrated"
expmat=GetAssayData(SeuratObj)

rogue.res2 <- rogue(expmat, labels = SeuratObj@meta.data$curated_anno, samples = SeuratObj@meta.data$Sample, platform = "UMI", span = 0.6)


#Process the result
Rogue_score<<-NULL
process_result<-function(temp_celltype){
  temp_data=cbind(names(rogue.res2)[temp_celltype],na.omit(rogue.res2[[temp_celltype]]))
  Rogue_score<<-rbind(Rogue_score,temp_data)
}
result=apply(matrix(1:length(rogue.res2)),1,process_result)

Rogue_score=as.data.frame(Rogue_score)
colnames(Rogue_score)<-c("Celltype","Rogue_score")
Rogue_score$Rogue_score=as.numeric(as.vector(Rogue_score$Rogue_score))

library(ggplot2)
p=ggplot(Rogue_score, aes(x=Celltype, y=Rogue_score,color=Celltype)) + 
  geom_boxplot(fill="white")+
  theme_classic()+theme(legend.position='none',axis.text.x = element_text(angle = 90))+scale_colour_manual(values=color_man)
ggsave(file="CD8TNK_Inte_Rogue.pdf",p,width = 5,height = 4)
saveRDS(Rogue_score,file="Rogue_score_celltype_sample.rds")
