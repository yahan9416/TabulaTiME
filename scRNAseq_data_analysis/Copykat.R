#Running copykat
#The input data of copykat is raw count 

#conda activate /mnt/Storage/home/yuejiali/miniconda3/envs/MAESTRO


library(copykat)
library(Seurat)

setwd("/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Without_P12OLP/CCA_result/Nfeatures5000_includingInteract_CCA/Seperate_cell_lineage/Epithelial")
Epi_Seurat=readRDS("CCA_Epithelial_repca_cluster.rds")
#DefaultAssay(Epi_Seurat)<-"integrated"
expmat=as.matrix(Epi_Seurat@assays$RNA@counts)
copykat.test <- copykat(rawmat=expmat, id.type="S", ngene.chr=1, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=10)

predict=copykat.test$prediction[,2]