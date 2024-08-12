
library(copykat)
library(Seurat)
library(MAESTRO)

##CopyKAT######
args <- commandArgs(trailingOnly = TRUE)
ST_seurat_path<- args[1]

ST_seurat=readRDS(ST_seurat_path)
expmat=GetAssayData(ST_seurat)
expmat=as.matrix(expmat)
copykat.test <- copykat(rawmat=expmat, id.type="S", ngene.chr=1, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=10)
predict=copykat.test$prediction[,2]