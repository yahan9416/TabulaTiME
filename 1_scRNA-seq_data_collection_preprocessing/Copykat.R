#Running copykat
#The input data of copykat is raw count 

#conda activate ./miniconda3/envs/MAESTRO


library(copykat)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)
rowcount_matrix <- args[1]

expmat=Read10X_h5(rowcount_matrix)#the path of raw count expression
copykat.test <- copykat(rawmat=expmat, id.type="S", ngene.chr=1, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=10)
predict=copykat.test$prediction[,2]