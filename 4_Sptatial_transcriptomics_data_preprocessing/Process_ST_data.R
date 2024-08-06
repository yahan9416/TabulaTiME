#ST data process
library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
Oral_ST_P13Tumor=Load10X_Spatial("./ST_data/CRC_34417225/ST/ST/liver1/",filename="filtered_feature_bc_matrix.h5")
Oral_ST_P13Tumor <- PercentageFeatureSet(Oral_ST_P13Tumor, "^MT-", col.name = "percent_mito")
Oral_ST_P13Tumor = Oral_ST_P13Tumor[, Oral_ST_P13Tumor$nFeature_Spatial > 200 & Oral_ST_P13Tumor$percent_mito < 25 ]
Oral_ST_P13Tumor <- Oral_ST_P13Tumor[!grepl("^mt-", rownames(Oral_ST_P13Tumor)), ]
####quality control
pdf("Oral_ST_P13_after_QC.pdf",width = 12,height = 12)
SpatialFeaturePlot(Oral_ST_P13Tumor, features = c("nCount_Spatial", "nFeature_Spatial","percent_mito"))
dev.off()

##CopyKAT######
library(copykat)
expmat=GetAssayData(Oral_ST_P13Tumor)
expmat=as.matrix(expmat)
copykat.test <- copykat(rawmat=expmat, id.type="S", ngene.chr=1, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="", n.cores=10)
predict=copykat.test$prediction[,2]

####Clustering
expmat=GetAssayData(Oral_ST_P13Tumor)
Oral_ST_P13Tumor=ral_ST_P13Tumor
Oral_ST_P13Tumor <- RunPCA(Oral_ST_P13Tumor, assay = "SCT", verbose = FALSE)
Oral_ST_P13Tumor <- FindNeighbors(Oral_ST_P13Tumor, reduction = "pca", dims = 1:30)
Oral_ST_P13Tumor <- FindClusters(Oral_ST_P13Tumor, verbose = FALSE)
Oral_ST_P13Tumor <- RunUMAP(Oral_ST_P13Tumor, reduction = "pca", dims = 1:30)
pdf("Oral_ST_P13Tumor_Seurat_cluster.pdf",width = 12,height = 5)
DimPlot(Oral_ST_P13Tumor, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()


######Marker gene Visualization
DefaultAssay(Oral_ST_P13Tumor)<-"Spatial"
setwd("./ST_data/CRC_34417225/ST/ST/liver1")
pdf("marker_gene_expression_seurat.pdf",width = 15,height = 30)
SpatialFeaturePlot(object = Oral_ST_P13Tumor, features = c("COL1A1","CTHRC1","SFRP1","KRT19","PTPRC","CCL5","SAA1"), ncol = 1)
dev.off()