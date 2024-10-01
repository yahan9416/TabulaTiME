#' Trajectory analysis
#' @param Seurat_obj_path Path of seurat obj
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Seurat_obj_path <- args[1]

library(monocle)
library(MAESTRO)
library(Seurat)
library(ggplot2)

Mye_Fibro_Endo_Seu=readRDS(Seurat_obj_path) #load the seurat obhect
SeuratObj.markers <- FindAllMarkers(Mye_Fibro_Endo_Seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DefaultAssay(Mye_Fibro_Endo_Seu)<-"integrated"
expmat=GetAssayData(Mye_Fibro_Endo_Seu)
cell_metadata=Mye_Fibro_Endo_Seu@meta.data
expmat=as.matrix(expmat)

Cell_type_DEG=SeuratObj.markers[which(SeuratObj.markers$avg_log2FC > 0.5),]
marker_gene=intersect(unique(Cell_type_DEG$gene),rownames(expmat))
index_marker=match(marker_gene,rownames(expmat))
expmat=expmat[index_marker,]
gene_metadata=rownames(expmat)
gene_metadata=as.data.frame(gene_metadata)
rownames(gene_metadata)=rownames(expmat)
colnames(gene_metadata)="gene_short_name"
colnames(expmat)=rownames(cell_metadata)

pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_metadata)
cds <- newCellDataSet(expmat,
                      phenoData = pd,
                      featureData = fd,
                      expressionFamily=negbinomial.size())




cds <- estimateSizeFactors(cds)
cds=estimateDispersions(cds)
#cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds)
cds <- orderCells(cds)
