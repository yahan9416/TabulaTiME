#Trajectory
library(monocle)
library(MAESTRO)
library(Seurat)
library(ggplot2)

Mye_Fibro_Endo_Seu #load the seurat obhect
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


setwd("/Storage/hanya/tme_Difference/Celltype_transition")

cds <- estimateSizeFactors(cds)
cds=estimateDispersions(cds)
#cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds)
cds <- orderCells(cds)

p=plot_cell_trajectory(cds, color_by = "Celltype",cell_size=0.5,show_tree = F,show_backbone=F,show_branch_points=F)+scale_color_manual(values=c("#726BAE",rep("grey",10),"#EA945A",rep("grey",5),"#276D9F","grey","grey","grey","grey","#60A897","#F5BC6E","grey","grey","grey","#86A667","#8AB6D6"))
ggsave("Monocel2_CelltypeDEG_trajectory_cell_type.pdf",p,width = 8,height = 6)

p=plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size=0.5,show_tree = F,show_backbone=F,show_branch_points=F)
ggsave("Monocel2_CelltypeDEG_trajectory_pseudotime_rev.pdf",p,width = 5,height = 4)