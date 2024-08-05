#Running InferCNV
##Creating an InferCNV object based on your three required inputs: the read count matrix, cell type annotations, and the gene ordering file:


library(Seurat)
library(infercnv)
setwd("/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/InferCNV")
expr = Read10X_h5("/mnt/Storage2/home/hanya/project/Carcinogenesis/scRNAseq/Carcinogenesis_merge_gene_count.h5")

cell.info <-  data.frame(row.names = colnames(expr))
cell.info$anno <-"cells"
cell.info$anno[grep("Normal",colnames(expr))]="Normal"
colnames(cell.info)=NULL
write.table(cell.info,file = "cell.anno.txt",col.names = FALSE,sep = "\t",quote = FALSE)
##gene order file
gene_order_file=read.table("/mnt/Storage2/home/hanya/project/Dream_challenge_immune_infiltration/scRNA_data_process/inferCNV_gene_order_file.txt",header=F,sep="\t")
gene_order_file=gene_order_file[-1*which(duplicated(gene_order_file[,1])),]
rownames(gene_order_file)=gene_order_file[,1]
gene_order_file=gene_order_file[,-1]
# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= expr,
                                    annotations_file= cell.info,
                                    gene_order_file= gene_order_file,
                                    ref_group_names="Normal")


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,  # use 1 for smart-seq/Fluidigm C1, 0.1 for 10x-genomics/SNRS/MARS-Seq/Microwell/Drop-seq
                             out_dir="k2.c1",  # dir is auto-created for storing outputs
                             cluster_by_groups=FALSE,   # cluster
                             k_obs_groups = 2,
                             denoise=T)



#after run InferCNV process the result
infercnv.group <- read.csv("./k2.c1/infercnv.observation_groupings.txt",sep = " ",row.names = 1)
infercnv.group <- data.frame(row.names = rownames(infercnv.group),infercnv=infercnv.group[,1])
group.1 <- rownames(infercnv.group)[which(infercnv.group$infercnv == 1)]
group.2 <- rownames(infercnv.group)[which(infercnv.group$infercnv == 2)]

##separate the normal and cancer cells
infercnv_obj = readRDS('./k2.c1/run.final.infercnv_obj')
infercnv = as.data.frame(t(infercnv_obj@expr.data))
var.infercnv <- as.data.frame(apply(infercnv,1,function(x)var(x)))
var.cnv.group1 <- var.infercnv[which(row.names(var.infercnv) %in% group.1),]
var.cnv.group2 <- var.infercnv[which(row.names(var.infercnv) %in% group.2),]
if (mean(var.cnv.group1)>mean(var.cnv.group2)){
  tmp=data.frame(sample=group.1,anno=rep("cancer",length(group.1))) 
  tmp2=data.frame(sample=group.2,anno=rep("normal",length(group.2))) 
  infercnv.anno=rbind(tmp,tmp2)
}else{
  tmp=data.frame(sample=group.1,anno=rep("normal",length(group.1))) 
  tmp2=data.frame(sample=group.2,anno=rep("cancer",length(group.2))) 
  infercnv.anno=rbind(tmp,tmp2)
}
write.table(infercnv.anno,"infercnv.anno.txt",col.names = F,row.names = FALSE,sep = "\t",quote = FALSE)
