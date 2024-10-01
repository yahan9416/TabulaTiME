#' Identify the meta-program of each cell lineage, and process the meta-programs (MPs) with a focus on quality control and annotation.
#' @param jaccrad_similarity_path Path of jaccrad similarity of each pragrams
#' @param Top_Feature_path Path of top features of each programs
#' @param gmt_path Path of Hallmark gene list 
#' @author Ya Han


args <- commandArgs(trailingOnly = TRUE)
jaccrad_similarity_path <- args[1]
Top_Feature <-args[2]
gmt_path<-args[3]

#First load the simliarity of each NMP derived programs
jaccrad_similarity=readRDS(jaccrad_similarity_path)
temp_result=apply(jaccrad_similarity,1,function(x){length(which(x >= 10))})
temp_result_more10=names(temp_result)[which(temp_result > 4)]
jaccrad_similarity= jaccrad_similarity[match(temp_result_more10,rownames(jaccrad_similarity)),match(temp_result_more10,colnames(jaccrad_similarity))]



total_cells=rownames(jaccrad_similarity)
batch_identify_MPs<-function(temp_index){
  temp_result=apply(jaccrad_similarity,1,function(x){length(which(x > 10))})
  if(length(total_cells) >1){
    Potiential_founder=names(which.max(temp_result))
    temp_similarity=jaccrad_similarity[match(Potiential_founder,rownames(jaccrad_similarity)),]
    cluster_cells=c(Potiential_founder,names(temp_similarity)[which(temp_similarity >10)])
    total_cells<<-setdiff(total_cells,cluster_cells)
    jaccrad_similarity<<-jaccrad_similarity[match(total_cells,rownames(jaccrad_similarity)),match(total_cells,colnames(jaccrad_similarity))]
    if(length(cluster_cells) >5){
      return(cluster_cells)
    }
  }
}
result<-apply(matrix(1:500),1,batch_identify_MPs)
lapply(result,function(x) length(x))
result=result[1:length(which(unlist(lapply(result,function(x) length(x))) >0))]
Identified_MPs=result

#Process reuslt
Cell_used=unique(unlist(result))
jaccrad_similarity<<-jaccrad_similarity[match(Cell_used,rownames(jaccrad_similarity)),match(Cell_used,colnames(jaccrad_similarity))]

annotation_col=data.frame(Sample=as.vector(unlist(lapply(strsplit(colnames(jaccrad_similarity)," "),function(x) x[1] ))))
annotation_col$Cancertype=as.vector(unlist(lapply(strsplit(annotation_col$Sample,"_"),function(x) x[1] )))
annotation_col$Dataset=as.vector(unlist(lapply(strsplit(annotation_col$Sample,"_"),function(x) x[2] )))
annotation_col$Sample=as.vector(unlist(lapply(strsplit(annotation_col$Sample,"_"),function(x) x[4] )))

rownames(annotation_col)=rownames(jaccrad_similarity)
MPs<<-rep(1,nrow(annotation_col))
temp_result=apply(matrix(1:length(result)),1,function(x){MPs[match(result[[x]],rownames(annotation_col))]<<-x })
annotation_col$MPs=MPs
annotation_col$MPs=as.character(annotation_col$MPs)


jaccrad_similarity=jaccrad_similarity[match(unlist(result),rownames(jaccrad_similarity)),match(unlist(result),rownames(jaccrad_similarity))]


############MPs quality control###########

MPs_dataset=table(annotation_col$MPs,annotation_col$Dataset)
MPs_dataset_sta=apply(MPs_dataset,1,function(x) length(which(x >0)))
unidataset_MPs=names(MPs_dataset_sta)[which(MPs_dataset_sta ==1)]

MPs_cancertype=table(annotation_col$MPs,annotation_col$Cancertype)
MPs_cancertype_sta=apply(MPs_cancertype,1,function(x) length(which(x >0)))
unicancertype_MPs=names(MPs_cancertype_sta)[which(MPs_cancertype_sta ==1)] 

# with a strong enrichment of either ribosomal protein genes or mitochondrial-encoded genes 是low-qulaity data
#Module function enrichment analysis
#generated genes in each MPs
Top_Feature=readRDS(Top_Feature_path)
batch_generate_gene_each_MPs<-function(temp_mps){
  #print(temp_mps)
  temp_program=rownames(annotation_col)[which(annotation_col$MPs == temp_mps)]
  gene_list=unlist(Top_Feature[,match(temp_program,colnames(Top_Feature))])
  if(length(grep("\\.",gene_list)) !=0){
    gene_list=gene_list[-1*grep("\\.",gene_list)]}
  if(length(grep("^LINC",gene_list)) !=0){
    gene_list=gene_list[-1*grep("^LINC",gene_list)]}
  appear_num=sort(table(gene_list),decreasing = TRUE)[50]
  if(as.numeric(appear_num) == 1){
    appear_num=2}
  
  if(length(names(table(gene_list))[which(table(gene_list) >= appear_num)]) > 50){
    return(names(sort(table(gene_list),decreasing = TRUE)[1:50]))
  }
  
  return(names(table(gene_list))[which(table(gene_list) >= appear_num)])
}
MPs_genelist=apply(matrix(unique(annotation_col$MPs)),1,batch_generate_gene_each_MPs)

range(unlist(lapply(MPs_genelist,function(x) length(x))))
names(MPs_genelist)=unique(annotation_col$MPs)

#####Functional enrichment analysis###############
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(org.Hs.eg.db)
args = commandArgs(T)
genelist=as.vector(unlist(args[1]))
output_path = args[2] #setting fold change

geneset_name=c("h.all.v7.0.entrez.gmt","c1.all.v7.0.entrez.gmt","c2.all.v7.0.entrez.gmt","c3.all.v7.0.entrez.gmt","c4.all.v7.0.entrez.gmt","c5.all.v7.0.entrez.gmt","c6.all.v7.0.entrez.gmt","c7.all.v7.0.entrez.gmt")
description=c("Hallmark gene sets","Positional gene sets","Curated gene sets","Motif gene sets","Computational gene sets","GO gene sets","Oncogenic signatures","Immunologic signatures")
genesets=matrix(c(geneset_name,description),ncol=2,byrow=F)


batch_do_enrich<-function(x){
  print(x)
  #gene_set=DEG_gene$gene[which(DEG_gene$cluster == x )]
  #gene_set=rownames(clusterP_N.markers)
  gene_set=x
  gene_set <- bitr(gene_set, fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
  
  all_enrichment_result<<-NULL
  enrichment_analysis<-function(y){
    current_gene_set=read.gmt(paste0(gmt_path,genesets[y,1]))
    enrichment_result=enricher(gene_set$ENTREZID, TERM2GENE=current_gene_set)
    if ( dim(enrichment_result)[1] >0 ) {
      enrichment_result=cbind(as.data.frame(enrichment_result)[,c(1,3,5,7)],collection=genesets[y,2])
      all_enrichment_result<<-rbind(all_enrichment_result,enrichment_result)}
  }
  result=apply(matrix(c(1,3,6)),1,enrichment_analysis)
  
  return(all_enrichment_result)
}
finall_result<-lapply(MPs_genelist,batch_do_enrich)

names(finall_result)=names(MPs_genelist)


annotation_result<<-NULL
process_annotation_result<-function(x){
  annotation_result<<-rbind(annotation_result,as.matrix(cbind(finall_result[[x]],names(finall_result)[x])))
}
result=apply(matrix(1:length(finall_result)),1,process_annotation_result)
rownames(annotation_result)=NULL

#write.table(annotation_result,file="Myeloid_NMF_10_Module_enrich_result.txt",col.names = TRUE,row.names = FALSE,sep="\t",quote = FALSE)
#####Remove the low-qulaity MPs#######
names(Identified_MPs)=1:length(Identified_MPs)
result=result[-1*Lowquality_MPs]
jaccrad_similarity=jaccrad_similarity[match(unlist(result),rownames(jaccrad_similarity)),match(unlist(result),rownames(jaccrad_similarity))]

annotation_col=annotation_col[-1*which(annotation_col$MPs %in% Lowquality_MPs ),]
MPs_cluster_result=result
annotation_result=annotation_result[-1*which(annotation_result[,6] %in% Lowquality_MPs),]
