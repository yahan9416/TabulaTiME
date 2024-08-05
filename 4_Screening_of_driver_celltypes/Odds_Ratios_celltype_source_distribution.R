meta_info$Treatment=Before_meta$Treatment[match(rownames(meta_info),rownames(Before_meta))]
treatment_meta=meta_info[which(meta_info$Treatment == "None"),]

celltype_source_sta=table(treatment_meta$Source,treatment_meta$curated_anno)
celltype_source_Fish<<-NULL
batch_character_meta_cluster_dist<-function(temp_celltype){
  batch_process_source<-function(temp_source){
    celltype_index=match(temp_celltype,colnames(celltype_source_sta))
    source_index=match(temp_source,rownames(celltype_source_sta))
    a=celltype_source_sta[source_index,celltype_index]
    b=rowSums(celltype_source_sta)[match(temp_source,names(rowSums(celltype_source_sta)))]-a
    c=colSums(celltype_source_sta)[match(temp_celltype,names(colSums(celltype_source_sta)))]-a
    d=sum(celltype_source_sta)-a-b-c
    temp_source_cell <- data.frame(
      "CTHRC1" = c(a, c),
      "non_CTHRC1" = c(b, d),
      row.names = c("Tumor", "Non-tumor"),
      stringsAsFactors = FALSE )
    test <- fisher.test(temp_source_cell)
    #return(temp_celltype,temp_source,unlist(test)[4],unlist(test)[1])
    celltype_source_Fish<<-rbind(celltype_source_Fish,c(temp_celltype,temp_source,unlist(test)[4],unlist(test)[1]))
  }
  result=apply(matrix(c("Normal","Precancerous","Tumor","Metastatic" )),1,batch_process_source)
}
result=apply(matrix(unique(meta_info$curated_anno)),1,batch_character_meta_cluster_dist)

celltype_source_Fish=as.data.frame(celltype_source_Fish)
colnames(celltype_source_Fish)=c("Celltype","Source","OR","P_value")
celltype_source_Fish$OR=as.numeric(celltype_source_Fish$OR)
celltype_source_Fish$P_value=as.numeric(celltype_source_Fish$P_value)
celltype_source_Fish$P_adjust=p.adjust(celltype_source_Fish$P_value,method="BH") 

heatmap_annotation<<-matrix(rep(0,length(celltype_source_Fish$OR)),ncol=length(unique(celltype_source_Fish$Source)))
celltype_source_Fish_heat=matrix(unlist(celltype_source_Fish$OR),byrow = TRUE,ncol=length(unique(celltype_source_Fish$Source)))
colnames(celltype_source_Fish_heat)=c("Blood","Normal","Precancerous","Tumor","Metastatic" )
rownames(celltype_source_Fish_heat)=unique(meta_info$curated_anno)

p_value=matrix(unlist(celltype_source_Fish$P_adjust),byrow = TRUE,ncol=length(unique(celltype_source_Fish$Source)))
colnames(p_value)=colnames(celltype_source_Fish_heat)
rownames(p_value)=rownames(celltype_source_Fish_heat)

library(ComplexHeatmap)
library(pheatmap)
library(circlize)

if (!is.null(p_value)){
  ssmt <- p_value< 0.01
  p_value[ssmt] <-'**'
  smt <- p_value >0.01& p_value <0.05
  p_value[smt] <- '*'
  p_value[!ssmt&!smt]<- ''
} else {
  p_value <- F
}

pdf("Celltype_source_Fish_heatmap_treatmentnaive.pdf",height = 7,width = 5)
Heatmap(celltype_source_Fish_heat,cluster_columns=FALSE,c( colorRampPalette(colors = c("white","#FAAC6E","#E35726"))(30)), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text( p_value[i, j], x, y, gp = gpar(fontsize = 17, col = "white"))
        })
dev.off()