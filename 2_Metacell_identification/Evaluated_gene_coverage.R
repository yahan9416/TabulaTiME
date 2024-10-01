#' Evaluating the number of genes detected in the dataset.
#' @param Metacell_info The RData file of each dataset corresponding metacell information
#' @author Ya Han

args <- commandArgs(trailingOnly = TRUE)
Metacell_info<- args[1]

load(Metacell_info)
# inclusing the expression of metacell Metacell_TPM
test_gene_coverage<-function(sample_gene_exp){
     return(length(which(sample_gene_exp != 0)))
}
gene_coverage=apply(Metacell_TPM,2,test_gene_coverage)
Coverage=data.frame(coverage=gene_coverage,cell_number=unlist(strsplit(temp_file,"\\/"))[10])

  