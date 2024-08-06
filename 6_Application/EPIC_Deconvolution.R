#Deconvolution


library(EPIC)
microarray=read.table("./data_mrna_agilent_microarray.txt",header = TRUE,sep="\t",row.names = 1)
res3 <- EPIC(bulk=microarray, reference=TRef)
names(res3)
head(res3$cellFractions)
head(res3$cellFractions$CD8_Tcells)
saveRDS(res3,file="./Estimate_Imme/EPIC_microarray_estimate_result.rds")
