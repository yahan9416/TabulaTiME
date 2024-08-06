#Load the estimated immune infiltration for TCGA patients
TCGA_estimate_Infi=read.table("./infiltration_estimation_for_tcga.csv",header = TRUE,sep=",")
TCGA_estimate_Infi[,1]=substring(TCGA_estimate_Infi[,1],1,12)

TabulaTIME_ecotype=readRDS("./cluster_information.rds")
TabulaTIME_ecotype$Patiene=substring(rownames(TabulaTIME_ecotype),1,12)
Tabul_cosample=intersect(TabulaTIME_ecotype$Patiene,TCGA_estimate_Infi[,1])
TabulaTIME_ecotype=TabulaTIME_ecotype[match(Tabul_cosample,TabulaTIME_ecotype$Patiene),]


Cancer_cell=read.table("./Cancercell_TMEsubtypes.txt",header = TRUE,fill=TRUE,sep="\t")
Cancer_cell=Cancer_cell[which(!is.na(Cancer_cell[,8])),]
cancercell_cosample=intersect(Cancer_cell[,1],TCGA_estimate_Infi[,1])
Cancer_cell=Cancer_cell[match(cancercell_cosample,Cancer_cell[,1]),]

Cell_ten=read.table("./Cell_Altas_TCGAsample_10CEs.txt",header = TRUE,fill=TRUE,sep="\t")
Cell_ten=Cell_ten[-1*which(Cell_ten$Discrete.assignment == "Unassigned"),]
Cell_ten$Sample=substring(Cell_ten$Sample,1,12)
Cell_cosample=intersect(Cell_ten[,1],TCGA_estimate_Infi[,1])
Cell_ten=Cell_ten[match(Cell_cosample,Cell_ten[,1]),]

Immunity_Six=read.table("./Immunity_SixImmunesubtype.txt",header = TRUE,fill=TRUE,sep="\t")
Immunity_Six=Immunity_Six[which(!is.na(Immunity_Six$Immune.Subtype)),]
Immunity_cosample=intersect(Immunity_Six[,1],TCGA_estimate_Infi[,1])
Immunity_Six=Immunity_Six[match(cancercell_cosample,Immunity_Six[,1]),]



batch_Immuneinfiltration_method<-function(temp_IEM){
  #TabulaTiME
  Temp_TCGA_infil=data.frame(Sample=TCGA_estimate_Infi[,1],Infil=rowSums(TCGA_estimate_Infi[,grep(temp_IEM,colnames(TCGA_estimate_Infi))]))
  Temp_TCGA_infil=Temp_TCGA_infil[match(Tabul_cosample,Temp_TCGA_infil[,1]),]
  Icutoff=median(na.omit(Temp_TCGA_infil$Infil))
  Temp_TCGA_infil$Infiltration_class="high"
  Temp_TCGA_infil$Infiltration_class[which(Temp_TCGA_infil$Infil < Icutoff)]="low"
  Temp_TCGA_infil$Tabluaclass=TabulaTIME_ecotype$Ecotypes[match(Temp_TCGA_infil$Sample,TabulaTIME_ecotype$Patiene)]
  Temp_TCGA_infil$TabluaBioclass="high"
  Temp_TCGA_infil$TabluaBioclass[which(Temp_TCGA_infil$Tabluaclass %in% c("C_DesHP","C_DesLP"))]="low"
  static=table(Temp_TCGA_infil$Infiltration_class,Temp_TCGA_infil$TabluaBioclass)
  Tabula_acc=(static[1,1]+static[2,2])/dim(Temp_TCGA_infil)[1]
  
  #Cancer cell
  Temp_TCGA_infil=data.frame(Sample=TCGA_estimate_Infi[,1],Infil=rowSums(TCGA_estimate_Infi[,grep(temp_IEM,colnames(TCGA_estimate_Infi))]))
  Temp_TCGA_infil=Temp_TCGA_infil[match(cancercell_cosample,Temp_TCGA_infil[,1]),]
  Icutoff=median(na.omit(Temp_TCGA_infil$Infil))
  Temp_TCGA_infil$Infiltration_class="high"
  Temp_TCGA_infil$Infiltration_class[which(Temp_TCGA_infil$Infil < Icutoff)]="low"
  Temp_TCGA_infil$Cancer_cellclass=Cancer_cell[match(Temp_TCGA_infil$Sample,Cancer_cell[,1]),8]
  Temp_TCGA_infil$Cancer_Biocellclass="high"
  Temp_TCGA_infil$Cancer_Biocellclass[which(Temp_TCGA_infil$Cancer_cellclass %in% c("D"))]="low"
  static=table(Temp_TCGA_infil$Infiltration_class,Temp_TCGA_infil$Cancer_Biocellclass)
  Cancercell_acc=(static[1,1]+static[2,2])/dim(Temp_TCGA_infil)[1]
  
  #Immunity
  Temp_TCGA_infil=data.frame(Sample=TCGA_estimate_Infi[,1],Infil=rowSums(TCGA_estimate_Infi[,grep(temp_IEM,colnames(TCGA_estimate_Infi))]))
  Temp_TCGA_infil=Temp_TCGA_infil[match(Immunity_cosample,Temp_TCGA_infil[,1]),]
  Icutoff=median(na.omit(Temp_TCGA_infil$Infil))
  Temp_TCGA_infil$Infiltration_class="high"
  Temp_TCGA_infil$Infiltration_class[which(Temp_TCGA_infil$Infil < Icutoff)]="low"
  Temp_TCGA_infil$Immunity_class=Immunity_Six[match(Temp_TCGA_infil$Sample,Immunity_Six[,1]),3]
  Temp_TCGA_infil$Immunity_Biocellclass="high"
  Temp_TCGA_infil$Immunity_Biocellclass[which(Temp_TCGA_infil$Immunity_class %in% c("C3","C4","C5 "))]="low"
  static=table(Temp_TCGA_infil$Infiltration_class,Temp_TCGA_infil$Immunity_Biocellclass)
  Immunity_acc=(static[1,1]+static[2,2])/dim(Temp_TCGA_infil)[1]
  
  
  #Cell
  Temp_TCGA_infil=data.frame(Sample=TCGA_estimate_Infi[,1],Infil=rowSums(TCGA_estimate_Infi[,grep(temp_IEM,colnames(TCGA_estimate_Infi))]))
  Temp_TCGA_infil=Temp_TCGA_infil[match(Cell_cosample,Temp_TCGA_infil[,1]),]
  Icutoff=median(na.omit(Temp_TCGA_infil$Infil))
  Temp_TCGA_infil$Infiltration_class="high"
  Temp_TCGA_infil$Infiltration_class[which(Temp_TCGA_infil$Infil < Icutoff)]="low"
  Temp_TCGA_infil$Cell_class=Cell_ten$Discrete.assignment[match(Temp_TCGA_infil$Sample,Cell_ten[,1])]
  Temp_TCGA_infil$Cell_Bioclass="high"
  Temp_TCGA_infil$Cell_Bioclass[which(Temp_TCGA_infil$Cell_class %in% c("CE7","CE8","CE4","CE5","CE6","CE2"))]="low"
  static=table(Temp_TCGA_infil$Infiltration_class,Temp_TCGA_infil$Cell_Bioclass)
  Cell_acc=(static[1,1]+static[2,2])/dim(Temp_TCGA_infil)[1]
  
  
  return(c(temp_IEM,Tabula_acc,Cancercell_acc,Immunity_acc,Cell_acc))
}
result=apply(matrix(matrix(c("TIMER","CIBERSORT.ABS","MCPCOUNTER","XCELL"))),1,batch_Immuneinfiltration_method)


data=data.frame(Accuracy=as.numeric(as.vector(unlist(result[-1,]))),Infiltration_method=rep(c("TIMER","CIBERSORT.ABS","MCPCOUNTER","XCELL"),each=4),Class_method=rep(c("TabulaTiMe","CancerCell","Immunity","Cell"),4))

library(ggplot2)
library(RColorBrewer)
p <- ggplot(data, aes(x=Infiltration_method, y=Accuracy,fill=Class_method)) + 
  geom_dotplot(binaxis='y', stackdir='center')+ scale_fill_brewer(palette="Dark2")+theme_linedraw()

ggsave(file="Four_Classmehtod_TCGA_immune_inflation.pdf",p,width = 5,height = 4)
