library(ggsankey)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

Immunity_Six=read.table("./Immunity_SixImmunesubtype.txt",header = TRUE,fill=TRUE,sep="\t")
Immunity_Six=Immunity_Six[which(!is.na(Immunity_Six$Immune.Subtype)),]
table(Immunity_Six$Immune.Subtype)

TabulaTIME_ecotype=readRDS("./cluster_information.rds")
TabulaTIME_ecotype$Patiene=substring(rownames(TabulaTIME_ecotype),1,12)

length(intersect(TabulaTIME_ecotype$Patiene,Immunity_Six[,1]))
coocurrence_sample=intersect(TabulaTIME_ecotype$Patiene,Immunity_Six[,1])
TabulaTIME_ecotype=TabulaTIME_ecotype[match(coocurrence_sample,TabulaTIME_ecotype$Patiene),]
Immunity_Six=Immunity_Six[match(coocurrence_sample,Immunity_Six[,1]),]
Ecotype_info=cbind(Immunity_Six$Immune.Subtype,TabulaTIME_ecotype$Ecotypes)
rownames(Ecotype_info)=coocurrence_sample
colnames(Ecotype_info)=c("Immunity_Six","TabulaTIME_ecotype")
Ecotype_info=as.data.frame(Ecotype_info)
#Ecotype_info$Immunity_Six=gsub("C","",Ecotype_info$Immunity_Six)
#Ecotype_info$TabulaTIME_ecotype=gsub("C_","",Ecotype_info$TabulaTIME_ecotype)


Ecotype_info <- Ecotype_info %>%
  make_long(Immunity_Six,TabulaTIME_ecotype)

p=ggplot(Ecotype_info, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  scale_fill_brewer(palette="Paired")+
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(plot.title = element_text(hjust = .5)) +
  ggtitle("Immunity & TabulaTIME ecotype")

ggsave(file="TabulaTIME_Immunity_Ecotype_Sankey_figure.pdf",p,width = 10,height = 9)


###Cell
Cell_ten=read.table("./Cell_Altas_TCGAsample_10CEs.txt",header = TRUE,fill=TRUE,sep="\t")
Cell_ten=Cell_ten[-1*which(Cell_ten$Discrete.assignment == "Unassigned"),]
Cell_ten$Sample=substring(Cell_ten$Sample,1,15)

TabulaTIME_ecotype=readRDS("./cluster_information.rds")
#TabulaTIME_ecotype$Patiene=substring(rownames(TabulaTIME_ecotype),1,12)

length(intersect(rownames(TabulaTIME_ecotype),Cell_ten$Sample))
coocurrence_sample=intersect(rownames(TabulaTIME_ecotype),Cell_ten$Sample)
TabulaTIME_ecotype=TabulaTIME_ecotype[match(coocurrence_sample,rownames(TabulaTIME_ecotype)),]
Cell_ten=Cell_ten[match(coocurrence_sample,Cell_ten$Sample),]
Ecotype_info=cbind(Cell_ten$Discrete.assignment,TabulaTIME_ecotype$Ecotypes)
rownames(Ecotype_info)=coocurrence_sample
colnames(Ecotype_info)=c("Cell_ten","TabulaTIME_ecotype")
Ecotype_info=as.data.frame(Ecotype_info)
table(Ecotype_info$TabulaTIME_ecotype,Ecotype_info$Cell_ten)

library(ggsankey)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
Ecotype_info <- Ecotype_info %>%
  make_long(Cell_ten,TabulaTIME_ecotype)

color_man=c( "#5D69B1", "#52BCA3", "#99C945", "#CC61B0", "#24796C","#DAA51B", "#2F8AC4", "#764E9F", "#ED645A", "#A5AA99", "#BCBD22", "#B279A2", "#EECA3B", "#17BECF", "#FF9DA6", "#778AAE", "#1B9E77","#A6761D", "#526A83", "#B82E2E", "#80B1D3", "#68855C", "#D95F02","#BEBADA", "#AF6458", "#D9AF6B", "#9C9C5E", "#625377", "#8C785D","#88CCEE", "#E73F74", "#FFFFB3", "#CCEBC5", "#332288", "#A65628","#0096FF", "#F3D4F4", "#FDCDAC", "#548235", "#9271CB", "#917F1D")
p=ggplot(Ecotype_info, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  scale_fill_manual(values=color_man)+
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(plot.title = element_text(hjust = .5)) +
  ggtitle("Immunity & TabulaTIME ecotype")
ggsave(file="TabulaTIME_Cell_Ecotype_Sankey_figure.pdf",p,width = 15,height = 15)


####Cancer cell#
Immunity_Six=read.table("./Cancercell_TMEsubtypes.txt",header = TRUE,fill=TRUE,sep="\t")
Immunity_Six=Immunity_Six[which(!is.na(Immunity_Six[,8])),]


TabulaTIME_ecotype=readRDS("./cluster_information.rds")
TabulaTIME_ecotype$Patiene=substring(rownames(TabulaTIME_ecotype),1,12)

length(intersect(TabulaTIME_ecotype$Patiene,Immunity_Six[,1]))
coocurrence_sample=intersect(TabulaTIME_ecotype$Patiene,Immunity_Six[,1])
TabulaTIME_ecotype=TabulaTIME_ecotype[match(coocurrence_sample,TabulaTIME_ecotype$Patiene),]
Immunity_Six=Immunity_Six[match(coocurrence_sample,Immunity_Six[,1]),]
Ecotype_info=cbind(Immunity_Six[,8],TabulaTIME_ecotype$Ecotypes)
rownames(Ecotype_info)=coocurrence_sample
colnames(Ecotype_info)=c("CancerCell_TMEsubtype","TabulaTIME_ecotype")
Ecotype_info=as.data.frame(Ecotype_info)
#Ecotype_info$Immunity_Six=gsub("C","",Ecotype_info$Immunity_Six)
#Ecotype_info$TabulaTIME_ecotype=gsub("C_","",Ecotype_info$TabulaTIME_ecotype)
Ecotype_info=as.data.frame(Ecotype_info)
Ecotype_info$CancerCell_TMEsubtype[which(Ecotype_info$CancerCell_TMEsubtype == "D")]="1_D"
Ecotype_info$CancerCell_TMEsubtype[which(Ecotype_info$CancerCell_TMEsubtype == "F")]="2_F"
Ecotype_info$CancerCell_TMEsubtype[which(Ecotype_info$CancerCell_TMEsubtype == "IE/F")]="3_IE/F"
Ecotype_info$CancerCell_TMEsubtype[which(Ecotype_info$CancerCell_TMEsubtype == "IE")]="4_IE"
Ecotype_info$TabulaTIME_ecotype[which(Ecotype_info$TabulaTIME_ecotype == "C_DesHP")]="1_Des"
Ecotype_info$TabulaTIME_ecotype[which(Ecotype_info$TabulaTIME_ecotype == "C_DesLP")]="1_Des"
Ecotype_info$TabulaTIME_ecotype[which(Ecotype_info$TabulaTIME_ecotype == "C_Naive")]="2_Naive"
Ecotype_info$TabulaTIME_ecotype[which(Ecotype_info$TabulaTIME_ecotype == "C_ActHigSt")]="3_ActHigSt"
Ecotype_info$TabulaTIME_ecotype[which(Ecotype_info$TabulaTIME_ecotype == "C_ActLowSt")]="4_ActLowSt"

library(ggsankey)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
Ecotype_info <- Ecotype_info %>%
  make_long(CancerCell_TMEsubtype,TabulaTIME_ecotype)

p=ggplot(Ecotype_info, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = .6) +
  geom_alluvial_text(size = 3, color = "white") +
  scale_fill_brewer(palette="Paired")+
  theme_alluvial(base_size = 18) +
  labs(x = NULL) +
  theme(plot.title = element_text(hjust = .5)) +
  ggtitle("CancerCell & TabulaTIME ecotype")

ggsave(file="TabulaTIME_Cancercell_Ecotype_Sankey_figure_cor.pdf",p,width = 10,height = 9)

