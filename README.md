<a href=""><img src="" alt="DOI"></a> <a href=""><img src="http://timer2.compbio.cn/TabulaTIME" alt="Website"></a> 
# Workflow for integrating tumor scRNA-seq data and constructing pan-cancer landscapes: TabulaTIME


Tumor microenvironment (TME) evolves during malignant cell recruiting and reprograming the non-cancerous cells and remodeling the vasculature and extracellular matrix (ECM) to orchestrate a tumor-supportive environment. The composition and functional state of TME can vary considerably between the organ in which the tumor arises, due to the unique tissue-resident cells. Here, we collected the expression of 4,254,586 cells from 735 patient samples across 36 cancer types, and construct a pan-cancer landscape to depict the diversity of TMEs.
## Environment 
    Ubuntu 9.3.0
    R version 4.0.5	
    Python version 3.8.10	

## Install software
### Install R package MAESTRO V1.4.1
    conda config --add channels defaults
    conda config --add channels liulab-dfci
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda install mamba -c conda-forge
    mamba create -n MAESTRO maestro=1.4.1 -c liulab-dfci
### Install R package Seurat v2.3.4 	
    source("https://z.umn.edu/archived-seurat")
### Install R package Monocle v2.8 	
    source("http://bioconductor.org/biocLite.R") 
    biocLite("monocle")	

