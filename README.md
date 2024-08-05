<a href=""><img src="" alt="DOI"></a> <a href=""><img src="http://timer2.compbio.cn/TabulaTIME" alt="Website"></a> 
# Workflow for integrating tumor scRNA-seq data and constructing pan-cancer landscapes: TabulaTIME

## 0 scRNA-seq data collection and preprocessing

### 0.1 Data collection

We collected the published cancer-associated scRNA-Seq datasets from 735 patients across 36 cancer types. These datasets were sourced from the GEO database (https://www.ncbi.nlm.nih.gov/geo/) and our previous work the Tumor Immune Single-cell Hub (TISCH) database (http://tisch.comp-genomics.org/home/). Additionally, we incorporated scRNA-Seq datasets derived from healthy donors, including 3 peripheral blood mononuclear cells (PBMC) and 6 datasets from normal tissues.

A total of 103 studies encompassing 4,479,563 cells were collated.

### 0.2 Data preprocessing

A standardized analysis workflow based on MAESTRO v1.1.0 for processing all the collected datasets, including quality control, batch effect removal, and cell clustering. The raw count or TPM (if available) served as input for the workflow.  

 To ascertain the detailed process of single-cell RNA sequencing datasets, we kindly direct you to the comprehensive instructions available within the following GitHub repository: https://github.com/DongqingSun96/TISCH/tree/master/code.

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

