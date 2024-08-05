<a href=""><img src="" alt="DOI"></a> <a href=""><img src="http://timer2.compbio.cn/TabulaTIME" alt="Website"></a> 
# Workflow for integrating tumor scRNA-seq data and constructing pan-cancer landscapes: TabulaTIME

## 1. scRNA-seq data collection and preprocessing
The increasing accumulation of scRNA-seq datasets in the public domain allows for the integration of tumor-associated datasets from different cancer types and stages. This integration helps to identify common characteristics of the TME and detect pan-cancer or organ-specific mechanisms. 
### 1.1 Data collection

We collected the published cancer-associated scRNA-Seq datasets from 735 patients across 36 cancer types. These datasets were sourced from the GEO database (https://www.ncbi.nlm.nih.gov/geo/) and our previous work the Tumor Immune Single-cell Hub (TISCH) database (http://tisch.comp-genomics.org/home/). Additionally, we incorporated scRNA-Seq datasets derived from healthy donors, including 3 peripheral blood mononuclear cells (PBMC) and 6 datasets from normal tissues.

A total of 103 studies encompassing 4,479,563 cells were collated.

### 1.2 Data preprocessing

A standardized analysis workflow based on MAESTRO v1.1.0 for processing all the collected datasets, including quality control, batch effect evaluation and correction, and cell clustering. The raw count or TPM (if available) served as input for the workflow.  

To ascertain the detailed process of single-cell RNA sequencing datasets, we kindly direct you to the comprehensive instructions available within the following GitHub repository: https://github.com/DongqingSun96/TISCH/tree/master/code.

## 2. MetaCell identification

To reduce the technical noises and computing resource costs, TabulaTIME grouped cells with similar expressions into metacells within each dataset. The average log TPM-transformed gene expression of all cells within each metacell was utilized to represent the metacell's expression.

### 2.1 Determining the optimal number of cells per metacell
Gene Coverage

Within metacell variation (GINI index)

### 2.2 Evaluating the performance of metacell
Local inverse Simpson’s Index (LISI)

Biological program signature scores

### 2.3 Generating metacell within each datasets
Each metacell was generated from a specific sample, the clinical information-including tissue origin, treatment, and response conditions-remained consistent with that of the corresponding samples from which they were derived

## 3. Metacell integration 
For the totality of metacells derived from all scRNA-seq datasets, we integrated the metacells, then evaluated and rectified any prevailing batch effects.

To gain more detailed insights into the metacell heterogeneity of specific cell types, we divide all cells into six lineages for downstream analysis, , including cytotoxic lymphocytes (CD8+ T and NK cells), conventional and regulatory lymphocytes (conventional CD4+ T and regulatory T cells), B lymphocytes (B and plasma cells), myeloid cells (monocyte, macrophage, mast, dendritic cells), fibroblasts (fibroblast and myofibroblasts), endothelial cells, and epithelial cells.
### 3.1 Identifying the optimal batch effect correction method
1-Adjusted rand index (1-ARI)

Local Inverse Simpson's Index (LISI)

Entropy

Average Silhouette Width (ASW)
### 3.2 Identifying the optimal clustering resolution

Silhouette score

Clustree

### 3.3 Metacell annotation

Manually curatition

ROGUE

## 4. Spatial transcriptomics (ST) data collection and preprocessing
To investigate the spatial localization of specific cell types and the relative positioning of two particular cell types, we utilized spatial transcriptomic data, encompassing gene-spot matrices, cell coordinates, and corresponding tissue images.

### 4.1 ST Data collection and pre-processing
From the GEO database, we collect the spatial transcriptomics (ST) data of 62 patients from six cancer types, including primary liver cancer (PLC), ovarian carcinoma (OV), pancreatic adenocarcinoma (PAAD), colorectal cancer (CRC), breast cancer (BRCA), and squamous cell carcinoma (SCC).
To safeguard data quality, we filtered the spots based on a minimum detection threshold of 250 genes and a maximum of 25% mitochondrial gene expression, while also excluding genes with fewer than 10 read counts or those expressed in fewer than two spots.

### 4.2 Malignant cell identification 
Malignant cell markers

Copy number variation

Annotation from the original studies

### 4.3 Spot deconvolutaion
STRIDE by matched scRNA-Seq data or selected scRNA-Seq data from TISCH

## 5. Screening of driver cell types
Screening of potential cancer-driven cell types by quantifying their relative abundance across different sources, cancer types, and spatial localization. Furthermore, it facilitated the investigation of cell type-specific functions, and the estimation of their effects on immune cell infiltration and prognosis
### 5.1 Distribution
#### Source preference analysis
Analysis of variance (ANOVA) test

odds ratios

#### Cancer types preference analysis
#### Spatial localization analysis
Colocation analysis

Distance measurement between Fibroblast to malignant cells

### 5.2 Function analysis
Functional associated signature score

Intratumor heterogeneity derived from Non-negative matrix factorization (NMF)

Functional enrichment analysis

Metabolic activity

Cell communication 

### 5.3 Clinical effect
Survival anlaysis

Immune infiltration 

## 6. Application

TabulaTIME defined 56 unique cell types across different cancer types using scRNA-seq. With this high-resolution reference, we then sought to investigate whether we could stratify patients into different tumor subtypes based on their expressed ecotypes.
Cell type annotation is vital for interpreting function phenotypes of cells when analyzing scRNA-seq datasets. A comprehensive and fully annotated dataset is highly needed for reference-based cell type annotation methods and could significantly improve annotation performance. We next tested whether the integrated blueprint from TabulaTIME could served as a reference map for pan-cancer single-cell annotation.
### 6.1 Bulk ecotype analysis
Cell type Ecotypes

Patient Tumor microenviroment subtypes

### 6.2 Automatic cell type annotation
Reference map

Marker gene list

## Environment 
    Ubuntu 9.3.0
    R version 4.0.5	
    Python version 3.8.10	

