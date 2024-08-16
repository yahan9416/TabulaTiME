<a href=""><img src="" alt="DOI"></a>
# TabulaTIME: A Workflow for Integrating Single-Cell RNA Sequencing Data of Tumors and Building Pan-Cancer Landscapes

To characterize the heterogeneous composition and evolution of TME during tumor initiation, progression, and metastasis across different cancer types, we present the Tabula of the Tumor Immune Microenvironment (TabulaTIME) framework. TabulaTIME characterized the pan-cancer TME landscape based on the integration of large-scale scRNA-seq datasets. The framework consists of five major modules: tumor-related scRNA-seq data collection, data pre-processing and MetaCell identification, integration of all lineages, lineage-specific integration, and characterization of cell sub-types. <a href=""><img src="http://timer2.compbio.cn/TabulaTIME" alt="TabulaTIME"></a>  

<img src="https://github.com/yahan9416/TabulaTiME/blob/main/Image/TabulaTIME_workflow.png" alt="Image Description" width="100%" />

## 1. scRNA-seq data collection and preprocessing
To illuminate the intricate dynamics of the tumor microenvironment (TME) throughout the genesis and progression of tumors across a spectrum of cancer types, we aggregated single-cell RNA sequencing (scRNA-Seq) datasets from 746 donors representing 36 distinct cancer varieties. These datasets were meticulously curated from the GEO database and our prior investigations within the Tumor Immune Single-cell Hub (TISCH) database.

All collected datasets underwent rigorous pre-processing through the MAESTRO workflow, which encompassed quality control, elimination of doublets and batch effects, cell clustering, and precise cell type annotation. Ultimately, a total of 103 studies were collated, comprising an impressive 4,479,563 cells.

## 2. MetaCell identification
To mitigate the impact of technical disturbances and minimize computing resource expenses, TabulaTIME employed a strategy of grouping cells with similar expressions into MetaCells within each dataset, ensuring that each MetaCell encompassed approximately 30 cells. The average log TPM-transformed gene expression of all cells within each MetaCell was utilized as a representative measure of the MetaCell's expression.

## 3. Metacell integration 
For the entirety of metacells derived from all scRNA-seq datasets, we conducted an integrative analysis, subsequently assessing and correcting any prevailing batch effects.

To delve deeper into the metacell heterogeneity within specific cell types, we categorized all cells into six distinct lineages for downstream analysis: cytotoxic lymphocytes (CD8+ T and NK cells), conventional and regulatory lymphocytes (conventional CD4+ T and regulatory T cells), B lymphocytes (B and plasma cells), myeloid cells (monocytes, macrophages, mast cells, dendritic cells), fibroblasts (fibroblasts and myofibroblasts), endothelial cells, and epithelial cells.


## 4. Spatial transcriptomics (ST) data collection and preprocessing
To investigate the spatial localization of specific cell types and the relative positioning of two particular cell types, we utilized spatial transcriptomic data, encompassing gene-spot matrices, cell coordinates, and corresponding tissue images.

## 5. Characterization of MetaCell subtypes
To facilitate the characterization of MetaCell subtypes, we examined the distribution of subtypes across various sources and cancer types, as well as the functions and phenotypes of subtypes, and their impact on survival.

## 6. Application

TabulaTIME defined 56 unique cell types across different cancer types using scRNA-seq. With this high-resolution reference, we then sought to investigate whether we could stratify patients into different tumor subtypes based on their expressed ecotypes.

Cell type annotation is vital for interpreting function phenotypes of cells when analyzing scRNA-seq datasets. A comprehensive and fully annotated dataset is highly needed for reference-based cell type annotation methods and could significantly improve annotation performance. We next tested whether the integrated blueprint from TabulaTIME could served as a reference map for pan-cancer single-cell annotation.

## Environment 
    Ubuntu 9.3.0
    R version 4.0.5	
    Python version 3.8.10	

