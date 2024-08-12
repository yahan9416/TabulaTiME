##  Run_scRNAseq_data_analysis.sh

### Data collection
Count table 
TPM table

### Quality control

 Cell number per dataset (> 1000)
 
 UMI count per cell (>1000)
 
 Gene number per cell (>500) 
 
 Mitochondrial genes per cell (< 15%)
 
 Doublets removal: Scrublet expected_doublet_rate=0.06
  
### Data Pre-processing

 #### Malignant cell identification
 
 Annotation from the original studies
 
 Copy number variation
 
 Malignant cell markers
  
 #### Batch effect evaluation and correction
 Entropy-based metric
 Canonical Correlation analysis
  
 #### Cell clustering
 Louvain algorithm
  
 #### Cell type annotation:
 MAESTRO based on DE genes
 Manually curatition

 ### 
 To ascertain the detailed process of single-cell RNA sequencing datasets, we kindly direct you to the comprehensive instructions available within the following GitHub repository: https://github.com/DongqingSun96/TISCH/tree/master/code.
