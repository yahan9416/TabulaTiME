##  Run_scRNAseq_data_analysis.sh

### Step 1. Data collection
All datasets were sourced from the GEO database and our previous work the Tumor Immune Single-cell Hub (TISCH) database.


### Quality control
    # Convert 10X mtx format matrix to HDF5 format.
    MAESTRO mtx-to-h5 --type Gene --matrix matrix.mtx --feature features.tsv --gene-column 2 --barcode barcodes.tsv --species GRCh38 --directory . --outprefix "Sample1"
    # Merge 10X HDF5 files
    MAESTRO merge-h5 --type Gene --h5 sample1.h5 sample2.h5
    # scRNAseq data analysis pipeline
    MAESTRO scrna-analysis --format h5 --matrix ./TabulaTIME_GSE1.h5 \
--count-cutoff 1000 --gene-cutoff 500 --assembly GRCh38 \
--directory . --outprefix TabulaTIME_GSE1

  
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
