## Application
### Bulk ecotype analysis
    #Patient Tumor microenviroment subtypes
    #Cell type Ecotypes
    Define_Ecosystem_of_TCGA_dataset.R DEG_gene_path TCGA_exp_path TCGA_clincal_path
    
    #Benchmark with similar studies
    Benchmark_TabulaTIME_similiary_stuides_Patient_TME_subtypes.R TCGA_estimate_Infi_path TabulaTIME_ecotype_path Cancer_cell_path Cell_ten_path Immunity_Six_path
    Benchmark_TabulaTIME_similiary_stuides_ColdHot_tumor.R TCGA_estimate_Infi_path TabulaTIME_ecotype_path Cancer_cell_path Cell_ten_path Immunity_Six_path
    Benchmark_TabulaTIME_similiary_stuides_Cellstates.R Seuratobj_path TabulaTIME_DEG_path TabulaTIME_DEG_SCINA_path Cell_69cellstates_path

### Automatic cell type annotation
    #pre-train Reference map
    Selina_Celltype_Prediction.sh
    TabulaTIME_pretrain_SELINA_params.pt
    TabulaTIME_pretrain_SELINA_meta.pkl
   
