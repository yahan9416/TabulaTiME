## characterization of cell sub-types

#### Step 1. Distribution
    #Source preference analysis
    #ANOVA
    ANOVA_celltype_proportion.R meta_info_path Datainfo_path
    
    #odds ratios
    Odds_Ratios_celltype_source_distribution.R meta_info_path
    
    #Cancer types preference analysis
    Cancer_type_preference_analysis.R meta_info_path Datainfo_path Cancer_type_info

### Step 2. Function analysis
    #Functional associated signature score
    Signature_score_T.R SeuratObj_path
    Signature_score_Myeloid.R SeuratObj_path
    
    #Intratumor heterogeneity
    #Investigating the functional programs of each sample within specific cell types using Non-negative Matrix Factorization (NMF)
    Biological_program_generated_NMF.R Seurat_obj_path Column_Majorcelltype
    
    #Identify the meta-programs (MPs) of each cell lineage, and process the MPs with a focus on quality control and annotation.
    NMF_MetaPrograms_identification.R jaccrad_similarity_path Top_Feature gmt_path
  
    
    #Functional enrichment analysis
    library(clusterProfiler)
    enricher()
    
    #Infer the metabolic activity of each metacell based on expression matrix.
    Metabilic_activity_each_cell.R Seurat_obj_path KEGG_Module_path MetaCyc_path Enzyme_relationship
    
    #Cell communication
    CellChat.R Seurat_obj_path
    
    cellphonedb method statistical_analysis ./meta_infp.txt ./exp_count.txt --threads=40 --counts-data gene_name --project-name=celltype_inter >log.txt 
      


### Step 3. Clinical effect
    #Survival anlaysis
    #Evaluate the clinical outcomes associated with each cell type in the TCGA (The Cancer Genome Atlas) dataset.
    Survival_analysis.R TCGA_exp_path TCGA_clinical_path DEG_path
    
    #Compute aggregate z-score and p_value using Stouffer's method.
    Stouffer_Combind_Zscore.r
