## characterization of cell sub-types

#### Step 1. Distribution
    #Source preference analysis
    #ANOVA
    ANOVA_celltype_proportion.R
    #odds ratios
    
    Odds_Ratios_celltype_source_distribution.R
    #Cancer types preference analysis
    Cancer_type_preference_analysis.R

### Step 2. Function analysis
    #Functional associated signature score
    Signature_score_T.R
    Signature_score_Myeloid.R
    
    #Intratumor heterogeneity
    #Non-negative matrix factorization (NMF)
    NMF_MetaPrograms_identification.R
    Biological_program_generated_NMF.R
    
    #Functional enrichment analysis
    library(clusterProfiler)
    enricher()
    
    #Metabolic activity
    Metabilic_activity_each_cell.R
    
    #Cell communication
    CellChat.R
    cellPhoneDB.sh


### Step 3. Clinical effect
    #Survival anlaysis
    Survival_analysis.R
    Stouffer_Combind_Zscore.r

