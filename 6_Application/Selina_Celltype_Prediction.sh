#Predicting cell type annotations in scRNA-seq datasets
#First step: pre-train the model
selina train --path-in /fs/home/hanya/Project/TME_Immune_difference/Selina_predict/Reference/Pretrain_result/Before_CCA_selected_DEG --path-out /fs/home/hanya/Project/TME_Immune_difference/Selina_predict/Reference/Pretrain_result/Before_CCA_selected_DEG  --outprefix Merger_beforeCCA_selectgene --disease

#Second: directly predicted
selina predict --query-expr ./res/query_expr.txt --model ./res/pre-trained_params.pt --seurat ./res/query_res.rds --path-out ./res/predict/disease --disease