get_deg_aucs <- function(degs, 
                         scores,
                         score_col,
                         gwas_only = FALSE) {
  
  if(gwas_only) {
    
    gwas_scores <- scores %>% 
      mutate(true_class = as.factor(ifelse(name %in% degs$gene, "P", "U"))) 
    
    gwas_auc <- roc_auc(gwas_scores, truth = true_class, p)$.estimate
    
    gwas_tib <- tibble(score_type = "gwas_p",
                       auc = gwas_auc)
    
    return(gwas_tib)
    
  }else{
    
    score_ord <- scores %>% 
      mutate(true_class = as.factor(ifelse(gene_name %in% degs$gene, "P", "U"))) %>% 
      select(gene_name, true_class, {{score_col}})
    
    auc <- roc_auc(score_ord, truth = true_class, {{score_col}})$.estimate
    
    return(auc)
    
    # score_cols <- colnames(score_ord)[-1:-2]
    # 
    # roc_list <- vector("list", length = length(score_cols))
    # 
    # names(roc_list) <- score_cols
    # 
    # for(i in 1:length(roc_list)) {
    #   
    #   roc_list[[i]] <- roc_auc(score_ord, 
    #                            truth = true_class, 
    #                            score_cols[[i]])$.estimate
    #   
    # }
    # 
    # roc_tib <- as_tibble(roc_list) %>% 
    #   pivot_longer(everything(),
    #                names_to = "score_type",
    #                values_to = "auc")
  } 
  
}