format_folds <- function(results_obj, learner) {
  
  filt_roc <- results_obj %>% 
    filter(!is.na(roc))
  
  res <- filt_roc %>% 
    mutate(roc = map(filt_roc$roc, dplyr::select, genes, class, .pred_P))
  
  noms <- function(roc_group) {
    
    names(roc_group$roc) <- 1:nrow(roc_group)
    
    tibble(cost = roc_group$cost[1],
           #degree = roc_group$degree[1],
           roc = list(eat(roc_group$roc,
                          .by = c("genes", "class"),
                          .mode = "full") %>% 
                        clean_names()))
    
  }
  
  if(learner == "svm_poly") {
    
    res_grouped <- res %>% 
      group_split(cost, degree) %>% 
      map(., noms)
    
  }else if(learner == "svm_rbf") {
    
    res_grouped <- res %>% 
      group_split(cost) %>% 
      map(., noms)
    
  }else if(learner == "svm_linear") {
    
    res_grouped <- res %>%
      group_split(cost) %>%
      map(., noms)
    
  }
  
  res_bound <- do.call(rbind, res_grouped)
  
  return(res_bound)
  
}
