calc_fold_median <- function(fold_obj) {
  
  res_final <- fold_obj %>%
    rowwise() %>% 
    mutate(mean_pred_val = exp(mean(log(c_across(contains("pred"))), na.rm = TRUE))) %>%
    ungroup() %>%
    select(-ends_with("p")) %>%
    mutate(class = as.factor(class)) %>% 
    arrange(class)
  
  auc_val <- roc_auc(res_final, truth = class, mean_pred_val)$.estimate
  
  return(auc_val)

}
