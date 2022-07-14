predict_models <- function(best_model, boot_num, splits, fold_dir) {
  
  set.seed(42)
  
  
  holdout_pred <- 
    predict(best_model, assessment(splits) %>% dplyr::select(-class), type = "prob") %>% 
    bind_cols(assessment(splits) %>% dplyr::select(class, genes)) %>% 
    mutate(class = as.factor(class)) %>% 
    select(genes, class, .pred_P)
  
  write_csv(holdout_pred, paste(fold_dir, paste0(boot_num, ".csv"), sep = "/"))
  
}