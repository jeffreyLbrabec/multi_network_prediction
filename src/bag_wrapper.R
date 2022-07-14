bag_wrapper <- function(splits, models, fold_id, res_dir, tissue) {
  
  if(!dir.exists(file.path(res_dir, fold_id))) {
    
    dir.create(file.path(res_dir, fold_id))
    fold_dir <- paste(res_dir, fold_id, sep = "/")
    
  }
  
  boot_preds <- map2(models$best_models, 
                     models$boots, 
                     predict_models, 
                     splits = splits,
                     fold_dir = fold_dir) 
  
  names(boot_preds) <- 1:length(boot_preds)
  
  comb_boot_preds <- eat(boot_preds,
                         .by = c("genes", "class"),
                         .mode = "full") %>% 
    clean_names()
  
  comb_boot_med_preds <- comb_boot_preds %>% 
    rowwise() %>% 
    mutate(mean_pred_val = exp(mean(log(c_across(contains("pred"))), na.rm = TRUE))) %>%
    ungroup() %>%
    select(-starts_with("pred_p")) %>%
    mutate(class = as.factor(class)) %>% 
    arrange(class)
  
  write_rds(comb_boot_med_preds, paste(res_dir, paste(fold_id, tissue, "comb_boot_med_preds.rds", sep = "_"), sep = "/"))
  
  auc_val <- roc_auc(comb_boot_med_preds, truth = class, mean_pred_val)$.estimate
  
  
}