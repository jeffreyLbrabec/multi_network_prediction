get_fold_aucs <- function(final_pred_dir) {
  
  dirs <- tibble(fold_dirs = list.files(final_pred_dir, full.names = TRUE),
                 fold_id = list.files(final_pred_dir)) %>% 
    filter(!str_detect(fold_id, ".rds"))
  
  res <- dirs %>% 
    mutate(results = map(fold_dirs, read_boots)) %>% 
    select(fold_id, results) %>% 
    unnest(results)
  
  return(res)
  
}

read_boots <- function(fold_dir) {
  
  boot_dirs <- list.files(fold_dir, full.names = TRUE)
  
  boot_preds <- map(boot_dirs, read_csv, col_types = list(col_double(), col_factor(), col_double()))
  
  names(boot_preds) <- 1:length(boot_preds)
  
  boot_aucs <- map(boot_preds, auc_wrapper)
  
  auc_tib <- as_tibble(boot_aucs) %>% 
    pivot_longer(cols = everything(), 
                 names_to = "boot_id",
                 values_to = "auc") %>% 
    mutate(boot_id = str_c("Bootstrap", boot_id, sep = " "))
  
}                

auc_wrapper <- function(boot_res) {
  
  roc_auc(boot_res, truth = class, .pred_P)$.estimate
  
}