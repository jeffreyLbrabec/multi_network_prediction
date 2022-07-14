get_final_scores <- function(res_dir) {
  
  res_files <- list.files(here(res_dir), 
                          full.names = TRUE)
  
  fold_dirs <- res_files[!str_detect(res_files, ".rds")]
  
  fold_ids <- str_extract(fold_dirs, "Fold..")
  
  res_tib <- tibble(fold_id = fold_ids)
  
  all_preds <- res_tib %>% 
    mutate(preds = map(fold_dirs,
                       pred_wrapper))
  
  names(all_preds$preds) <- 1:nrow(all_preds)
  
  joined_preds <- eat(all_preds$preds,
                      .by = c("genes", "class"),
                      .mode = "full") %>% 
    clean_names()
  
  all_fold_mean_preds <- joined_preds %>% 
    rowwise() %>% 
    mutate(final_mean_pred = exp(mean(log(c_across(contains("pred"))), na.rm = TRUE))) %>% 
    ungroup() %>% 
    select(-ends_with("mean_pred_val")) %>% 
    mutate(class = as.factor(class)) %>% 
    arrange(class)
  
  return(all_fold_mean_preds)
  
}

pred_wrapper <- function(fold_dir) {
  
  boot_files <- list.files(fold_dir, full.names = TRUE)
  
  boot_list <- map(boot_files,
                   read_csv,
                   show_col_types = FALSE)
  
  names(boot_list) <- 1:length(boot_list)
  
  preds <- eat(boot_list,
               .by = c("genes", "class"),
               .mode = "full") %>% 
    clean_names()
  
  comb_boot_med_preds <- preds %>% 
    rowwise() %>% 
    mutate(mean_pred_val = exp(mean(log(c_across(contains("pred"))), na.rm = TRUE))) %>%
    ungroup() %>%
    select(-ends_with("pred_p")) %>%
    mutate(class = as.factor(class)) %>% 
    arrange(class)
  
}



