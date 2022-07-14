modeling_roc <- function(object,
                         boot_dir,
                         cost = 1, 
                         degree = 1, 
                         learner = c("svm_poly", "svm_rbf", "svm_linear")) {
  
  
  #set.seed(42)
  
  recipe <- recipe(class ~ ., data = analysis(object)) %>% 
    update_role(genes, new_role = "gene_id") %>%
    step_center(all_predictors()) %>% 
    step_scale(all_predictors()) %>% 
    step_zv(all_predictors()) %>% 
    step_corr(all_predictors(), threshold = 0.4) %>% 
    themis::step_downsample(class) 
  
  juiced_analysis_object <- prep(recipe, retain = TRUE) %>% 
    juice() %>% 
    select(class, genes, everything())
  
  if(learner == "svm_poly") {
    
    mod <-svm_poly(cost = cost, degree = degree) %>% 
      set_engine("kernlab") %>% 
      set_mode("classification") %>% 
      fit(class ~ ., data = juiced_analysis_object)
    
    mod_file <- paste0(paste("cost", cost, "degree", degree, sep = "_"), ".rds")
    
    write_rds(mod, paste(boot_dir, mod_file, sep = "/"))
    
    holdout_pred <- 
      predict(mod, assessment(object) %>% dplyr::select(-class), type = "prob") %>% 
      bind_cols(assessment(object) %>% dplyr::select(class, genes))
    
  }else if(learner == "svm_rbf") {
    
    mod <- svm_rbf(cost = cost) %>% 
      set_engine("kernlab") %>% 
      set_mode("classification") %>% 
      fit(class ~ ., data = juiced_analysis_object)
    
    mod_file <- paste0(paste("cost", cost, sep = "_"), ".rds")
    
    write_rds(mod, paste(boot_dir, mod_file, sep = "/"))
    
    holdout_pred <- 
      predict(mod, assessment(object) %>% dplyr::select(-class), type = "prob") %>% 
      bind_cols(assessment(object) %>% dplyr::select(class, genes))
    
  }else if(learner == "svm_linear") {
    
    mod <- svm_linear(cost = cost) %>% 
      set_engine("kernlab") %>% 
      set_mode("classification") %>% 
      fit(class ~ ., data = juiced_analysis_object)
    
    mod_file <- paste0(paste("cost", cost, sep = "_"), ".rds")
    
    write_rds(mod, paste(boot_dir, mod_file, sep = "/"))
    
    holdout_pred <- 
      predict(mod, assessment(object) %>% dplyr::select(-class), type = "prob") %>% 
      bind_cols(assessment(object) %>% dplyr::select(class, genes))
    
    
  }
}