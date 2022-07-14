pick_boot_model <- function(boot, cost = 1, degree = 1, rbf_sigma = 1, learner, fold_id, results_dir) {
  
  model_files <- list.files(paste(results_dir, fold_id, boot, sep = "/"))
  
  if(learner == "svm_poly") {
    
    best_model <- model_files[str_detect(string = model_files, 
                                         pattern = paste0(paste("cost", cost, "degree", degree, sep = "_"), ".rds"))]
    
  }else if(learner == "svm_rbf") {
    
    best_model <- model_files[str_detect(string = model_files, 
                                         pattern = paste0(paste("cost", cost, sep = "_"), ".rds"))]
    
  }else if(learner == "svm_linear") {
    
    best_model <- model_files[str_detect(string = model_files, 
                                         pattern = paste0(paste("cost", cost, sep = "_"), ".rds"))]
    
  }
  
  best_mod <- read_rds(paste(results_dir, fold_id, boot, best_model, sep = "/"))
  
}