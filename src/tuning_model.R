tuning_model <- function(object,
                         boot_id,
                         fold_dir,
                         learner, 
                         cost = 1, 
                         degree = 1) {
  
  dir.create(file.path(fold_dir, boot_id))
  boot_dir <- paste(fold_dir, boot_id, sep = "/")
  
  if(learner == "svm_poly") {
    
    hypers <- tidyr::crossing(cost, degree)
    
    hypers %>%
      mutate(roc = map2(cost, degree, possibly(roc_wrapping, otherwise = NA_real_), object = object, boot_dir = boot_dir, learner = learner))
    
  }else if(learner == "svm_rbf") {
    
    tibble(cost = cost) %>% 
      mutate(roc = map(cost, possibly(roc_wrapping, otherwise = NA_real_), object = object, boot_dir = boot_dir, learner = learner))
    
  }else if(learner == "svm_linear") {
    
    tibble(cost = cost) %>% 
      mutate(roc = map(cost, possibly(roc_wrapping, otherwise = NA_real_), object = object, boot_dir = boot_dir, learner = learner))
    
  }
  
}