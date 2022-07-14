pick_model <- function(cost = 1, degree = 1, rbf_sigma = 1, learner, fold_id, results_dir) {
  
  boots <- list.files(paste(results_dir, fold_id, sep = "/"))
  
  if(learner == "svm_poly") {
    
    best_models <- map(boots, pick_boot_model, cost = cost, degree = degree, learner = learner, fold_id = fold_id, results_dir = results_dir)
    
  }else if(learner == "svm_rbf") {
    
    best_models <- map(boots, pick_boot_model, cost = cost, learner = learner, fold_id = fold_id, results_dir = results_dir)
    
  }else if(learner == "svm_linear") {
    
    best_models <- map(boots, pick_boot_model, cost = cost, fold_id = fold_id, learner = learner, results_dir = results_dir)
    
  }
  
  fold_tib <- tibble(boots = boots, best_models = best_models)
  
}