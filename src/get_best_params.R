get_best_params <- function(object, learner = c("svm_poly", "svm_rbf")) {
  
  if(learner == "svm_poly") {
    
    pooled_inner_object <- object %>% bind_rows
    
    pooled_inner_object[which.max(pooled_inner_object$mean_roc),]
    
  }
  
}