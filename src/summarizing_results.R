summarizing_results <- function(object, 
                                fold_id,
                                learner,
                                res_dir,
                                cost = 1, 
                                degree = 1) {
  
  dir.create(file.path(res_dir, fold_id))
  fold_dir <- paste(res_dir, fold_id, sep = "/")
  
  if(learner == "svm_poly") {
    
    map2_df(object$splits, object$id, tuning_model, fold_dir = fold_dir, learner = learner, cost = cost, degree = degree) #%>% 
      # group_by(cost, degree) %>%
      # summarize(mean_roc = median(roc, na.rm = TRUE),
      #           n = length(roc),
      #           .groups = "drop")
    
  }else if(learner == "svm_rbf") {
    
    map2_df(object$splits, object$id, tuning_model, fold_dir = fold_dir, learner = learner, cost = cost) #%>% 
      # group_by(cost, rbf_sigma) %>% 
      # summarize(mean_roc = mean(roc, na.rm = TRUE),
      #           n = length(roc),
      #           .groups = "drop")
    
  }else if(learner == "svm_linear") {
    
    map2_df(object$splits, object$id, tuning_model, fold_dir = fold_dir, learner = learner, cost = cost) #%>% 
    # group_by(cost, rbf_sigma) %>% 
    # summarize(mean_roc = mean(roc, na.rm = TRUE),
    #           n = length(roc),
    #           .groups = "drop")
  }
  
}