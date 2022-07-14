#needs to be object$inner_resamples
bag_boots <- function(object) {
  
  make_bags <- function(object) {
    
    pos <- object %>% 
      filter(class == "P")
    neg <- object %>% 
      filter(class == "U")
    
    boot_df <- bind_rows(pos, neg %>% slice_sample(n = nrow(pos), replace = TRUE))
    
  }
  
  get_indices <- function(training_list, test_list) {
    
    analysis_nrows <- map(training_list, nrow)
    assessment_nrows <- map(test_list, nrow)
    
    ind <- vector("list", 25)
    
    for(i in 1:length(analysis_nrows)) {
      
      ind[[i]] <- list(analysis = seq(analysis_nrows[[i]]),
                       assessment = analysis_nrows[[i]] + seq(assessment_nrows[[i]]))
      
    }
    
    return(ind)
    
  }
  
  
  test_boots <- map(object$splits, assessment)
  train_boots <- map(object$splits, analysis)
  
  bagged_boots <- map(train_boots, make_bags)
  
  combined_boots <- map2(bagged_boots, test_boots, bind_rows)
  
  indices <- get_indices(bagged_boots, test_boots)
  
  final_splits <- map2(indices, combined_boots, make_splits, class = "boot_split") 
  list(splits = final_splits,
       id = names0(length(final_splits), "Bootstrap"))
  
}


