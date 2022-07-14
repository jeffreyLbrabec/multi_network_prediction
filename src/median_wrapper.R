median_wrapper <- function(fold) {
  
  fold %>% 
    mutate(auc = map_dbl(fold$roc, calc_fold_median))
  
}
