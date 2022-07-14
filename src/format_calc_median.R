format_calc_median <- function(results_obj) {
  
  res <- results_obj %>% 
    mutate(roc = map(results_obj$roc, dplyr::select, genes, class, .pred_P, -.pred_U))
  
  noms <- function(roc_group) {
    
    names(roc_group$roc) <- 1:25
    
    tibble(cost = roc_group$cost[1],
           degree = roc_group$degree[1],
           roc = list(eat(roc_group$roc,
                          .by = c("genes", "class"),
                          .mode = "full") %>% 
                        clean_names()))
    
  }
  
  
  res_grouped <- res %>%
    group_split(cost, degree) %>% 
    map(., noms)
  
  res_bound <- do.call(rbind, res_grouped)
  
  res_final <- res_bound %>%
    rowwise() %>% 
    mutate(median_pred_val = median(c(pred_p,
                                      x2_pred_p,
                                      x3_pred_p,
                                      x4_pred_p,
                                      x5_pred_p,
                                      x6_pred_p,
                                      x7_pred_p,
                                      x8_pred_p,
                                      x9_pred_p,
                                      x10_pred_p,
                                      x11_pred_p,
                                      x12_pred_p,
                                      x13_pred_p,
                                      x14_pred_p,
                                      x15_pred_p,
                                      x16_pred_p,
                                      x17_pred_p,
                                      x18_pred_p,
                                      x19_pred_p,
                                      x20_pred_p,
                                      x21_pred_p,
                                      x22_pred_p,
                                      x23_pred_p,
                                      x24_pred_p,
                                      x25_pred_p), 
                                    na.rm = TRUE)) %>%
    ungroup() %>%
    select(-ends_with("p")) %>%
    arrange(class)
  
  roc_val <- roc_auc(res_final, truth = class, median_pred_val)$.estimate
  return(roc_val)
  
}