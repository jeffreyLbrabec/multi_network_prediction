get_fpr <- function(model_res) {
  
  arr_res <- model_res %>% 
    dplyr::mutate(class = as.factor(class)) %>% 
    dplyr::arrange(final_mean_pred) 
  
  mod_curve <- yardstick::roc_curve(data = arr_res, truth = class, final_mean_pred) %>% 
    dplyr::filter(!is.infinite(.threshold))
  
  sens_spec_res <- bind_cols(arr_res, mod_curve) %>% 
    dplyr::mutate(fpr = 1 - specificity,
                  log_fpr = -log10(fpr)) %>% 
    dplyr::select(gene_name, class, final_mean_pred, sensitivity, fpr, log_fpr)
  
  return(sens_spec_res)
  
}