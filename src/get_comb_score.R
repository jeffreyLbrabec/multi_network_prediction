get_comb_score <- function(tib, scores) {
  
  score_tib <- tib %>% 
    inner_join(scores, by = c("gene" = "gene_name")) %>% 
    mutate(log_p = -log10(model_p)) %>% 
    rowwise() %>% 
    mutate(comb_score = sum(log_fpr > .$log_fpr & log_p > .$log_p)) %>% 
    ungroup() %>% 
    mutate(comb_score = comb_score/n()) %>% 
    arrange(desc(comb_score)) %>% 
    select(gene, log_p, log_fpr, comb_score, cell_type, direction)
  
}