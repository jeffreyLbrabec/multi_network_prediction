get_plot_cs <- function(tib, direction, scores, cell_type, comp_group, title) {
  
  score_tib <- tib %>% 
    inner_join(scores, by = c("gene" = "gene_name")) %>% 
    mutate(log_p = -log10(model_p)) %>% 
    rowwise() %>% 
    mutate(comb_score = sum(log_fpr > .$log_fpr & log_p > .$log_p)) %>% 
    ungroup() %>% 
    mutate(comb_score = comb_score/n()) %>% 
    arrange(desc(comb_score))
  
  cdf_plot <- gene_rank_plotr(score_tib, 
                              title = paste(title, direction, sep = " "), 
                              comp_group = comp_group)
  
  ggsave(paste(here("results/figures/cdf_plots"), 
               comp_group, 
               paste0(cell_type, direction, ".png"), 
               sep = "/"))
  
}

get_cs_ont <- function(tib,
                       scores,
                       organism, 
                       ordered_query = TRUE, 
                       significant = TRUE) {
  
  score_tib <- tib %>% 
    inner_join(scores, by = c("gene" = "gene_name")) %>% 
    mutate(log_p = -log10(model_p)) %>% 
    rowwise() %>% 
    mutate(comb_score = sum(log_fpr > .$log_fpr & log_p > .$log_p)) %>% 
    ungroup() %>% 
    mutate(comb_score = comb_score/n()) %>% 
    arrange(desc(comb_score))
  
  module_ontology <- gost(score_tib$gene, 
                          organism = organism, 
                          ordered_query = ordered_query,
                          significant = significant,
                          sources = "GO:BP",
                          correction_method = "fdr")
  
  module_res <- module_ontology$result %>% 
    as_tibble() %>% 
    select(term_id, p_value, term_name)
  
  return(module_res)
  
}