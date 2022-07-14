get_ont <- function(data,
                    organism, 
                    ordered_query = FALSE, 
                    significant = TRUE) {
  
  module_ontology <- gost(data$gene, 
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