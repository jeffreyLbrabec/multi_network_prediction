pick_gene_score <- function(tib) {
  
  tib %>% 
    select(gene, contains("mixed_model_p")) %>% 
    rowwise() %>% 
    mutate(max_p = max(c_across(contains("mixed_model_p")), na.rm = TRUE)) %>% 
    ungroup() %>% 
    select(gene, model_p = max_p)
  
}

eat_wrapper <- function(first, second, third, fourth) {
  
  set_list <- list(first, second, third, fourth) %>% 
    `names<-`(c(1, 2, 3, 4))
  
  eat(set_list, .mode = "full", .by = "gene", .check = "")
  
}
