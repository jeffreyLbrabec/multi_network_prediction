# Convert Gene IDs to Gene Names ------------------------------------------

convert_to_names <- function(data, species, target) {
  
  gene_names <- gconvert(data, 
                         organism = species, 
                         target = target,
                         numeric_ns = "ENTREZGENE_ACC") %>% 
    as_tibble()
  
  return(gene_names)
  
  # converted_data <- data %>%
  #   mutate(gene_names = gene_names$names) %>% 
  #   select(gene_names, class, final_mean_pred)
  # 
  # return(converted_data)
  
}
