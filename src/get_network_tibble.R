get_network_tibble <- function(network_dir = here("data/networks"),
                               tissue_name,
                               sig_genes) {
  
  adj_mat <- read_rds(paste(network_dir, paste(tissue_name, "adj.rds", sep = "_"), sep = "/")) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "genes") %>% 
    select(genes, everything())
  
  # Add Class Info From sig_genes -------------------------------------------
  
  adj_tib <- adj_mat %>% 
    as_tibble() %>% 
    mutate(genes = factor(genes),
           class = ifelse(genes %in% sig_genes, "P", "U"),
           class = as.factor(class)) %>% 
    select(class, genes, everything()) %>% 
    clean_names()
}