dim_red_and_merge <- function(network_dir = here("data/networks"),
                              merged_net_name = "merged_net",
                              thresh = .99,
                              sig_genes) { #add downsamp here?
  
  adj_mats <- list.files(network_dir, full.names = FALSE)
  adj_mats <- adj_mats[!str_detect(adj_mats, ".RData") & !str_detect(adj_mats, merged_net_name)]
  
  adj_list <- list()
  
  for(i in 1:length(adj_mats)) {
    
    adj_list[[i]] <- read_rds(paste(network_dir, adj_mats[i], sep = "/")) %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "tp_genes") %>% 
      select(tp_genes, everything())
    
  }
  
  intersect_of_colnames <- map(adj_list, colnames) %>% 
    reduce(intersect)
  
  filt_adj_list <- map(adj_list, magrittr::extract, intersect_of_colnames) %>% 
    map(column_to_rownames, var = "tp_genes")
  
  
  # Feature Engineering -----------------------------------------------------
  
  pcaed <- map(filt_adj_list, t) %>% 
    map(scale) %>%
    map(prcomp) %>% 
    map(thresh_pca, thresh = thresh) %>% 
    map(as.data.frame) %>% 
    map(`rownames<-`, colnames(filt_adj_list[[1]])) %>% 
    map(t) %>% 
    map(as.matrix)
  
  names(pcaed) <- c(1:length(pcaed))
  
  
  # Bind Networks Together --------------------------------------------------
  
  merged_net <- do.call("rbind", pcaed)
  
  
  # Add Class Info From sig_genes -------------------------------------------
  
  merged_net_transp <- t(merged_net)
  
  genes <- rownames(merged_net_transp)
  feats <- colnames(merged_net_transp)
  
  merged_tib <- merged_net_transp %>% 
    as_tibble() %>% 
    mutate(genes = factor(genes),
           class = ifelse(genes %in% sig_genes, "P", "U"),
           class = as.factor(class)) %>% 
    select(class, genes, everything()) %>% 
    clean_names()
  
  
  # Save Merged Network Tibble -----------------------------------------------------
  
  write_rds(merged_tib, paste(network_dir, paste(merged_net_name, ".rds", sep = ""), sep = "/")) 
  
  return(merged_tib)
}