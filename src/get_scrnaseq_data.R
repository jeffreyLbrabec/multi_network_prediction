#Extract the desired values from the excel dataset. Returns a list of 2, UP and DOWN

get_scrnaseq_data <- function(filename, sheet, range) {
  
  genes <- read_xlsx(path = filename, 
                     sheet = sheet,
                     range = range) %>% 
    clean_names() %>% 
    rename(gene = x1, degs_ind_model = de_gs_ind_model, degs_ind_mix_models = de_gs_ind_mix_models)
  
  up_degs <- genes %>% 
    filter(degs_ind_model == TRUE, degs_ind_mix_models == TRUE) %>% 
    select(gene, mixed_model_p, ind_model_fc, mixed_model_z) %>% 
    filter(ind_model_fc > 0) %>% 
    mutate(across(contains("mixed_model"), as.numeric))
  
  down_genes <- genes %>% 
    filter(degs_ind_model == TRUE, degs_ind_mix_models == TRUE) %>% 
    select(gene, mixed_model_p, ind_model_fc, mixed_model_z) %>% 
    filter(ind_model_fc < 0) %>% 
    mutate(across(contains("mixed_model"), as.numeric))
  
  all_genes <- genes %>% 
    filter(degs_ind_model == TRUE, degs_ind_mix_models == TRUE) %>% 
    select(gene, mixed_model_p, ind_model_fc, mixed_model_z) %>% 
    mutate(across(contains("mixed_model"), as.numeric))
  
  deg_list <- list(UP = up_degs, DOWN = down_genes, ALL = all_genes)
  
}