library(tidyverse)
library(tidymodels)
library(here)
library(janitor)
library(gprofiler2)
library(magrittr)
library(cluster)
library(fpc)
library(powerjoin)

#Source in Functions
fun_dir <- list.files(here("src"), full.names = TRUE)

for(i in 1:length(fun_dir)) {
  
  source(fun_dir[i])
  
}

#Get Significant Genes
mwas_genes <- read_table(here("data/sum_stats_genes.genes.out"),
                         col_names = TRUE) %>% 
  mutate(across(starts_with("GENE"), as.character)) %>% 
  clean_names()

mwas_gconvert <- gconvert(mwas_genes$gene,
                          organism = "hsapiens",
                          target = "ENTREZGENE",
                          numeric_ns = "ENTREZGENE_ACC") %>% 
  filter(!duplicated(name), !duplicated(input)) %>% 
  as_tibble() %>% 
  dplyr::select(gene = input, name) %>%  
  mutate(across(gene, as.character))

mwas_genes %<>% inner_join(mwas_gconvert, by = "gene")

sig_genes <- mwas_genes %>% 
  select(gene, p) %>% 
  filter(p < 0.01)

#Make Network Tibble
neuron_net <- read_rds(here("data/human_networks/neuron_top.RData"))

neuron_adj_mat <- tissue.adj.mat(tissue.net = neuron_net,
                                 gene.list = sig_genes$gene, 
                                 verbose = TRUE)

write_rds(neuron_adj_mat, here("data/human_networks/neuron_adj.rds"))
neuron_adj_mat <- read_rds(here("data/human_networks/neuron_adj.rds"))

neuron_tibble <- get_network_tibble(network_dir = here("data/human_networks"),
                                    tissue_name = "neuron",
                                    sig_genes = sig_genes$gene)
write_csv(neuron_tibble, here("data/neuron_tibble.csv"))

neuron_mat <- neuron_tibble %>%
  filter(class == "P") %>% 
  select(-class) %>% 
  column_to_rownames(var = "genes") %>% 
  as.matrix()

neuron_kmeds <- pam(x = neuron_mat, 
                    k = 500)
neuron_kmeds_keep <- str_c("x", rownames(neuron_kmeds$medoids))

neuron_kmeds_tibble <- neuron_tibble %>% 
  select(class, genes, all_of(neuron_kmeds_keep))

write_csv(neuron_kmeds_tibble, here("data/neuron_kmeds_tibble.csv"))

#--------------
glia_net <- read_rds(here("data/human_networks/glia_top.RData"))
glia_adj_mat <- tissue.adj.mat(tissue.net = glia_net,
                               gene.list = sig_genes$gene,
                               verbose = TRUE)
write_rds(glia_adj_mat, here("data/human_networks/glia_adj.rds"))
glia_adj_mat <- read_rds(here("data/human_networks/glia_adj.rds"))

glia_tibble <- get_network_tibble(network_dir = here("data/human_networks"),
                                  tissue_name = "glia",
                                  sig_genes = sig_genes$gene)
write_csv(glia_tibble, here("data/glia_tibble.csv"))

glia_mat <- glia_tibble %>% 
  filter(class == "P") %>% 
  select(-class) %>% 
  column_to_rownames(var = "genes") %>% 
  as.matrix()

glia_kmeds <- pam(x = glia_mat,
                  k = 500)

glia_kmeds_keep <- str_c("x", rownames(glia_kmeds$medoids))
glia_kmeds_tibble <- glia_tibble %>% 
  select(class, genes, all_of(glia_kmeds_keep))

write_csv(glia_kmeds_tibble, here("data/glia_kmeds_tibble.csv"))

#--------------
astro_net <- read_rds(here("data/human_networks/astrocyte_top.RData"))
astro_adj_mat <- tissue.adj.mat(tissue.net = astro_net,
                                gene.list = sig_genes$gene,
                                verbose = TRUE)
write_rds(astro_adj_mat, here("data/human_networks/astrocyte_adj.rds"))
astro_adj_mat <- read_rds(here("data/human_networks/astrocyte_adj.rds"))

astro_tibble <- get_network_tibble(network_dir = here("data/human_networks"),
                                   tissue_name = "astrocyte",
                                   sig_genes = sig_genes$gene)
write_csv(astro_tibble, here("data/astro_tibble.csv"))

astro_mat <- astro_tibble %>% 
  filter(class == "P") %>% 
  select(-class) %>% 
  column_to_rownames(var = "genes") %>% 
  as.matrix()

astro_kmeds <- pam(x = astro_mat,
                   k = 500)
astro_kmeds_keep <- str_c("x", rownames(astro_kmeds$medoids))
astro_kmeds_tibble <- astro_tibble %>% 
  select(class, genes, all_of(astro_kmeds_keep))

write_csv(astro_kmeds_tibble, here("data/astro_kmeds_tibble.csv"))



multi_net_mat <- astro_kmeds_tibble %>% 
  select(-class) %>% 
  inner_join(glia_kmeds_tibble, by = "genes") %>% 
  select(-class) %>% 
  inner_join(neuron_kmeds_tibble, by = "genes") %>% 
  select(class, genes, everything()) 

write_rds(multi_net_mat, here("kmeds_analysis/multi_net_mat.RData"))

multi_net_tib <- astro_kmeds_tibble %>% 
  select(-class) %>% 
  inner_join(glia_kmeds_tibble, by = "genes") %>% 
  select(-class) %>% 
  inner_join(neuron_kmeds_tibble, by = "genes") %>% 
  select(class, genes, everything())
any(is.na(multi_net_tib))

write_csv(multi_net_tib, here("data/multi_net_kmeds_tibble.csv"))
huh <- read_csv(here("data/multi_net_kmeds_tibble.csv"))


# Trying to take the medians instead of merging the networks --------------

merge_ids_feats <- function(feats, genes, ids) {
  
  genes <- enquo(genes)
  
  merged <- bind_cols(ids, !!genes := feats)
  
  return(merged)
  
}

split_feats <- function(net_tib) {
  
  net_list <- as.list(net_tib)
  
  ids <- bind_cols(net_list[1:2])
  feats <- net_list[-c(1:2)]
  genes <- names(feats)
  
  split_net <- map2(feats, genes, merge_ids_feats, ids)
  
}

split_neuron <- split_feats(neuron_tibble)
split_astro <- split_feats(astro_tibble)
split_glia <- split_feats(glia_tibble)

huh <- eat(split_neuron[[1]], split_astro[[1]], split_glia[[1]], .by = c("genes", "class"), .mode = "full")


boop <- map2(split_neuron, split_astro, full_join, by = c("genes", "class")) %>% 
  map2(split_glia, full_join, by = c("genes", "class"))


boop

get_med_edge <- function(split_edges) {
  
  split_edges %>% 
    rowwise() %>% 
    mutate(med_edge = median(c_across(starts_with("x")), na.rm = TRUE)) %>% 
    ungroup() %>% 
    select(class, genes, med_edge)
  
}

huh <- map(boop, get_med_edge)

edge_names <- names(huh)
huh_named <- map2(huh, edge_names, name_cols)

name_cols <- function(tib, edge_names) {
  
  renamed_tib <- tib %>% rename_with(~ str_replace(.x, "med_edge", edge_names), .cols = med_edge)
}

multi_net_med <- huh_named %>% 
  purrr::reduce(full_join, by = c("genes", "class"))


write_csv(multi_net_med, here("data/new_multi_net_med.csv"))
