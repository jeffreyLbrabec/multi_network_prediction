library(igraph)
library(here)
library(tidyverse)

astro_kmeds_tibble <- read_csv(here("data/astro_kmeds_tibble.csv"))
glia_kmeds_tibble <- read_csv(here("data/glia_kmeds_tibble.csv"))
neuron_kmeds_tibble <- read_csv(here("data/neuron_kmeds_tibble.csv"))
multi_net_kmeds_tibble <- read_csv(here("data/multi_net_kmeds_tibble.csv"))

astro_kmeds_mat <- astro_kmeds_tibble %>% 
  select(-class) %>% 
  slice_head(n = 500) %>% 
  column_to_rownames("genes") %>% 
  as.matrix()

astro_graph <- graph_from_adjacency_matrix(astro_kmeds_mat, 
                                           mode = "undirected", 
                                           weighted = TRUE)
cohesion(astro_graph) #385

glia_kmeds_mat <- glia_kmeds_tibble %>% 
  select(-class) %>% 
  slice_head(n = 500) %>% 
  column_to_rownames("genes") %>% 
  as.matrix()

glia_graph <- graph_from_adjacency_matrix(glia_kmeds_mat, 
                                           mode = "undirected", 
                                           weighted = TRUE)
cohesion(glia_graph) #305

neuron_kmeds_mat <- neuron_kmeds_tibble %>% 
  select(-class) %>% 
  slice_head(n = 500) %>% 
  column_to_rownames("genes") %>% 
  as.matrix()

neuron_graph <- graph_from_adjacency_matrix(neuron_kmeds_mat, 
                                          mode = "undirected", 
                                          weighted = TRUE)
cohesion(neuron_graph) #93

multi_net_kmeds_mat <- multi_net_kmeds_tibble %>% 
  select(-class) %>% 
  slice_head(n = 1500) %>% 
  column_to_rownames("genes") %>% 
  as.matrix()

multi_net_graph <- graph_from_adjacency_matrix(multi_net_kmeds_mat, 
                                            mode = "undirected", 
                                            weighted = TRUE)
cohesion(multi_net_graph) #93
