library(pheatmap)
library(igraph)
library(here)
library(tidyverse)
library(gprofiler2)

el = as.matrix(read.csv(here("results/gephi/adj_mats/glia_adj_mat_fpr2.csv")))
#clust = as.numeric(read.csv(here("results/gephi/adj_mats/glia_fpr_2 default node.csv"))[ , 1])

set.seed(42)
g = graph.edgelist(el[ , 1:2], directed = FALSE)
E(g)$weight = as.numeric(el[ , 3])
com = spinglass.community(graph = g, weights = E(g)$weight)
#write_rds(com, here("results/network_analysis/glia_communities.rds"))
com <- read_rds(here("results/network_analysis/glia_communities.rds"))

ucom = sort(unique(com$membership))
adj = as.matrix(get.adjacency(g, attr = "weight"))
perm = NULL
for(i in 1:max(ucom)){
  curr_idx = which(com$membership == ucom[i])
  if(length(curr_idx) > 1){
    curr_deg = rowSums(adj[curr_idx, curr_idx])
    curr_perm = order(curr_deg, decreasing = TRUE)
    perm = c(perm, curr_idx[curr_perm])
    print(rownames(adj)[curr_idx[curr_perm[1]]])
  }
  else{
    perm = c(perm, curr_idx)
  }
}

png(here("results/figures/glia_adj_plot.png"), width = 10, height = 10, units = "in", res = 600)
pheatmap(adj[perm, perm], show_rownames = FALSE, show_colnames = FALSE,
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

glia_member_one <- gost(query = com[[1]], 
                        organism = "hsapiens", 
                        significant = TRUE, 
                        correction_method = "fdr", 
                        sources = "GO:BP")
glia_member_one$result %>% select(term_name)

glia_member_two <- gost(query = com[[2]], 
                        organism = "hsapiens", 
                        significant = TRUE, 
                        correction_method = "fdr", 
                        sources = "GO:BP")
glia_member_two$result %>% select(term_name)

glia_member_three <- gost(query = com[[3]], 
                          organism = "hsapiens", 
                          significant = TRUE, 
                          correction_method = "fdr", 
                          sources = "GO:BP")
glia_member_three$result %>% select(term_name)


