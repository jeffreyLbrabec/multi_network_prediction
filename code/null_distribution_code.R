#Null distribution for overlap
num_perm = 1e5
null_overlap = matrix(NA, nrow = num_perm, ncol = 1)
for(n in 1:num_perm){
  null_score_fun = function(x){
    rand_samp = sample(1:length(gene_score), length(x), replace = FALSE)
    as.matrix(gene_score[rand_samp, , drop = FALSE])
  }
  null_dist = sapply(new_snp_window_genes, null_score_fun)
  names(null_dist) = names(new_snp_window_genes)
  max_fun = function(x){
    c(rownames(x)[which.max(x)], x[which.max(x)])
  }
  null_max = t(sapply(null_dist, max_fun))
  choice_gene = matrix(NA, nrow = length(snp_info$gene_name), ncol = 2)
  for(i in 1:dim(choice_gene)[1]){
    curr_name = snp_info$gene_name[i]
    if(curr_name %in% rownames(all_dist[[i]])){
      choice_gene[i, 1] = curr_name
      choice_gene[i, 2] = all_dist[[i]][curr_name, ]
    }
  }
  keep_idx = !is.na(choice_gene[ , 1]) & !is.na(null_max[ , 1])
  null_overlap[n] = sum(choice_gene[keep_idx, 1] == null_max[keep_idx, 1])
}
true_overlap = sum(choice_gene[keep_idx, 1] == max_gene[keep_idx, 1])

true_overlap = sum(choice_gene[keep_idx, 1] %in% gene_tib_mn$our_max_gene)

# Get scores for genes in each window
new_snp_window_genes = readRDS(here("data/new_snp_window_genes.rds"))
linear_mn_final_scores <- read_csv(here("results/final_scores/linear_mn_final_scores.csv"))
linear_mn_fpr <- get_fpr(linear_mn_final_scores)
gene_score = as.matrix(linear_mn_fpr[ , "log_fpr"])
rownames(gene_score) = as.matrix(linear_mn_fpr[ , "gene_name"])
max_score = function(x){
  max(gene_score[x, ])
}
true_dist = sapply(new_snp_window_genes, max_score)
all_score_fun = function(x){
  as.matrix(gene_score[x, ])
}
all_dist = sapply(new_snp_window_genes, all_score_fun)
names(all_dist) = names(new_snp_window_genes)
# Step through each window and look at sorted scores
snp_info <- list(snp_id = snp_locations$refsnp_id,
                 snp_loc = snp_locations$chrom_start,
                 gene_name = snp_locations$gene,
                 chr_name = snp_locations$chr_name)
k = 1
snp_info$gene_name[k]
curr_scores = all_dist[[snp_info$snp_id[k]]]
curr_scores[order(curr_scores, decreasing = TRUE), , drop = FALSE]

gene_tib_mn <- tibble(rsid = names(new_snp_windows),
                      interval = paste0(snp_info$chr_name, ":", if_else(snp_info$snp_loc - 1000000 < 0, 0, snp_info$snp_loc - 1000000), "-", snp_info$snp_loc + 1000000),
                      bellenguez_gene = snp_info$gene_name,
                      our_max_gene = map_chr(new_snp_windows, get_our_max_gene),
                      our_max_score = map_dbl(new_snp_windows, get_our_max_score),
                      our_genes_over_one = map_chr(new_snp_windows, get_our_sig_genes)) %>% 
  arrange(rsid) %>% 
  mutate(tier = snp_locations$tier)
