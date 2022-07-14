library(tidyverse)
library(here)
library(Seurat)
library(SingleCellExperiment)
library(readxl)

load("~/Downloads/subclustering_out/subclust_markers_DE_list.RData")
load("~/Downloads/subclustering_out/subclust.out.by.celltype.RData")

filtered_mat <- Matrix::readMM(here("data/single_cell_data/tsne_data/filtered_count_matrix.mtx"))
gene_rows <- read_table(here("data/single_cell_data/tsne_data/filtered_gene_row_names.txt"),
                        col_names = FALSE)
col_data <- read_table(here("data/single_cell_data/tsne_data/filtered_column_metadata.txt"))

id_map <- read_csv(here("data/single_cell_data/tsne_data/snRNAseqPFC_BA10_id_mapping.csv"))
samp_path <- read_xlsx(here("data/single_cell_data/tsne_data/41586_2019_1195_MOESM5_ESM.xlsx")) %>% 
  arrange(Subject)

path_data <- id_map %>%
  summarize(projid = unique(projid), Subject = unique(Subject)) %>% 
  arrange(Subject) %>% 
  full_join(samp_path, by = "Subject") %>% 
  select(projid, subject = Subject, path_group = pathology.group) %>% 
  mutate(projid = as.factor(projid))

tsne_data <- col_data %>% 
  mutate(projid = as.factor(projid)) %>% 
  left_join(path_data, by = "projid")

group_labs <- tsne_data %>% 
  group_by(broad.cell.type) %>% 
  summarize(mean_tsne1 = mean(tsne1),
            mean_tsne2 = mean(tsne2))

tsne_plot <- tsne_data %>% 
  ggplot() +
  geom_point(aes(x = tsne1, y = tsne2, color = broad.cell.type)) +
  geom_label_repel(data = group_labs, 
                   aes(x = mean_tsne1, 
                       y = mean_tsne2, 
                       label = broad.cell.type),
                   size = 10) +
  scale_color_viridis_d() +
  labs(x = "TSNE-1",
       y = "TSNE-2",
       color = "Cell-type") +
  theme_light() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 25))

ggsave(here("results/figures/tsne_plot.png"),
       tsne_plot,
       height = 10,
       width = 10,
       units = "in",
       dpi = 600)




ex_mn_genes

sig_ex_mat <- filtered_mat[which(gene_rows$X1 %in% ex_mn_genes$gene), ]

jaccard <- function(a, b) {
  intersection <-  length(intersect(a, b))
  union <-  length(a) + length(b) - intersection
  return(intersection/union)
}

func_sig_mn <- linear_mn_fpr %>% 
  filter(log_fpr > 1)

jaccard(ex_neu_np_ep_degs$ALL$gene, func_sig_mn$gene_name)
jaccard(in_neu_np_ep_degs$ALL$gene, func_sig_mn$gene_name)