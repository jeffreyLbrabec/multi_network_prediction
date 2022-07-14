library(tidyverse)

amyg_vol_res <- read_csv("~/Documents/Research/Manuscripts/NBFP combined score paper/BRABEC_ET_AL_FRONTIERS/frontiers_submission_materials/Supplementals/suppl_file_8_final_amygdala_scores.csv")
hipp_vol_res <- read_csv("~/Documents/Research/Manuscripts/NBFP combined score paper/BRABEC_ET_AL_FRONTIERS/frontiers_submission_materials/Supplementals/suppl_file_7_final_hippocampus_scores.csv")

ex_neu_degs <- read_rds(here("data/single_cell_data/no_path_v_early_path/ex_neu_np_ep_degs.rds"))
in_neu_degs <- read_rds(here("data/single_cell_data/no_path_v_early_path/in_neu_np_ep_degs.rds"))
ast_degs <- read_rds(here("data/single_cell_data/no_path_v_early_path/ast_np_ep_degs.rds"))
olig_degs <- read_rds(here("data/single_cell_data/no_path_v_early_path/olig_np_ep_degs.rds"))

comb_amyg_auc <- get_deg_aucs(ex_neu_degs$ALL, amyg_vol_res, comb_score_amyg)
comb_hipp_auc <- get_deg_aucs(ex_neu_degs$ALL, hipp_vol_res, comb_score_hipp)

func_amyg_auc <- get_deg_aucs(ex_neu_degs$ALL, amyg_vol_res, mwas_log_amyg_fp)
func_hipp_auc <- get_deg_aucs(ex_neu_degs$ALL, hipp_vol_res, mwas_log_hipp_fp)

mwas_genes <- read_table(here("data/sum_stats_genes.genes.out")) %>% 
  janitor::clean_names() %>% 
  mutate(across(gene, as.character))

mwas_names <- gprofiler2::gconvert(query = mwas_genes$gene, organism = "hsapiens", target = "ENTREZGENE", numeric_ns = "ENTREZGENE_ACC")

mwas_score <- mwas_names %>% 
  filter(!duplicated(input)) %>% 
  select(input, target) %>% 
  inner_join(mwas_genes, by = c("input" = "gene")) %>% 
  select(gene_name = target, p) %>% 
  as_tibble() %>% 
  mutate(log_p = -log10(p))

mwas_auc <- get_deg_aucs(ex_neu_degs$ALL, mwas_score, log_p)

ex_neu_aucs <- tibble(score = c("CS Amygdala", "CS Hippocampus", "FS Amygdala", "FS Hippocampus", "Meta-GWAS"),
                      auc = c(comb_amyg_auc, comb_hipp_auc, func_amyg_auc, func_hipp_auc, mwas_auc)) %>% 
  ggplot(aes(x = fct_reorder(score, auc), y = auc, fill = score)) +
  geom_col() +
  scale_fill_viridis_d() +
  expand_limits(y = 1) +
  labs(x = "",
       y = "AUC") +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   vjust = .5,
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))
ggsave(filename = here("results/figures/ex_neu_aucs_aim2.png"))


# Inhibitory Neuron -------------------------------------------------------

comb_amyg_auc <- get_deg_aucs(in_neu_degs$ALL, amyg_vol_res, comb_score_amyg)
comb_hipp_auc <- get_deg_aucs(in_neu_degs$ALL, hipp_vol_res, comb_score_hipp)

func_amyg_auc <- get_deg_aucs(in_neu_degs$ALL, amyg_vol_res, mwas_log_amyg_fp)
func_hipp_auc <- get_deg_aucs(in_neu_degs$ALL, hipp_vol_res, mwas_log_hipp_fp)

mwas_auc <- get_deg_aucs(in_neu_degs$ALL, mwas_score, log_p)

in_neu_aucs <- tibble(score = c("CS Amygdala", "CS Hippocampus", "FS Amygdala", "FS Hippocampus", "Meta-GWAS"),
                      auc = c(comb_amyg_auc, comb_hipp_auc, func_amyg_auc, func_hipp_auc, mwas_auc)) %>% 
  ggplot(aes(x = fct_reorder(score, auc), y = auc, fill = score)) +
  geom_col() +
  scale_fill_viridis_d() +
  expand_limits(y = 1) +
  labs(x = "",
       y = "AUC",
       title = "Inhibitory Neuron DEG AUCs") +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   vjust = .5,
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))
ggsave(filename = here("results/figures/in_neu_aucs_aim2.png"), plot = in_neu_aucs)


# Astrocyte ---------------------------------------------------------------

comb_amyg_auc <- get_deg_aucs(ast_degs$ALL, amyg_vol_res, comb_score_amyg)
comb_hipp_auc <- get_deg_aucs(ast_degs$ALL, hipp_vol_res, comb_score_hipp)

func_amyg_auc <- get_deg_aucs(ast_degs$ALL, amyg_vol_res, mwas_log_amyg_fp)
func_hipp_auc <- get_deg_aucs(ast_degs$ALL, hipp_vol_res, mwas_log_hipp_fp)

mwas_auc <- get_deg_aucs(ast_degs$ALL, mwas_score, log_p)

ast_aucs <- tibble(score = c("CS Amygdala", "CS Hippocampus", "FS Amygdala", "FS Hippocampus", "Meta-GWAS"),
                   auc = c(comb_amyg_auc, comb_hipp_auc, func_amyg_auc, func_hipp_auc, mwas_auc)) %>% 
  ggplot(aes(x = fct_reorder(score, auc), y = auc, fill = score)) +
  geom_col() +
  scale_fill_viridis_d() +
  expand_limits(y = 1) +
  labs(x = "",
       y = "AUC",
       title = "Astrocyte DEG AUCs") +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   vjust = .5,
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))
ggsave(filename = here("results/figures/ast_aucs_aim2.png"), plot = ast_aucs)


# Oligodendrocyte ---------------------------------------------------------

comb_amyg_auc <- get_deg_aucs(olig_degs$ALL, amyg_vol_res, comb_score_amyg)
comb_hipp_auc <- get_deg_aucs(olig_degs$ALL, hipp_vol_res, comb_score_hipp)

func_amyg_auc <- get_deg_aucs(olig_degs$ALL, amyg_vol_res, mwas_log_amyg_fp)
func_hipp_auc <- get_deg_aucs(olig_degs$ALL, hipp_vol_res, mwas_log_hipp_fp)

mwas_auc <- get_deg_aucs(olig_degs$ALL, mwas_score, log_p)

olig_aucs <- tibble(score = c("CS Amygdala", "CS Hippocampus", "FS Amygdala", "FS Hippocampus", "Meta-GWAS"),
                   auc = c(comb_amyg_auc, comb_hipp_auc, func_amyg_auc, func_hipp_auc, mwas_auc)) %>% 
  ggplot(aes(x = fct_reorder(score, auc), y = auc, fill = score)) +
  geom_col() +
  scale_fill_viridis_d() +
  expand_limits(y = 1) +
  labs(x = "",
       y = "AUC",
       title = "Oligodendrocyte DEG AUCs") +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   vjust = .5,
                                   size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))
ggsave(filename = here("results/figures/olig_aucs_aim2.png"), plot = olig_aucs)

