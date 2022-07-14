source(here("src/fold_auc_processing.R"))

astro_line_res <- get_fold_aucs(here("results/astro_linear_final_predictions"))

ggplot(astro_line_res, aes(fold_id, auc, color = fold_id, fill = fold_id)) +
  geom_boxplot() +
  geom_jitter() +
  scale_color_viridis_d() +
  scale_fill_viridis_d(alpha = 0.5) +
  labs(x = "",
       y = "AUC",
       title  = "Linear SVM Astrocyte Results") +
  theme_light() +
  theme(legend.position = "false")
  
astro_rbf_dirs <- tibble
