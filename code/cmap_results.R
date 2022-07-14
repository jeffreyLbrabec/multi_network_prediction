#Post-process Cmap results
library(cmapR)
library(here)
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(ggrepel)

comb_score_drugs <- parse_gctx(here("results/cmap/ex_neu_comb_score_drugs.gct"))
ex_neu_deg_drugs <- parse_gctx(here("results/cmap/ex_neu_diff_exp_drugs.gct"))
late_ex_neu_drugs <- parse_gctx(here("results/cmap/late_ad_signal_drugs.gct"))
late_comb_drugs <- parse_gctx(here("results/cmap/late_combined_drug_signal.gct"))

comb_score_drugs_tib <- as_tibble(comb_score_drugs@rdesc) %>% 
  mutate(norm_cs = comb_score_drugs@mat[,1])
ex_neu_deg_drugs_tib <- as_tibble(ex_neu_deg_drugs@rdesc) %>% 
  mutate(norm_cs = ex_neu_deg_drugs@mat[,1])
late_ex_neu_drugs_tib <- as_tibble(late_ex_neu_drugs@rdesc)
late_comb_drugs_tib <- as_tibble(late_comb_drugs@rdesc)

cs_drugs <- comb_score_drugs_tib %>% 
  select(pert_iname, cell_iname, pert_idose, pert_itime, moa, target_name, fdr_q_nlog10, norm_cs) %>% 
  filter(cell_iname %in% c("NEU", "NPC")) %>% 
  #arrange(desc(norm_cs)) %>% 
  mutate(score_type = "comb_score")

deg_drugs <- ex_neu_deg_drugs_tib %>% 
  select(pert_iname, cell_iname, pert_idose, pert_itime, moa, target_name, fdr_q_nlog10, norm_cs) %>% 
  filter(cell_iname %in% c("NEU", "NPC")) %>% 
  #arrange(desc(norm_cs)) %>% 
  mutate(score_type = "deg")

late_drugs <- late_ex_neu_drugs_tib %>% 
  select(pert_iname, cell_iname, pert_idose, pert_itime, moa, target_name, fdr_q_nlog10, raw_cs) %>% 
  filter(moa != "-666", cell_iname %in% c("NEU", "NPC")) %>% 
  arrange(desc(raw_cs)) %>% 
  mutate(score_type = "late")

late_cs_drugs <- late_comb_drugs_tib %>% 
  select(pert_iname, cell_iname, pert_idose, pert_itime, moa, target_name, fdr_q_nlog10, raw_cs) %>% 
  filter(moa != "-666", cell_iname %in% c("NEU", "NPC")) %>% 
  arrange(desc(raw_cs)) %>% 
  mutate(score_type = "late_cs")


plot_tib <- dplyr::full_join(cs_drugs %>% select(pert_iname, cell_iname, pert_idose, pert_itime, moa, cs_norm_cs = norm_cs), 
                             deg_drugs %>% select(pert_iname, cell_iname, pert_idose, pert_itime, moa, deg_norm_cs = norm_cs), 
                             by = c("pert_iname", "cell_iname", "pert_idose", "pert_itime")) %>% 
  mutate(cs_norm_cs = replace_na(cs_norm_cs, 0),
         deg_norm_cs = replace_na(deg_norm_cs, 0))

late_plot_tib <- dplyr::full_join(late_cs_drugs %>% select(pert_iname, cell_iname, pert_idose, pert_itime, moa, cs_raw_cs = raw_cs),
                                  late_drugs %>% select(pert_iname, cell_iname, pert_idose, pert_itime, moa, late_raw_cs = raw_cs),
                                  by = c("pert_iname", "cell_iname", "pert_idose", "pert_itime"))

top_both <- filter(plot_tib, cs_raw_cs < 0 & deg_raw_cs < 0) %>% 
  rowwise() %>% 
  mutate(comb_score = sum(cs_raw_cs < .$cs_raw_cs & deg_raw_cs < .$deg_raw_cs)) %>% 
  ungroup() %>% 
  mutate(comb_score = comb_score/n()) %>% 
  arrange(desc(comb_score)) %>% 
  slice_head(n = 20) %>% 
  select(top_both_drugs = pert_iname, top_both_dose = pert_idose)

top_cs <- plot_tib %>% 
  arrange(cs_raw_cs) %>% 
  slice_head(n = 20) %>% 
  select(top_cs_drugs = pert_iname, top_cs_dose = pert_idose)

top_deg <- plot_tib %>% 
  arrange(deg_raw_cs) %>% 
  slice_head(n = 20) %>% 
  select(top_deg_drugs = pert_iname, top_deg_dose = pert_idose)

bottom_both <- filter(plot_tib, cs_raw_cs > 0 & deg_raw_cs > 0) %>% 
  rowwise() %>% 
  mutate(comb_score = sum(cs_raw_cs > .$cs_raw_cs & deg_raw_cs > .$deg_raw_cs)) %>% 
  ungroup() %>% 
  mutate(comb_score = comb_score/n()) %>% 
  arrange(desc(comb_score)) %>% 
  slice_head(n = 20) %>% 
  select(bottom_both_drugs = pert_iname, bottom_both_dose = pert_idose)

bottom_cs <- plot_tib %>% 
  arrange(desc(cs_raw_cs)) %>% 
  slice_head(n = 20) %>% 
  select(bottom_cs_drugs = pert_iname, bottom_cs_dose = pert_idose)

bottom_deg <- plot_tib %>% 
  arrange(desc(deg_raw_cs)) %>% 
  slice_head(n = 20) %>% 
  select(bottom_deg_drugs = pert_iname, bottom_deg_dose = pert_idose)

bind_cols(top_both, top_cs, top_deg, bottom_both, bottom_cs, bottom_deg) %>% 
  write_csv(here("results/cmap/drug_list.csv"))


top_drug_plot <- ggplot(plot_tib) +
  geom_point(data = filter(plot_tib, deg_raw_cs < 0 & cs_raw_cs < 0),
             aes(x = deg_raw_cs,
                 y = cs_raw_cs)) +
  geom_point(data = filter(plot_tib, pert_iname %in% c("pargyline", 
                                                       "zaldaride",
                                                       "glibenclamide",
                                                       "ebelactone-b",
                                                       "lamotrigine",
                                                       "lypressin",
                                                       "TER-14687",
                                                       "papaverine",
                                                       "brivanib",
                                                       "AS-703026",
                                                       "FR-122047",
                                                       "tetrabenazine",
                                                       "kenpaullone",
                                                       "CHIR-99021",
                                                       "BRD-K33396764",
                                                       "BRD-K68202742",
                                                       "piroxicam",
                                                       "retinol",
                                                       "dipyridamole",
                                                       "olanzapine") & deg_raw_cs < 0 & cs_raw_cs < 0) %>% 
               filter(!cs_raw_cs > -0.3) %>% 
               filter(!(deg_raw_cs > -0.3 & cs_raw_cs > -0.35)),
             aes(x = deg_raw_cs,
                 y = cs_raw_cs),
             color = "darkgreen",
             size = 10) +
  geom_label_repel(data = filter(plot_tib, pert_iname %in% c("pargyline", 
                                                             "zaldaride",
                                                             "glibenclamide",
                                                             "ebelactone-b"#,
                                                             # "lamotrigine",
                                                             # "lypressin",
                                                             # "TER-14687",
                                                             # "papaverine",
                                                             # "brivanib",
                                                             # "AS-703026",
                                                             # "FR-122047",
                                                             # "tetrabenazine",
                                                             # "kenpaullone",
                                                             # "CHIR-99021",
                                                             # "BRD-K33396764",
                                                             # "BRD-K68202742",
                                                             # "piroxicam",
                                                             # "retinol",
                                                             # "dipyridamole",
                                                             #"olanzapine"
                                                             ) & deg_raw_cs < 0 & cs_raw_cs < 0) %>% 
                     filter(!cs_raw_cs > -0.3) %>% 
                     filter(!(deg_raw_cs > -0.3 & cs_raw_cs > -0.35)),
                   aes(x = deg_raw_cs,
                       y = cs_raw_cs,
                       label = paste("Drug: ", pert_iname, '\n MOA:', moa.x))) +
  labs(x = "Differential Expression Score",
       y = "Combined Ranking Score") +
  theme_light()

ggsave(filename = here("results/figures/top_drugs_plot.png"),
       height = 10,
       width = 10,
       units = "in",
       dpi = 600)

full_drug_plot <- ggplot(plot_tib) +
  geom_point(aes(x = deg_norm_cs,
                 y = cs_norm_cs)) +
  geom_label_repel(data = filter(plot_tib, pert_iname == "donepezil"),
                  aes(x = deg_norm_cs,
                      y = cs_norm_cs,
                      label = pert_iname),
                  nudge_y = -0.05) +
  geom_point(data = filter(plot_tib, pert_iname == "donepezil"),
             aes(x = deg_norm_cs,
                 y = cs_norm_cs),
             color = "orange") +
  # geom_hline(yintercept = -1.63, linetype = "dashed", color = "sienna") +
  # geom_vline(xintercept = -1.68, linetype = "dashed", color = "sienna") +
  geom_label_repel(data = filter(plot_tib, deg_norm_cs == 0) %>% slice_tail(n = 3),
                  aes(x = deg_norm_cs,
                      y = cs_norm_cs,
                      label = pert_iname),
                  nudge_x = 0.07) +
  geom_label_repel(data = filter(plot_tib, cs_norm_cs == 0) %>% arrange(deg_norm_cs) %>% slice_head(n = 3),
                  aes(x = deg_norm_cs,
                      y = cs_norm_cs,
                      label = pert_iname),
                  nudge_y = 0.06,
                  nudge_x = 0.02) +
  geom_point(data = filter(plot_tib, deg_norm_cs < -1.55 & cs_norm_cs < -1.55),
             aes(x = deg_norm_cs,
                 y = cs_norm_cs),
             color = "darkgreen") +
  geom_label_repel(data = filter(plot_tib, deg_norm_cs < -1.55 & cs_norm_cs < -1.55),
                  aes(x = deg_norm_cs,
                      y = cs_norm_cs,
                      label = pert_iname)) +
  labs(x = "Differential Expression Score",
       y = "Combined Ranking Score") +
  theme_light()

ggsave(filename = here("results/figures/full_drug_plot.png"),
       height = 8,
       width = 8,
       units = "in",
       dpi = 600)

late_full_drug_plot <- ggplot(late_plot_tib) +
  geom_point(aes(x = late_raw_cs,
                 y = cs_raw_cs)) +
  geom_label_repel(data = filter(late_plot_tib, pert_iname == "donepezil"),
                   aes(x = late_raw_cs,
                       y = cs_raw_cs,
                       label = pert_iname),
                   nudge_y = -0.05) +
  geom_point(data = filter(late_plot_tib, pert_iname == "donepezil"),
             aes(x = late_raw_cs,
                 y = cs_raw_cs),
             color = "orange") +
  geom_hline(yintercept = -0.33, linetype = "dashed", color = "sienna") +
  geom_vline(xintercept = -0.42, linetype = "dashed", color = "sienna") +
  geom_label_repel(data = filter(late_plot_tib, late_raw_cs == 0) %>% slice_tail(n = 3),
                   aes(x = late_raw_cs,
                       y = cs_raw_cs,
                       label = pert_iname),
                   nudge_x = 0.07) +
  geom_label_repel(data = filter(late_plot_tib, cs_raw_cs == 0) %>% arrange(late_raw_cs) %>% slice_head(n = 3),
                   aes(x = late_raw_cs,
                       y = cs_raw_cs,
                       label = pert_iname),
                   nudge_y = 0.06,
                   nudge_x = 0.02) +
  geom_point(data = filter(late_plot_tib, late_raw_cs < -0.42 & cs_raw_cs < -0.33),
             aes(x = late_raw_cs,
                 y = cs_raw_cs),
             color = "darkgreen") +
  geom_label_repel(data = filter(late_plot_tib, late_raw_cs < -0.42 & cs_raw_cs < -0.33),
                   aes(x = late_raw_cs,
                       y = cs_raw_cs,
                       label = pert_iname)) +
  labs(x = "Late Differential Expression Score",
       y = "Combined Ranking Score") +
  theme_light()

plot_tib %>% 
  arrange(desc(cs_raw_cs)) %>% 
  slice_head(n = 20) %>% 
  select(pert_iname, moa.y)

plot_tib %>% 
  arrange(cs_raw_cs) %>% 
  slice_head(n = 20) %>% 
  select(pert_iname, moa.y)

fig <- plot_ly(data = plot_tib, 
               x = ~deg_norm_cs, 
               y = ~cs_norm_cs,
               # Hover text:
               text = ~paste("Drug: ", pert_iname))
saveWidget(fig, here("results/cmap/drug_comparison.html"), selfcontained = F, libdir = "lib")
zip(here("results/cmap/drug_comparison.zip"), c(here("results/cmap/drug_comparison.html"), here("results/cmap/lib")))


late_fig <- plot_ly(data = late_plot_tib,
                    x = ~late_raw_cs,
                    y = ~cs_raw_cs,
                    #Hover text:
                    text = ~paste("Drug: ", pert_iname, 'MOA:', moa.x))
