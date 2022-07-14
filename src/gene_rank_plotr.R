gene_rank_plotr <- function(score_tab = NULL, title, comp_group) {
  
  top_10_comb_score <- score_tab %>% 
    slice_max(comb_score, n = 10)
  
  score_tab_plot <-
    ggplot(score_tab,
           aes(x = log_p, y = log_fpr, label = gene)) +
    geom_point(aes(color = comb_score), size = 5) +
    scale_color_viridis() +
    labs(x = "-Log10 P-Value",
         y = "Multi-Network Predicted Score",
         color = "Combined Score",
         title = title,
         subtitle = str_to_title(str_replace_all(comp_group, "_", " ")))
  
  score_tab_plot <- score_tab_plot + geom_text_repel(data = top_10_comb_score, aes(label = gene, fontface = "italic"))
  return(score_tab_plot)
}