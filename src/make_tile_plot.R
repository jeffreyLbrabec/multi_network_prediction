#Function to wrap up the tile plots

make_tile_plot <- function(auc_tibble, title, subtitle) {
  
  grid_plot <- auc_tibble %>% 
    ggplot(aes(tissue, cell)) +
    geom_tile(aes(fill = auc)) +
    geom_label(aes(label = round(auc, digits = 2))) +
    scale_fill_viridis_c() +
    labs(x = "Tissue Network",
         y = "RNAseq Cells",
         fill = "AUC",
         title = title,
         subtitle = subtitle) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid = element_blank())
  
  return(grid_plot)
  
}