#function to read in REVIGO semantic space csv file and plot

semantic_plot <- function(file, title = "My Semantic Plot", subtitle, disp = 0.15) {
  
  revigo_tib <- read_csv(file, na = c("", "NA", "null")) %>% 
    janitor::clean_names() %>% 
    mutate(eliminated = as.logical(eliminated)) %>%
    filter(eliminated == FALSE) %>% 
    mutate(across(contains("plot_"), as.numeric))
  
  ex <- revigo_tib %>%
    filter(dispensability < disp)
  
  plot_x_range = max(revigo_tib$plot_x) - min(revigo_tib$plot_x)
  plot_y_range = max(revigo_tib$plot_y) - min(revigo_tib$plot_y)
  
  
  revigo_tib %>% 
    ggplot() +
    geom_point(aes(plot_x, plot_y, color = value, size = log_size), alpha = (0.6)) +
    scale_color_viridis_c() +
    geom_point(aes(plot_x, plot_y, size = log_size), shape = 21, fill = "transparent", color = "black") + 
    scale_size(range=c(5, 15)) + 
    geom_label_repel(data = ex, aes(plot_x, plot_y, label = name), color = "black", size = 3) +
    labs(x = "Semantic Space X", 
         y = "Semantic Space Y",
         size = "Log P-Value Size",
         color = "Log P-Value",
         title = title,
         subtitle = subtitle) +
    xlim(min(revigo_tib$plot_x)-plot_x_range/10,max(revigo_tib$plot_x)+plot_x_range/10) +
    ylim(min(revigo_tib$plot_y)-plot_y_range/10,max(revigo_tib$plot_y)+plot_y_range/10) +
    theme_light() +
    theme(legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key.size = unit(0.75, 'lines'))
  
  
}

