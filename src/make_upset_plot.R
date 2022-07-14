make_upset_plot <- function(deg_list, title = NULL, subtitle = NULL, fill = "cornflowerblue") {
  
  deg_tib <- eat(deg_list, 
                 .by = "gene",
                 .mode = "full") %>% 
    nest(cell_expr = contains("cell_expr")) %>% 
    mutate(cell_expr = map(cell_expr, 
                           pivot_longer, 
                           cols = everything(), 
                           names_to = 'names', 
                           values_to = 'values'),
           cell_expr = map(cell_expr, filter, !is.na(values)),
           cell_expr = map(cell_expr, select, values),
           cell_expr = map(cell_expr, pull))
  
  deg_tib %>% 
    ggplot(aes(x = cell_expr)) +
    geom_bar(fill = fill) +
    scale_x_upset(order_by = "freq") +
    scale_y_log10() +
    labs(x = "",
         y = "Number of Genes") +
    theme_light()
  
}
