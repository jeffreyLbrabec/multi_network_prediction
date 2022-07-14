my_inside_resampler <- function(src, df, times = 25) {
  
  pos_df <- neuron_juiced_tib %>% 
    filter(class == "P")
  neg_df <- neuron_juiced_tib %>% 
    filter(class == "U")
  
  num_pos <- nrow(pos_df)
  
  pos_list <- vector("list", times)
  pos_list <- map(pos_list, function(x) x <- pos_df)
  
  neg_list <- vector("list", times)
  neg_list <- map(neg_list, function(x) x <- neg_df)
  
  down_samp_negs <- map(neg_list, slice_sample, n = num_pos, replace = TRUE)
  
  boot_df <- map2(pos_list, down_samp_negs)
  
  df$data <- quote(as.data.frame(src))
  eval(df)
  
  src$inner_resamples
}