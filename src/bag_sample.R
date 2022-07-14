bag_sample <- function(x, num_pos, prop, times, ...) {
  
  num_neg <- nrow(x) - num_pos
  
  idx_pos <- purrr::map(rep(num_pos, times), 
                        sample, 
                        size = floor(num_pos * prop), ...)
  
  idx_neg <- purrr::map(rep(num_neg, times), 
                        sample, 
                        size = floor(num_pos * prop), ...)
  
  idx <- purrr::map2(idx_pos, idx_neg, c)
  
  out <- purrr::map_df(idx, function(ind, x) x[sort(ind), "idx"], x = x)
  
  out$rs_id <- rep(1:times, each = 2 * floor(num_pos * prop))
  
  out
}
