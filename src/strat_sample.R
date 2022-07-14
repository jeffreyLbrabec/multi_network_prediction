strat_sample <- function(x, prop, times, ...) {
  n <- nrow(x)
  idx <- purrr::map(rep(n, times), sample, size = floor(n*prop), ...)
  out <- purrr::map_df(idx, function(ind, x) x[sort(ind), "idx"], x = x)
  out$rs_id <- rep(1:times, each = floor(n*prop))
  out
}
