split_unnamed <- function(x, f) {
  out <- split(x, f)
  unname(out)
}
