strata_check <- function(strata, data) {
  if (!is.null(strata)) {
    if (!is.character(strata) | length(strata) != 1) {
      rlang::abort("`strata` should be a single character value.")
    }
    if (inherits(data[, strata], "Surv")) {
      rlang::abort("`strata` cannot be a `Surv` object. Use the time or event variable directly.")
    }
    if (!(strata %in% names(data))) {
      rlang::abort(strata, " is not in `data`.")
    }
  }
  invisible(NULL)
}
