#' Linearly rescale a numeric outcome to `[0, 1]`
#'
#' Returns the rescaled vector together with the `scale_max` and `scale_min`
#' used so they can be passed downstream to the variance-unscaling step.
#' If `min(y) >= 0` then `scale_min = 0`; if `max(y) <= 1` then `scale_max = 1`.
#' These shortcuts mean that binary or already-bounded outcomes are passed
#' through unchanged.
#'
#' @param y Numeric vector.
#'
#' @return A list with components `y` (rescaled vector), `scale_min`, `scale_max`.
#'
#' @keywords internal
scale_outcome <- function(y) {
  scale_max <- if (max(y) > 1) max(y) else 1
  scale_min <- if (min(y) < 0) min(y) else 0
  y_scaled <- (y - scale_min) / (scale_max - scale_min)
  list(y = y_scaled, scale_min = scale_min, scale_max = scale_max)
}

#' Coerce a data argument to a data.table copy without mutating the caller's input
#'
#' @param data A data.frame, data.table, or tibble.
#' @keywords internal
as_dt_copy <- function(data) {
  if (data.table::is.data.table(data)) {
    data.table::copy(data)
  } else {
    data.table::as.data.table(data)
  }
}

#' Internal: stop with a formatted message
#' @keywords internal
stopf <- function(fmt, ...) stop(sprintf(fmt, ...), call. = FALSE)

#' Internal: emit a message gated by `verbose`
#' @keywords internal
vmsg <- function(verbose, fmt, ...) {
  if (isTRUE(verbose)) message(sprintf(fmt, ...))
}
