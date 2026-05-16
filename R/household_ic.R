#' Aggregate individual-level influence-curve components to household level
#'
#' Sums influence-curve values within each household (cluster) and scales by
#' `num_households / num_records` to keep the aggregated IC on the original
#' scale. This is a reimplementation of `ltmle:::HouseholdIC` so that the
#' package does not rely on a non-exported function.
#'
#' @param record_ic Numeric matrix (or vector) of influence-curve values with
#'   one row per individual record.
#' @param id Integer or character vector indicating cluster membership. Must
#'   have length equal to the number of rows of `record_ic`. If `NULL`,
#'   `record_ic` is returned unchanged.
#'
#' @return A numeric matrix with one row per household (unique `id`), scaled
#'   by `num_households / num_records`. If `record_ic` is a vector, the
#'   result is a one-column matrix.
#'
#' @keywords internal
household_ic <- function(record_ic, id) {
  if (is.null(id)) {
    return(record_ic)
  }
  if (is.null(dim(record_ic))) {
    record_ic <- as.matrix(record_ic)
  }
  if (length(id) != nrow(record_ic)) {
    stop("`id` length must equal the number of rows of `record_ic`.")
  }
  household_ic <- as.matrix(
    stats::aggregate(record_ic, list(id = id), sum)[, -1, drop = FALSE]
  )
  num_records <- nrow(record_ic)
  num_households <- nrow(household_ic)
  household_ic * num_households / num_records
}
