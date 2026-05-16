#' Add an `alpha` column of analysis weights
#'
#' Computes the weight `alpha` for each row of the input data per the
#' weighting scheme of Benitez et al. (2023). The scheme depends on
#' (i) whether the data are at the individual or cluster level, and
#' (ii) whether the target effect is cluster-level or individual-level.
#'
#' Let `J` denote the number of clusters and `N_tot = sum_j N_j` the total
#' number of individuals across clusters. Then:
#'
#' | Data level | target = "clust" | target = "indv"    |
#' |------------|------------------|--------------------|
#' | individual | `alpha = 1 / N_j`  | `alpha = J / N_tot` |
#' | cluster    | `alpha = 1`        | `alpha = J / N_tot * N_j` |
#'
#' Weights sum to `J` so that arithmetic means over rows give cluster-
#' equally-weighted (or individual-equally-weighted) averages.
#'
#' @param data A data frame / data.table with an `id` column (cluster id) and
#'   an `nIndv` column giving cluster size for that row.
#' @param target One of `"clust"` (cluster-level effect) or `"indv"`
#'   (individual-level effect).
#'
#' @return A `data.table` (copy) with an added or overwritten `alpha` column.
#'
#' @export
add_weights <- function(data, target = c("clust", "indv")) {
  target <- match.arg(target)
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  } else {
    data <- data.table::copy(data)
  }
  if (!all(c("id", "nIndv") %in% names(data))) {
    stop("`data` must have columns 'id' and 'nIndv'.")
  }

  data_clust <- data[, .(nIndv = mean(nIndv)), by = id]
  N_tot <- sum(data_clust$nIndv)
  J <- nrow(data_clust)

  is_indv_level <- nrow(data) > J

  data[, alpha := if (is_indv_level) {
    if (target == "clust") 1 / nIndv else J / N_tot
  } else {
    if (target == "clust") 1 else nIndv * J / N_tot
  }]
  data[]
}
