#' Build cross-validation folds, optionally stratified on a binary outcome
#'
#' Adapted from `.cvFolds` in the `cvAUC` package (Erin LeDell). With
#' `stratify = TRUE` and binary `Y`, folds preserve the marginal Y=0/Y=1
#' proportions; if either class has fewer than `V` observations, `V` is
#' shrunk to the minimum class size.
#'
#' @param V Integer number of folds requested.
#' @param Y Outcome vector (length equals length of `ids`). Used only when
#'   `stratify = TRUE`.
#' @param ids Vector of independent-unit identifiers (one per row of the
#'   data being split).
#' @param stratify Logical. If `TRUE` and `Y` is binary, stratify folds on
#'   `Y`. Default `FALSE`.
#'
#' @return A list of length `V` (or fewer after shrinkage), each element a
#'   vector of ids comprising one validation fold.
#'
#' @keywords internal
get_folds <- function(V, Y, ids, stratify = FALSE) {
  if (stratify && length(unique(Y)) == 2) {
    classes <- tapply(seq_along(Y), INDEX = Y, FUN = split, 1)
    ids_Y1 <- ids[classes$`1`]
    ids_noY1 <- ids[classes$`0`]
    if (length(ids_Y1) < V || length(ids_noY1) < V) {
      V <- min(length(ids_Y1), length(ids_noY1))
    }
    ids_Y1_split <- split(sample(length(ids_Y1)), rep(seq_len(V), length.out = length(ids_Y1)))
    ids_noY1_split <- split(sample(length(ids_noY1)), rep(seq_len(V), length.out = length(ids_noY1)))
    folds <- vector("list", V)
    for (v in seq_len(V)) {
      folds[[v]] <- c(ids_Y1[ids_Y1_split[[v]]], ids_noY1[ids_noY1_split[[v]]])
    }
  } else {
    ids_split <- split(sample(length(ids)), rep(seq_len(V), length.out = length(ids)))
    folds <- vector("list", V)
    for (v in seq_len(V)) {
      folds[[v]] <- ids[ids_split[[v]]]
    }
  }
  folds
}
