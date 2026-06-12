# Layered Stage-1 API: make easy things easy, hard things possible.
#   stage1_apply  -- general engine: map a per-cluster estimator over clusters
#   stage1_mean   -- preset: one measurement-adjusted cluster mean
#   stage1_ratio  -- preset: ratio of two cluster means
# stage1_mean / stage1_ratio are thin presets over stage1_apply; all three
# return the canonical cluster-level table consumed by add_weights()/stage2_tmle().


#' Stage 1: apply a per-cluster estimator and assemble the cluster table
#'
#' General Stage-1 engine. For each cluster (a distinct value of `id_col`), calls
#' `estimator_fun()` on that cluster's rows and collects the result into the
#' canonical cluster-level `data.table` consumed by [add_weights()] and
#' [stage2_tmle()]. Use this when neither [stage1_mean()] nor [stage1_ratio()]
#' fits — e.g. composite or longitudinal endpoints — by supplying your own
#' `estimator_fun` (which may call [ltmle::ltmle()] directly).
#'
#' @param data Individual-level data, one row per participant. Must contain
#'   `id_col`, `arm_col`, the columns named in `clust_cov`, and (optionally) a
#'   `pair` column for pair-matched designs.
#' @param estimator_fun Function of one argument (a cluster's rows, as a
#'   `data.table`) returning either a single number or a named numeric vector of
#'   one-or-more cluster-level endpoints (e.g. `c(unadj = ., tmle = .)`).
#' @param id_col Name of the cluster id column. Default `"id"`.
#' @param arm_col Name of the (cluster-constant) treatment column. Default
#'   `"A"`. Carried to the output as `A`.
#' @param clust_cov Character vector of cluster-level covariate columns to carry
#'   through to the cluster table (assumed constant within cluster).
#' @param verbose Logical.
#'
#' @return A `data.table` with one row per cluster: `id`, `nIndv`, `U` = 1, `A`,
#'   `pair` (if present), the `clust_cov` columns, plus the estimator's
#'   column(s). A single-number estimator yields a column named `est`.
#'
#' @seealso [stage1_mean()], [stage1_ratio()]
#' @export
stage1_apply <- function(data, estimator_fun, id_col = "id", arm_col = "A",
                         clust_cov = NULL, verbose = FALSE) {
  data <- as_dt_copy(data)
  if (!id_col %in% names(data)) stopf("`id_col` '%s' is not a column of `data`.", id_col)
  if (!arm_col %in% names(data)) stopf("`arm_col` '%s' is not a column of `data`.", arm_col)
  missing_cc <- setdiff(clust_cov, names(data))
  if (length(missing_cc)) stopf("clust_cov not in data: %s", paste(missing_cc, collapse = ", "))

  has_pair <- "pair" %in% names(data)
  clusters <- sort(unique(data[[id_col]]))
  n_clust  <- length(clusters)

  out <- vector("list", n_clust)
  for (j in seq_len(n_clust)) {
    oc  <- data[data[[id_col]] == clusters[j]]
    est <- estimator_fun(oc)
    if (is.null(names(est))) {
      if (length(est) != 1L) stopf("estimator_fun must return a named vector or a single value.")
      est <- c(est = unname(est))
    }
    meta <- list(id = clusters[j], nIndv = nrow(oc), U = 1L, A = oc[[arm_col]][1L])
    if (has_pair) meta$pair <- oc[["pair"]][1L]
    for (cc in clust_cov) meta[[cc]] <- oc[[cc]][1L]
    out[[j]] <- data.table::as.data.table(c(meta, as.list(est)))
    vmsg(verbose, "Stage 1 cluster %d/%d", j, n_clust)
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}


#' Stage 1: one measurement-adjusted cluster mean
#'
#' For each cluster, estimates \eqn{E[Y]} adjusting for differential measurement
#' on `adj_covars` via a censored-mean TMLE (see [fit_censored_mean()]). Returns
#' the canonical cluster table with a single endpoint column `est`. For the
#' **unadjusted** estimate, call with `adj_covars = NULL` (the default).
#'
#' `Y` may contain `NA` for unmeasured participants; missingness is modelled as
#' the censoring node. A fully-observed `Y` (no `NA`) reduces to its plain
#' `mean()`.
#'
#' @param data Individual-level data (one row per participant).
#' @param y_col Name of the outcome column (binary or bounded; `NA` = unmeasured).
#' @param id_col,arm_col,clust_cov,verbose See [stage1_apply()].
#' @param adj_covars Adjustment covariates: a character vector of column names,
#'   `NULL` (default; unadjusted), or a function `(cluster_data) -> character`
#'   resolved per cluster (e.g. to shrink the adjustment set for small clusters).
#' @param det_q_fun,det_cols Optional deterministic-Q closure and the column(s)
#'   it references (added to the modelling frame). See [fit_censored_mean()].
#' @param sl_library,variance_method Passed to [fit_censored_mean()].
#'
#' @return A `data.table`; see [stage1_apply()]. Endpoint column `est`.
#' @seealso [stage1_ratio()], [stage1_apply()]
#' @export
stage1_mean <- function(data, y_col, id_col = "id", adj_covars = NULL,
                        det_q_fun = NULL, det_cols = NULL,
                        sl_library = "glm", variance_method = "tmle",
                        arm_col = "A", clust_cov = NULL, verbose = FALSE) {
  stage1_apply(
    data,
    function(oc) .censored_mean(oc, y_col, adj_covars, det_q_fun, det_cols,
                                sl_library, variance_method),
    id_col = id_col, arm_col = arm_col, clust_cov = clust_cov, verbose = verbose)
}


#' Stage 1: ratio of two cluster means
#'
#' For each cluster, estimates `numerator / denominator`, where each side is a
#' measurement-adjusted cluster mean (a censored-mean TMLE; see
#' [fit_censored_mean()]). Returns the canonical cluster table with a single
#' endpoint column `est`. For the **unadjusted** ratio, call with
#' `num_adj = NULL, den_adj = NULL` (the defaults).
#'
#' This covers both the coverage/yield shape (a fully-observed numerator over a
#' partially-missing, adjusted denominator) and the care-cascade shape (both
#' sides partially missing and adjusted). A side whose column has **no** `NA`
#' reduces to its plain `mean()` (with no censoring the censored-mean TMLE
#' returns the empirical mean), so pre-code a fully-observed numerator as desired
#' (e.g. `NA` recoded to 0).
#'
#' @param data Individual-level data (one row per participant).
#' @param num_col,den_col Numerator and denominator column names (`NA` =
#'   unmeasured on that side).
#' @param num_adj,den_adj Adjustment covariates for each side: a character
#'   vector, `NULL` (unadjusted), or a function `(cluster_data) -> character`
#'   resolved per cluster (e.g. to shrink the set for small clusters).
#' @param num_det_q,den_det_q,num_det_cols,den_det_cols Optional deterministic-Q
#'   closures and their referenced columns, per side.
#' @param sl_library,variance_method Passed to [fit_censored_mean()].
#' @param id_col,arm_col,clust_cov,verbose See [stage1_apply()].
#'
#' @return A `data.table`; see [stage1_apply()]. Endpoint column `est`.
#' @seealso [stage1_mean()], [stage1_apply()]
#' @export
stage1_ratio <- function(data, num_col, den_col, id_col = "id",
                         num_adj = NULL, den_adj = NULL,
                         num_det_q = NULL, den_det_q = NULL,
                         num_det_cols = NULL, den_det_cols = NULL,
                         sl_library = "glm", variance_method = "tmle",
                         arm_col = "A", clust_cov = NULL, verbose = FALSE) {
  est_fun <- function(oc) {
    num <- .censored_mean(oc, num_col, num_adj, num_det_q, num_det_cols,
                          sl_library, variance_method)
    den <- .censored_mean(oc, den_col, den_adj, den_det_q, den_det_cols,
                          sl_library, variance_method)
    num / den
  }
  stage1_apply(data, est_fun, id_col = id_col, arm_col = arm_col,
               clust_cov = clust_cov, verbose = verbose)
}


# Internal: censored-mean estimate of E[Y] for one cluster's column, adjusting
# for `adj` (NULL = unadjusted). Builds the [adj..., det_cols..., C, Y] frame and
# calls fit_censored_mean(). `adj` may be a character vector, NULL, or a
# function(cluster_data)->character (resolved per cluster). A fully-observed
# column with no deterministic-Q short-circuits to the exact mean (with no
# censoring the censored-mean TMLE returns the empirical mean for any adjustment
# set; this is exact and skips a needless ltmle fit).
.censored_mean <- function(oc, col, adj, det_q, det_cols,
                           sl_library, variance_method) {
  if (is.function(adj)) adj <- adj(oc)
  y <- oc[[col]]
  if (is.null(det_q) && !anyNA(y)) return(mean(y))
  base <- data.frame(C = ltmle::BinaryToCensoring(is.censored = is.na(y)), Y = y)
  if (!is.null(det_cols)) base <- cbind(as.data.frame(oc[, det_cols, with = FALSE]), base)
  frame <- if (is.null(adj)) base else cbind(as.data.frame(oc[, adj, with = FALSE]), base)
  fit_censored_mean(frame, sl_library, det_q, variance_method)
}
