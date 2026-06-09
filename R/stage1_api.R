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
#' the canonical cluster table with two endpoint columns:
#' * `unadj` — censored mean with no covariate adjustment, and
#' * `tmle`  — censored mean adjusting for `adj_covars` (equals `unadj` when
#'   `adj_covars` is `NULL`).
#'
#' `Y` may contain `NA` for unmeasured participants; missingness is modelled as
#' the censoring node.
#'
#' @param data Individual-level data (one row per participant).
#' @param y_col Name of the outcome column (binary or bounded; `NA` = unmeasured).
#' @param id_col,arm_col,clust_cov,verbose See [stage1_apply()].
#' @param adj_covars Character vector of individual-level adjustment covariates
#'   for the `tmle` endpoint. `NULL` (default) gives an unadjusted estimate.
#' @param det_q_fun,det_cols Optional deterministic-Q closure and the column(s)
#'   it references (added to the modelling frame). See [fit_censored_mean()].
#' @param sl_library,variance_method Passed to [fit_censored_mean()].
#'
#' @return A `data.table`; see [stage1_apply()]. Endpoint columns `unadj`, `tmle`.
#' @seealso [stage1_ratio()], [stage1_apply()]
#' @export
stage1_mean <- function(data, y_col, id_col = "id", adj_covars = NULL,
                        det_q_fun = NULL, det_cols = NULL,
                        sl_library = "glm", variance_method = "tmle",
                        arm_col = "A", clust_cov = NULL, verbose = FALSE) {
  est_fun <- function(oc) {
    out <- .censored_mean_pair(oc, y_col, adj_covars, det_q_fun, det_cols,
                               sl_library, variance_method)
    c(unadj = out[["unadj"]], tmle = out[["tmle"]])
  }
  stage1_apply(data, est_fun, id_col = id_col, arm_col = arm_col,
               clust_cov = clust_cov, verbose = verbose)
}


#' Stage 1: ratio of two cluster means
#'
#' For each cluster, estimates `numerator / denominator`, where each side is a
#' cluster mean that may be measurement-adjusted (a censored-mean TMLE) or a
#' plain empirical mean. Returns the canonical cluster table with:
#' * `unadj` — numerator(unadjusted) / denominator(unadjusted), and
#' * `tmle`  — numerator(adjusted)   / denominator(adjusted).
#'
#' Each side is configured independently, so this covers both the coverage/yield
#' shape (plain-mean numerator over an adjusted denominator) and the care-cascade
#' shape (both sides censored and adjusted). When a side is a plain mean
#' (`*_censored = FALSE`), pre-code its column so the mean is the intended
#' numerator/denominator (e.g. `NA` recoded to 0).
#'
#' @param data Individual-level data (one row per participant).
#' @param num_col,den_col Numerator and denominator column names.
#' @param num_censored,den_censored If `TRUE` (default) the side is a censored
#'   mean (`NA` = unmeasured); if `FALSE` it is the plain `mean()` of the column.
#' @param num_adj,den_adj Adjustment covariates for each side's `tmle` endpoint.
#' @param num_det_q,den_det_q,num_det_cols,den_det_cols Optional deterministic-Q
#'   closures and their referenced columns, per side.
#' @param sl_library,variance_method Passed to [fit_censored_mean()].
#' @param id_col,arm_col,clust_cov,verbose See [stage1_apply()].
#'
#' @return A `data.table`; see [stage1_apply()]. Endpoint columns `unadj`, `tmle`.
#' @seealso [stage1_mean()], [stage1_apply()]
#' @export
stage1_ratio <- function(data, num_col, den_col, id_col = "id",
                         num_adj = NULL, den_adj = NULL,
                         num_det_q = NULL, den_det_q = NULL,
                         num_det_cols = NULL, den_det_cols = NULL,
                         num_censored = TRUE, den_censored = TRUE,
                         sl_library = "glm", variance_method = "tmle",
                         arm_col = "A", clust_cov = NULL, verbose = FALSE) {
  side <- function(oc, col, adj, det_q, det_cols, censored) {
    if (!censored) {
      m <- mean(oc[[col]])
      return(list(unadj = m, tmle = m))
    }
    .censored_mean_pair(oc, col, adj, det_q, det_cols, sl_library, variance_method)
  }
  est_fun <- function(oc) {
    num <- side(oc, num_col, num_adj, num_det_q, num_det_cols, num_censored)
    den <- side(oc, den_col, den_adj, den_det_q, den_det_cols, den_censored)
    c(unadj = num[["unadj"]] / den[["unadj"]],
      tmle  = num[["tmle"]]  / den[["tmle"]])
  }
  stage1_apply(data, est_fun, id_col = id_col, arm_col = arm_col,
               clust_cov = clust_cov, verbose = verbose)
}


# Internal: estimate the (unadjusted, adjusted) censored-mean pair for one
# cluster's column. Builds the [det_cols..., C, Y] frame and calls
# fit_censored_mean() once with no covariates and (if adj given) once with them.
.censored_mean_pair <- function(oc, col, adj, det_q, det_cols,
                                sl_library, variance_method) {
  y     <- oc[[col]]
  Cnode <- ltmle::BinaryToCensoring(is.censored = is.na(y))
  base  <- data.frame(C = Cnode, Y = y)
  if (!is.null(det_cols)) {
    base <- cbind(as.data.frame(oc[, det_cols, with = FALSE]), base)
  }
  unadj <- fit_censored_mean(base, sl_library, det_q, variance_method)
  if (is.null(adj)) {
    tmle <- unadj
  } else {
    adjf <- cbind(as.data.frame(oc[, adj, with = FALSE]), base)
    tmle <- fit_censored_mean(adjf, sl_library, det_q, variance_method)
  }
  list(unadj = unadj, tmle = tmle)
}
