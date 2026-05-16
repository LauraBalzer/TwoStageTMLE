#' Stage 1: individual-level TMLE for cluster-specific endpoints
#'
#' For each cluster separately, runs an individual-level TMLE to estimate
#' \eqn{Y^c = \mathbb{E}[\mathbb{E}(Y \mid \Delta = 1, W, M)]} — the
#' cluster-specific expected outcome controlling for differential outcome
#' measurement on baseline `W` and post-baseline `M` covariates. Returns a
#' cluster-level dataset suitable for [stage2_tmle()].
#'
#' Implements step 9.2 of Balzer et al. (2021): the individual-level TMLE is
#' fit using [ltmle::ltmle()] with `Anodes = Delta`, `Ynodes = Y`, `abar = 1`,
#' a user-supplied Super Learner library, and `stratify = TRUE`.
#'
#' @param data Individual-level observed data, one row per participant. Must
#'   contain columns:
#'   * `id` — cluster id
#'   * `A` — cluster-level treatment indicator
#'   * `Delta` — outcome measurement indicator (1 if measured, 0 if missing)
#'   * `Y` — outcome (or NA if `Delta = 0`)
#'   * the columns named in `ind_cov`
#'   * the columns named in `clust_cov`
#'   * `pair` (cluster pair id) if pair-matched
#' @param ind_cov Character vector of individual-level covariate names
#'   (baseline plus post-baseline) used by Stage 1 to adjust for differential
#'   measurement. Typically `c("W1", "W2", "M")` in the paper notation.
#' @param clust_cov Character vector of cluster-level covariates. Passed
#'   through to the returned cluster-level dataset; not used in the Stage 1
#'   fit (Stage 1 stratifies on cluster so cluster covariates are constant
#'   within a fit).
#' @param sl_library Super Learner library passed to [ltmle::ltmle()].
#'   Default `c("SL.mean", "SL.glm", "SL.gam")`, matching the paper.
#' @param verbose Logical.
#'
#' @return A `data.table` with one row per cluster and columns
#'   `id`, `pair` (if present in input), `nIndv` (cluster size), `U` = 1,
#'   `A`, all `clust_cov`, plus:
#'   * `unadj` — empirical mean among measured (unadjusted Stage 1 estimator)
#'   * `tmle` — individual-level TMLE estimate of \eqn{Y^c}
#'
#' @references
#' Balzer L.B. et al. (2021). Two-Stage TMLE to Reduce Bias and Improve
#' Efficiency in Cluster Randomized Trials. arXiv:2106.15737.
#'
#' @export
stage1_tmle <- function(data, ind_cov, clust_cov,
                        sl_library = c("SL.mean", "SL.glm", "SL.gam"),
                        verbose = FALSE) {
  data <- as_dt_copy(data)
  validate_stage1_input(data, ind_cov, clust_cov)

  has_pair <- "pair" %in% names(data)
  clusters <- sort(unique(data$id))
  n_clust <- length(clusters)
  clust_meta_cols <- c("id", if (has_pair) "pair" else NULL,
                       "nIndv", "U", "A", clust_cov)

  out <- vector("list", n_clust)
  for (j in seq_len(n_clust)) {
    these <- data$id == clusters[j]
    oc <- data[these]
    nIndv_j <- nrow(oc)

    unadj <- mean(oc[Delta == 1, Y])

    oc_sub <- as.data.frame(oc[, c(ind_cov, "Delta", "Y"), with = FALSE])
    suppressMessages(suppressWarnings({
      yc <- ltmle::ltmle(
        data = oc_sub, Anodes = "Delta", Ynodes = "Y", abar = 1,
        SL.library = sl_library, estimate.time = FALSE, stratify = TRUE
      )
    }))
    tmle <- yc$estimates["tmle"]

    meta <- as.list(oc[1L, intersect(clust_meta_cols, names(oc)), with = FALSE])
    meta$nIndv <- nIndv_j
    meta$U <- 1L
    out[[j]] <- data.table::as.data.table(c(meta, list(unadj = unadj, tmle = unname(tmle))))
    vmsg(verbose, "Stage 1 cluster %d/%d", j, n_clust)
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}


#' @keywords internal
validate_stage1_input <- function(data, ind_cov, clust_cov) {
  required <- c("id", "A", "Delta", "Y", ind_cov, clust_cov)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stopf("Stage 1 input is missing required columns: %s.",
          paste(missing_cols, collapse = ", "))
  }
  if (!all(data$A %in% c(0, 1))) {
    stopf("Column 'A' must be binary (0/1).")
  }
  if (!all(data$Delta %in% c(0, 1))) {
    stopf("Column 'Delta' must be binary (0/1).")
  }
  invisible(TRUE)
}
