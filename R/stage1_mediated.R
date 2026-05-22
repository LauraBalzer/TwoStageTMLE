#' Stage 1: ratio-based cluster endpoint for the multilevel-mediation-missing
#' data problem
#'
#' For each cluster separately, returns the cluster-specific endpoint
#' \eqn{Y^c = \mathbb{P}(Y_2 = 1) / \mathbb{E}[\mathbb{P}(Y_1 = 1 \mid M = 1, W)]}
#' from Nakato et al. (2026). This is the Stage 1 estimator for cluster
#' randomized trials in which the *eligibility / risk-status* variable `Y1`
#' is partially missing because it is only observed among those screened
#' (`M = 1`). The outcome `Y2` is assumed completely measured.
#'
#' Within each cluster:
#' * **Numerator** — empirical \eqn{\mathbb{P}(Y_2 = 1)}.
#' * **Denominator** — an individual-level TMLE for
#'   \eqn{\mathbb{E}[\mathbb{P}(Y_1 = 1 \mid M = 1, W)]} fit via
#'   [ltmle::ltmle()] with `Anodes = M`, `Ynodes = Y1`, `abar = 1`,
#'   `stratify = TRUE`, and the user-supplied Super Learner library.
#' * **Endpoint** — numerator divided by denominator.
#'
#' The returned cluster-level dataset can be passed directly to
#' [stage2_tmle()] (or via the [two_stage_tmle_mediated()] driver).
#'
#' Stage 2 treats this endpoint as a fixed observation at the cluster level;
#' the Stage 1 contribution to overall variance is assumed asymptotically
#' negligible when within-cluster sample sizes are large. See Box 2 of
#' Balzer et al. for the formal conditions.
#'
#' @param data Individual-level observed data, one row per participant. Must
#'   contain columns:
#'   * `id` — cluster id
#'   * `A` — cluster-level treatment indicator (0/1)
#'   * `M` — screening indicator (1 if screened, 0 if missed)
#'   * `Y1` — observed risk status (binary; by the `Y1 = M * Y1_star`
#'     convention `Y1 = 0` whenever `M = 0`)
#'   * `Y2` — outcome of interest, assumed fully observed (binary)
#'   * the columns named in `ind_cov`
#'   * the columns named in `clust_cov`
#'   * `pair` (cluster pair id) if pair-matched
#' @param ind_cov Character vector of individual-level covariate names used to
#'   fit the denominator. Per Box 1 of Balzer et al., this is `W` plus any
#'   time-varying covariates `L`.
#' @param clust_cov Character vector of cluster-level covariates. Passed
#'   through to the returned cluster-level dataset; not used in the Stage 1
#'   fit (Stage 1 stratifies on cluster so cluster covariates are constant
#'   within a fit).
#' @param sl_library Super Learner library passed to [ltmle::ltmle()].
#'   Default `c("SL.mean", "SL.glm", "SL.gam")`, matching the simulation
#'   library used in Nakato et al. (2026).
#' @param verbose Logical.
#'
#' @return A `data.table` with one row per cluster and columns
#'   `id`, `pair` (if present in input), `nIndv` (cluster size), `U` = 1,
#'   `A`, all `clust_cov`, plus:
#'   * `screened` — \eqn{\mathbb{P}(Y_2 = 1 \mid M = 1)} (biased; comparator only)
#'   * `eligible` — \eqn{\mathbb{P}(Y_2 = 1 \mid Y_1 = 1)} (biased; comparator only)
#'   * `unadj` — ratio with empirical denominator \eqn{\mathbb{P}(Y_1 = 1 \mid M = 1)}
#'   * `tmle` — ratio with TMLE-adjusted denominator (the recommended endpoint)
#'
#' @references
#' Nakato J., Balzer L.B. (2026). Methods for evaluating overall effectiveness
#' of cluster-level strategies aiming to improve reach and health outcomes.
#'
#' Balzer L.B. et al. (2021). Two-Stage TMLE to Reduce Bias and Improve
#' Efficiency in Cluster Randomized Trials. arXiv:2106.15737.
#'
#' @seealso [stage1_tmle()] for the missing-outcomes-only version;
#'   [two_stage_tmle_mediated()] for the end-to-end driver.
#'
#' @export
stage1_mediated <- function(data, ind_cov, clust_cov,
                            sl_library = c("SL.mean", "SL.glm", "SL.gam"),
                            verbose = FALSE) {
  data <- as_dt_copy(data)
  validate_stage1_mediated_input(data, ind_cov, clust_cov)

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

    # Numerator: empirical P(Y2 = 1)
    numerator <- mean(oc$Y2)

    # Biased comparators (returned for downstream sensitivity analyses)
    screened <- mean(oc[M == 1, Y2])
    eligible <- mean(oc[Y1 == 1, Y2])

    # Unadjusted denominator: empirical P(Y1 = 1 | M = 1)
    den_unadj <- mean(oc[M == 1, Y1])

    # TMLE denominator: P*(Y1 = 1 | M = 1, W)
    oc_sub <- as.data.frame(oc[, c(ind_cov, "M", "Y1"), with = FALSE])
    suppressMessages(suppressWarnings({
      den_fit <- ltmle::ltmle(
        data = oc_sub, Anodes = "M", Ynodes = "Y1", abar = 1,
        SL.library = sl_library, estimate.time = FALSE, stratify = TRUE
      )
    }))
    den_tmle <- unname(den_fit$estimates["tmle"])

    meta <- as.list(oc[1L, intersect(clust_meta_cols, names(oc)), with = FALSE])
    meta$nIndv <- nIndv_j
    meta$U <- 1L
    out[[j]] <- data.table::as.data.table(c(
      meta,
      list(
        screened = screened,
        eligible = eligible,
        unadj    = numerator / den_unadj,
        tmle     = numerator / den_tmle
      )
    ))
    vmsg(verbose, "Stage 1 (mediated) cluster %d/%d", j, n_clust)
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}


#' @keywords internal
validate_stage1_mediated_input <- function(data, ind_cov, clust_cov) {
  required <- c("id", "A", "M", "Y1", "Y2", ind_cov, clust_cov)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stopf("Stage 1 (mediated) input is missing required columns: %s.",
          paste(missing_cols, collapse = ", "))
  }
  if (!all(data$A %in% c(0, 1))) {
    stopf("Column 'A' must be binary (0/1).")
  }
  if (!all(data$M %in% c(0, 1))) {
    stopf("Column 'M' must be binary (0/1).")
  }
  if (!all(data$Y1 %in% c(0, 1))) {
    stopf(paste0(
      "Column 'Y1' must be binary (0/1). Use the Y1 = M * Y1_star convention ",
      "so that Y1 = 0 whenever M = 0."
    ))
  }
  if (!all(data$Y2 %in% c(0, 1))) {
    stopf("Column 'Y2' must be binary (0/1).")
  }
  invisible(TRUE)
}
