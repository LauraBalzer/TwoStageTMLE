#' Two-Stage TMLE for the multilevel-mediation-missing data problem
#'
#' End-to-end driver that runs [stage1_mediated()] then [stage2_tmle()] on
#' individual-level observed data, returning the Stage 2 effect estimate
#' together with the intermediate cluster-level dataset.
#'
#' This is the parallel of [two_stage_tmle()] for the setting where the
#' *eligibility* variable `Y1` is partially missing because it is observed
#' only among those screened (`M = 1`). The cluster-level endpoint is the
#' ratio estimator of Nakato et al. (2026):
#' \deqn{Y^c = \frac{\mathbb{P}(Y_2 = 1)}{\mathbb{E}[\mathbb{P}(Y_1 = 1 \mid M = 1, W)]}.}
#'
#' Stage 2 is the cluster-level TMLE of Balzer et al. (2021) and is
#' unchanged from the foundational driver — the ratio endpoint is handed in
#' as the outcome column. Variance is computed from the Stage 2 cluster-level
#' influence curve and treats the cluster-level endpoint as fixed. This is
#' asymptotically valid when within-cluster sample sizes are large; see
#' Box 2 of Balzer et al. (2021) for the formal conditions.
#'
#' @inheritParams stage1_mediated
#' @param target Stage 2 target effect (`"clust"` or `"indv"`).
#' @param ... Additional arguments forwarded to [stage2_tmle()] (e.g.,
#'   `goal`, `break_match`, `do_data_adapt`, `cand_q_adj`, `cand_g_adj`).
#'
#' @return A list with components:
#'   * `est` — the Stage 2 estimate `data.table`
#'   * `stage1` — the Stage 1 cluster-level dataset (with the four
#'     candidate endpoints; the `tmle` column is used as the Stage 2 outcome)
#'
#' @references
#' Nakato J., Balzer L.B. (2026). Methods for evaluating overall
#' effectiveness of cluster-level strategies aiming to improve reach and
#' health outcomes.
#'
#' Balzer L.B. et al. (2021). Two-Stage TMLE to Reduce Bias and Improve
#' Efficiency in Cluster Randomized Trials. arXiv:2106.15737.
#'
#' @seealso [two_stage_tmle()] for the missing-outcomes-only foundational
#'   driver; [stage1_mediated()], [stage2_tmle()].
#'
#' @export
two_stage_tmle_mediated <- function(data, ind_cov, clust_cov,
                                    sl_library = c("SL.mean", "SL.glm", "SL.gam"),
                                    target = c("clust", "indv"),
                                    ...,
                                    verbose = FALSE) {
  target <- match.arg(target)

  stage1 <- stage1_mediated(
    data = data, ind_cov = ind_cov, clust_cov = clust_cov,
    sl_library = sl_library, verbose = verbose
  )

  # Default Stage 2 input: take Stage 1's TMLE ratio endpoint as Y.
  stage2_in <- data.table::copy(stage1)
  stage2_in[, Y := tmle]
  stage2_in <- add_weights(stage2_in, target = target)

  est <- stage2_tmle(data = stage2_in, target = target, verbose = verbose, ...)

  list(est = est, stage1 = stage1)
}
