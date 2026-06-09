#' Two-Stage TMLE: end-to-end driver
#'
#' Convenience wrapper that runs [stage1_tmle()] then [stage2_tmle()] on
#' individual-level observed data, returning the Stage 2 effect estimate
#' together with the intermediate cluster-level dataset.
#'
#' @param data Individual-level observed data. See [stage1_tmle()] for the
#'   required columns.
#' @param ind_cov,clust_cov,sl_library See [stage1_tmle()].
#' @param ... Additional arguments forwarded to [stage2_tmle()] (e.g.,
#'   `goal`, `target`, `break_match`, `do_data_adapt`, `cand_q_adj`).
#' @param target Stage 2 target effect (`"clust"` or `"indv"`).
#' @param verbose Logical.
#'
#' @return A list with components:
#'   * `est` — the Stage 2 estimate `data.table`
#'   * `stage1` — the Stage 1 cluster-level dataset
#'
#' @details Internal reference implementation of the Balzer et al. (2021)
#'   end-to-end driver. Not part of the public API — compose the workflow from
#'   [stage1_apply()] + [add_weights()] + [stage2_tmle()] (see the
#'   `two-stage-tmle` vignette).
#' @keywords internal
two_stage_tmle <- function(data, ind_cov, clust_cov,
                           sl_library = c("SL.mean", "SL.glm", "SL.gam"),
                           target = c("clust", "indv"),
                           ...,
                           verbose = FALSE) {
  target <- match.arg(target)

  stage1 <- stage1_tmle(
    data = data, ind_cov = ind_cov, clust_cov = clust_cov,
    sl_library = sl_library, verbose = verbose
  )

  # Default Stage 2 input: take Stage 1's TMLE cluster endpoint as Y.
  stage2_in <- data.table::copy(stage1)
  stage2_in[, Y := tmle]
  stage2_in <- add_weights(stage2_in, target = target)

  est <- stage2_tmle(data = stage2_in, target = target, verbose = verbose, ...)

  list(est = est, stage1 = stage1)
}
