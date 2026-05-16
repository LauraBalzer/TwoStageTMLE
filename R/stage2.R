#' Stage 2: cluster-level TMLE for the intervention effect
#'
#' Estimates the intervention effect given a cluster-level (or
#' individual-level) dataset whose outcome has already been adjusted for
#' missingness (e.g., by [stage1_tmle()] in a Two-Stage workflow, or by the
#' user directly). Supports risk difference, arithmetic risk ratio, and odds
#' ratio targets; pair-matched and unmatched designs; cluster-level and
#' individual-level target effects.
#'
#' The TMLE is implemented as in § 9.4 of Balzer et al. (2021): an initial
#' outcome regression updated with a two-dimensional (or one-dimensional for
#' RD) clever covariate. Adaptive Pre-specification (Balzer et al. 2016) is
#' available via `do_data_adapt = TRUE` and selects the adjustment set that
#' minimises the CV variance of the resulting TMLE.
#'
#' @param data A `data.frame` / `data.table` with at minimum columns
#'   `id` (cluster id), `A` (treatment indicator, 0/1), `Y` (outcome),
#'   `U` (dummy `= 1` for the unadjusted estimator), and `alpha` (analysis
#'   weights — use [add_weights()] if unsure). If the design is pair-matched
#'   and `break_match = FALSE`, also include a `pair` column. If the data are
#'   at the individual level, an `nIndv` column with cluster sizes is required.
#' @param goal One of `"RD"` (risk difference), `"aRR"` (arithmetic risk
#'   ratio), `"OR"` (odds ratio).
#' @param target One of `"clust"` (cluster-level effect) or `"indv"`
#'   (individual-level effect).
#' @param unadjusted If `TRUE`, fits the unadjusted treatment-specific
#'   means directly without running the TMLE update. Faster and more stable
#'   for rare outcomes. Should not be used with `target = "indv"` on
#'   cluster-level data because weights are not incorporated.
#' @param q_adj,q_form,g_adj,g_form Pre-specified adjustment sets / algorithms.
#'   If `do_data_adapt = TRUE` they default to `NULL` and are chosen by APS.
#' @param do_data_adapt If `TRUE`, run Adaptive Pre-specification to select
#'   `q_adj` / `g_adj` from the candidate sets `cand_q_adj` / `cand_g_adj`.
#' @param cand_q_adj,cand_q_form,cand_g_adj,cand_g_form Candidate adjustment
#'   specifications. Use [cand_adj()] to build them.
#' @param V Integer. CV folds for APS. Defaults to `5`; ignored if there are
#'   fewer than 41 independent units (LOO is used).
#' @param remove_pscore If `TRUE`, drop the variable selected for `q_adj`
#'   from the propensity-score candidates. Useful when candidate sets are
#'   identical and main-terms GLMs are the only algorithm.
#' @param do_cv_variance If `TRUE`, also report CV-based inference.
#' @param break_match If `TRUE` (default), treat clusters as independent
#'   units; otherwise treat pairs as independent units (variance over the
#'   N/2 pair sums; t-distribution with N/2 - 1 df).
#' @param one_sided If `TRUE`, the reported p-value is one-sided. Must
#'   supply `alt_smaller`.
#' @param alt_smaller Direction of the one-sided alternative. `TRUE` means
#'   the intervention is hypothesised to *decrease* the mean outcome.
#' @param verbose Logical.
#' @param psi True parameter value, if known (e.g., in simulation).
#' @param return_ic If `TRUE`, return the full TMLE fit object alongside the
#'   estimate data.frame.
#'
#' @return A one-row `data.table` with point estimate, 95% CI, SE, p-value,
#'   bias, coverage, rejection indicator, and the chosen adjustment indices
#'   for Q and g. When `return_ic = TRUE`, the return value is a list with
#'   the data.table in `est` and the full fit in `ic`.
#'
#' @references
#' Balzer L.B., van der Laan M., Ayieko J., Kamya M., Chamie G., Schwab J.,
#' Havlir D., Petersen M. (2021). Two-Stage TMLE to Reduce Bias and Improve
#' Efficiency in Cluster Randomized Trials. arXiv:2106.15737.
#'
#' Balzer L.B., van der Laan M., Petersen M., SEARCH Consortium (2016).
#' Adaptive pre-specification in randomized trials with and without
#' pair-matching. *Statistics in Medicine* 35(10):4528-4545.
#'
#' @export
stage2_tmle <- function(data,
                        goal = c("aRR", "RD", "OR"),
                        target = c("clust", "indv"),
                        unadjusted = FALSE,
                        q_adj = NULL, q_form = "glm",
                        g_adj = NULL, g_form = "glm",
                        do_data_adapt = FALSE,
                        cand_q_adj = NULL, cand_q_form = "glm",
                        cand_g_adj = NULL, cand_g_form = "glm",
                        V = 5, remove_pscore = FALSE,
                        do_cv_variance = FALSE,
                        break_match = TRUE,
                        one_sided = FALSE, alt_smaller = NULL,
                        verbose = FALSE, psi = NA_real_,
                        return_ic = FALSE) {

  goal <- match.arg(goal)
  target <- match.arg(target)
  data <- as_dt_copy(data)

  validate_stage2_input(data, break_match)

  # Rescale Y to `[0, 1]` (no-op for already-bounded outcomes).
  rescale <- scale_outcome(data$Y)
  data[, Y := rescale$y]
  scale_min <- rescale$scale_min
  scale_max <- rescale$scale_max

  if (do_data_adapt) {
    select <- adaptive_prespec(
      goal = goal, target = target, break_match = break_match,
      data = data, V = V,
      cand_q_adj = cand_q_adj, cand_q_form = cand_q_form,
      cand_g_adj = cand_g_adj, cand_g_form = cand_g_form,
      remove_pscore = remove_pscore,
      q_adj = q_adj, g_adj = g_adj,
      scale_min = scale_min, scale_max = scale_max, verbose = verbose
    )
    q_index <- select$q_index; q_adj <- select$q_adj; q_form <- select$q_form
    g_index <- select$g_index; g_adj <- select$g_adj; g_form <- select$g_form
  } else {
    q_index <- 1L; g_index <- 1L
    select <- NULL
  }

  if (unadjusted) {
    est <- do_unadjusted(train = data, goal = goal,
                         scale_min = scale_min, scale_max = scale_max,
                         verbose = verbose)
  } else {
    est <- do_tmle(
      train = data, goal = goal, target = target,
      q_adj = q_adj, q_form = q_form,
      g_adj = g_adj, g_form = g_form,
      scale_min = scale_min, scale_max = scale_max,
      doing_cv = FALSE, verbose = verbose
    )
  }

  n_clust <- length(unique(data$id))
  vmsg(verbose, "Independent units: %d", n_clust)

  R1 <- est$R1; R0 <- est$R0
  Txt <- get_inference(psi_hat = R1, se = sqrt(est$var_R1),
                       df = n_clust - 2L)[, .(est, ci_lo, ci_hi, se)]
  Con <- get_inference(psi_hat = R0, se = sqrt(est$var_R0),
                       df = n_clust - 2L)[, .(est, ci_lo, ci_hi, se)]

  psi_hat <- switch(goal,
    aRR = log(R1 / R0),
    RD  = R1 - R0,
    OR  = log(R1 / (1 - R1) * (1 - R0) / R0)
  )

  if (break_match) {
    df <- n_clust - 2L
    var_hat <- est$var_break
  } else {
    df <- length(unique(data$pair)) - 1L
    var_hat <- est$var_pair
  }

  inference <- get_inference(goal = goal, psi = psi, psi_hat = psi_hat,
                             se = sqrt(var_hat), df = df,
                             one_sided = one_sided, alt_smaller = alt_smaller)

  est_dt <- data.table::data.table(
    txt_est = Txt$est, txt_ci_lo = Txt$ci_lo, txt_ci_hi = Txt$ci_hi, txt_se = Txt$se,
    con_est = Con$est, con_ci_lo = Con$ci_lo, con_ci_hi = Con$ci_hi, con_se = Con$se,
    psi = psi
  )
  est_dt <- cbind(est_dt, inference)
  # Pull these out before assigning into est_dt to avoid NSE collision with
  # the column `est_dt$est`.
  q_form_val <- est$q_form
  g_form_val <- est$g_form
  est_dt[, `:=`(q_adj = q_index, q_form = q_form_val,
                g_adj = g_index, g_form = g_form_val)]

  if (do_cv_variance && do_data_adapt) {
    inference_cv <- get_inference(
      goal = goal, psi = psi, psi_hat = psi_hat,
      se = sqrt(select$var_cv), df = df,
      one_sided = one_sided, alt_smaller = alt_smaller
    )
    data.table::setnames(inference_cv,
                         old = names(inference_cv),
                         new = paste0("cv_", names(inference_cv)))
    est_dt <- cbind(est_dt, inference_cv)
  }

  if (return_ic) {
    list(est = est_dt, ic = est)
  } else {
    est_dt
  }
}


#' Validate Stage 2 input
#' @keywords internal
validate_stage2_input <- function(data, break_match) {
  required <- c("id", "A", "Y", "U", "alpha")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stopf("Stage 2 input is missing required columns: %s.",
          paste(missing_cols, collapse = ", "))
  }
  if (!all(data$U == 1)) {
    stopf("Column 'U' must equal 1 for every row.")
  }
  if (!all(data$A %in% c(0, 1))) {
    stopf("Column 'A' must be binary (0/1).")
  }
  if (!break_match && !"pair" %in% names(data)) {
    stopf("`break_match = FALSE` requires a `pair` column in `data`.")
  }
  if (nrow(data) > length(unique(data$id)) && !"nIndv" %in% names(data)) {
    stopf("Individual-level data requires an `nIndv` column with cluster sizes.")
  }
  invisible(TRUE)
}
