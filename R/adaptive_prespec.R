# Adaptive Pre-specification for Stage 2 TMLE.
#
# Implements Balzer et al. "Adaptive pre-specification in randomized trials
# with and without pair-matching" (Statistics in Medicine 2016).
#
# The procedure:
#   1. Build candidate adjustment specifications (variables + algorithm).
#   2. For each candidate, estimate the CV variance of the resulting TMLE
#      using leave-one-(independent-unit)-out CV (or V-fold CV if there are
#      more than 40 independent units).
#   3. Select the candidate with the smallest CV variance.
#   4. Do the same selection for the propensity score, given the chosen
#      outcome-regression adjustment.
#
# Bug fix vs. from_analysis/aps.R:78 — the original code crashed when the
# user-supplied gAdj was NULL (the typical APS-both path) because
# `gAdj == 'U'` returns logical(0). The fix guards with !is.null.


#' Adaptive Pre-specification: select adjustment for Qbar and pscore
#'
#' @param goal,target,break_match See [stage2_tmle()].
#' @param data Data ready for Stage 2 (column `Y` already scaled to `[0,1]`,
#'   `alpha` set, `id`, `pair` if applicable, plus all candidate columns).
#' @param V Integer. Number of CV folds; ignored when the number of
#'   independent units is ≤ 40 (LOO is used instead).
#' @param cand_q_adj,cand_q_form,cand_g_adj,cand_g_form Candidate
#'   adjustment specifications as built by [cand_adj()].
#' @param remove_pscore If `TRUE`, the variable selected for `q_adj` is
#'   removed from the propensity-score candidates. Useful for very small
#'   trials when `cand_q_adj == cand_g_adj` and only main-terms GLMs are
#'   considered. See `from_analysis/aps.R` for the original implementation.
#' @param q_adj,g_adj Pre-specified adjustment sets. If supplied, the
#'   corresponding selection step is skipped.
#' @param scale_min,scale_max Outcome-rescaling parameters.
#' @param verbose Logical.
#'
#' @return A list with components `q_index`, `q_adj`, `q_form`, `g_index`,
#'   `g_adj`, `g_form`, `var_cv` (CV variance of the chosen TMLE).
#'
#' @keywords internal
adaptive_prespec <- function(goal, target = "clust", break_match = TRUE,
                             data, V = 5,
                             cand_q_adj, cand_q_form,
                             cand_g_adj, cand_g_form,
                             remove_pscore = FALSE,
                             q_adj = NULL, g_adj = NULL,
                             scale_min, scale_max, verbose = FALSE) {

  data <- as_dt_copy(data)
  data[, indpt_unit := if (break_match) id else pair]
  unique_ids <- unique(data$indpt_unit)

  if (length(unique_ids) > 40) {
    folds <- get_folds(V = V, Y = data$Y, ids = unique_ids)
  } else {
    folds <- lapply(unique_ids, function(u) u)
  }

  q_index <- NA_integer_
  q_form_sel <- NA_character_
  g_index <- NA_integer_
  g_form_sel <- NA_character_
  var_cv <- NA_real_

  if (is.null(q_adj)) {
    select_q <- suppressWarnings(cv_selector(
      goal = goal, target = target, break_match = break_match, data = data,
      cand_adj = cand_q_adj, cand_form = cand_q_form, for_q = TRUE,
      scale_min = scale_min, scale_max = scale_max, folds = folds
    ))
    if (verbose) message("APS Q-selection: ", paste(unlist(select_q$adj), collapse = ","),
                          " (", select_q$form, ")")
    q_index <- select_q$adj_index
    q_adj <- select_q$adj
    q_form_sel <- select_q$form

    # Short-circuit: if the unadjusted estimator is chosen for Qbar (or the
    # user pre-specified the unadjusted pscore), skip pscore selection.
    if (sum(q_adj == "U") || (!is.null(g_adj) && all(g_adj == "U"))) {
      g_index <- -99L
      g_adj <- "U"
      g_form_sel <- "glm"
      var_cv <- select_q$var_cv
    } else if (remove_pscore) {
      # Drop the Q-selected variable from pscore candidates if the candidate
      # sets are identical and gform is glm-only.
      if (identical(unlist(cand_g_adj), unlist(cand_q_adj)) &&
          all(cand_g_form == "glm")) {
        cand_g_adj <- cand_g_adj[-q_index]
        cand_g_form <- cand_g_form[-q_index]
      }
    }
  }

  if (is.null(g_adj)) {
    select_g <- suppressWarnings(cv_selector(
      goal = goal, target = target, break_match = break_match, data = data,
      cand_adj = cand_g_adj, cand_form = cand_g_form, for_q = FALSE,
      q_adj = q_adj, q_form = q_form_sel,
      scale_min = scale_min, scale_max = scale_max, folds = folds
    ))
    if (verbose) message("APS g-selection: ", paste(unlist(select_g$adj), collapse = ","),
                          " (", select_g$form, ")")
    g_index <- select_g$adj_index
    g_adj <- select_g$adj
    g_form_sel <- select_g$form
    var_cv <- select_g$var_cv
  }

  list(q_index = q_index, q_adj = q_adj, q_form = q_form_sel,
       g_index = g_index, g_adj = g_adj, g_form = g_form_sel,
       var_cv = var_cv)
}


#' CV selector: pick the candidate with the smallest CV risk
#'
#' @keywords internal
cv_selector <- function(goal, target, break_match, data,
                        cand_adj, cand_form, for_q,
                        q_adj = NULL, q_form = NULL,
                        scale_min, scale_max, folds) {

  if (length(cand_form) == 1L) {
    cand_form <- rep(cand_form, length(cand_adj))
  }
  if (length(cand_form) != length(cand_adj)) {
    stopf("Candidate-form vector length (%d) does not match candidate-adj list length (%d).",
          length(cand_form), length(cand_adj))
  }

  num_tmles <- length(cand_form)
  cv_risk <- numeric(num_tmles)
  var_cv <- numeric(num_tmles)
  for (k in seq_len(num_tmles)) {
    if (for_q) {
      ic <- get_ic_cv(
        goal = goal, target = target, break_match = break_match, data = data,
        q_adj = cand_adj[[k]], q_form = cand_form[k],
        g_adj = NULL, g_form = "glm",
        scale_min = scale_min, scale_max = scale_max, folds = folds
      )
    } else {
      ic <- get_ic_cv(
        goal = goal, target = target, break_match = break_match, data = data,
        q_adj = q_adj, q_form = q_form,
        g_adj = cand_adj[[k]], g_form = cand_form[k],
        scale_min = scale_min, scale_max = scale_max, folds = folds
      )
    }
    cv_risk[k] <- ic$cv_risk
    var_cv[k] <- ic$var_cv
  }

  adj_index <- which.min(cv_risk)
  list(cv_risk = cv_risk, adj_index = adj_index,
       adj = cand_adj[[adj_index]], form = cand_form[adj_index],
       var_cv = var_cv[adj_index])
}


#' Cross-validated influence-curve estimate for a single candidate
#'
#' For each fold, fit the TMLE on the training portion and compute the IC on
#' the validation portion using the *training* fits (no re-fit on validation).
#' Returns the average CV risk (mean of IC^2) and CV variance.
#'
#' @keywords internal
get_ic_cv <- function(goal, target, break_match, data,
                      q_adj, q_form, g_adj = NULL, g_form = "glm",
                      scale_min, scale_max, folds) {

  n_folds <- length(folds)
  DY_cv <- numeric()
  cv_risk_acc <- numeric()
  for (i in seq_len(n_folds)) {
    in_valid <- data$indpt_unit %in% folds[[i]]
    valid <- data[in_valid]
    train <- data[!in_valid]

    train_out <- do_tmle(
      train = train, goal = goal, target = target,
      q_adj = q_adj, q_form = q_form,
      g_adj = g_adj, g_form = g_form,
      scale_min = scale_min, scale_max = scale_max,
      doing_cv = TRUE, verbose = FALSE
    )

    valid_out <- do_tmle_validset(
      goal = goal, target = target, valid = valid, train_out = train_out,
      scale_min = scale_min, scale_max = scale_max
    )

    if (break_match) {
      DY_cv <- c(DY_cv, valid_out$DY)
      cv_risk_acc <- c(cv_risk_acc, mean(valid_out$DY^2))
    } else {
      DY_cv <- c(DY_cv, valid_out$DY_paired)
      cv_risk_acc <- c(cv_risk_acc, mean(valid_out$DY_paired^2))
    }
  }
  list(cv_risk = mean(cv_risk_acc), var_cv = stats::var(DY_cv) / length(DY_cv))
}


#' Evaluate TMLE fits from training on a held-out validation fold
#'
#' Used inside [get_ic_cv()]. Applies the training Q.out and p.out to the
#' validation rows and computes the IC.
#'
#' @keywords internal
do_tmle_validset <- function(goal, target, valid, train_out,
                             scale_min, scale_max) {
  q <- do_init_qbar(train = valid, q_adj = train_out$q_adj,
                    q_form = train_out$q_form, q_out = train_out$q_out)
  valid <- q$train
  g <- get_clever_cov(train = valid, g_adj = train_out$g_adj,
                      g_form = train_out$g_form, p_out = train_out$p_out)
  valid <- g$train
  valid <- do_targeting(train = valid, eps = train_out$eps, goal = goal)
  get_ic_variance(
    goal = goal, target = target, vdata = valid,
    R1 = train_out$R1, R0 = train_out$R0,
    scale_min = scale_min, scale_max = scale_max,
    doing_cv = TRUE
  )
}
