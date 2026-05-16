#' Build candidate adjustment specifications for Adaptive Pre-specification
#'
#' Convenience helper that builds the four arguments
#' (`cand_q_adj`, `cand_q_form`, `cand_g_adj`, `cand_g_form`) consumed by the
#' Stage 2 APS routine. The "simple" variant — used by the paper — considers
#' only main-terms GLMs adjusting for one covariate at a time, plus the
#' unadjusted estimator (`'U'`). The "fancy" variant adds candidate
#' estimators that adjust for the full covariate set with a chosen prediction
#' algorithm (e.g., `'step'`, `'step.interaction'`).
#'
#' @param all_cand Character vector of candidate covariate names.
#' @param cand_q_form_fancy Optional character vector of Qbar algorithms that
#'   should be evaluated against the full covariate set. Default `NULL`
#'   (simple APS).
#' @param cand_g_form_fancy Optional character vector of propensity-score
#'   algorithms for the full covariate set.
#'
#' @return A list with four elements:
#' * `cand_q_adj` — list of candidate adjustment-variable sets for Qbar(A,W)
#' * `cand_q_form` — character vector of corresponding algorithms
#' * `cand_g_adj` — list of candidate adjustment-variable sets for g(A|W)
#' * `cand_g_form` — character vector of corresponding algorithms
#'
#' The lengths of `cand_q_adj` and `cand_q_form` match; same for the g side.
#' Each entry of `cand_q_adj` is itself a character vector of variable names,
#' or the special value `'U'` which means "unadjusted (intercept only)".
#'
#' @export
cand_adj <- function(all_cand,
                     cand_q_form_fancy = NULL,
                     cand_g_form_fancy = NULL) {

  all_cand_u_list <- as.list(c("U", all_cand))

  # Default: main-terms GLM with each candidate alone.
  cand_q_form <- rep("glm", length(all_cand_u_list))
  cand_g_form <- rep("glm", length(all_cand_u_list))

  if (is.null(cand_q_form_fancy)) {
    cand_q_adj <- all_cand_u_list
  } else {
    cand_q_adj <- c(all_cand_u_list, rep(list(all_cand), length(cand_q_form_fancy)))
    cand_q_form <- c(cand_q_form, cand_q_form_fancy)
  }

  if (is.null(cand_g_form_fancy)) {
    cand_g_adj <- all_cand_u_list
  } else {
    cand_g_adj <- c(all_cand_u_list, rep(list(all_cand), length(cand_g_form_fancy)))
    cand_g_form <- c(cand_g_form, cand_g_form_fancy)
  }

  list(
    cand_q_adj = cand_q_adj, cand_q_form = cand_q_form,
    cand_g_adj = cand_g_adj, cand_g_form = cand_g_form
  )
}
