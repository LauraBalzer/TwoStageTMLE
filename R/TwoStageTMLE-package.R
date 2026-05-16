#' TwoStageTMLE: Two-Stage TMLE for Cluster Randomized Trials
#'
#' Implements the Two-Stage Targeted Minimum Loss-Based Estimator of
#' Balzer et al. (2021) for cluster randomized trials. See
#' <https://arxiv.org/abs/2106.15737> for the methodology paper.
#' Stage 1 uses an individual-level TMLE within each cluster separately to
#' adjust for differential outcome measurement. Stage 2 uses a
#' cluster-level TMLE with Adaptive Pre-specification to select adjustment
#' variables that optimise precision while preserving Type-I error control.
#'
#' @section Main functions:
#' * [two_stage_tmle()] — end-to-end driver: data in, effect estimate out.
#' * [stage1_tmle()] — Stage 1 only, returns cluster-specific endpoint estimates.
#' * [stage2_tmle()] — Stage 2 only, takes cluster-specific endpoints and
#'   returns the intervention-effect estimate with inference.
#' * [cand_adj()] — build candidate adjustment specifications for APS.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom data.table data.table as.data.table copy is.data.table setDT
#'   setnames := .SD .N .I
#' @importFrom stats glm predict plogis qlogis qt qnorm pt pnorm var model.matrix
#'   step binomial gaussian
NULL
