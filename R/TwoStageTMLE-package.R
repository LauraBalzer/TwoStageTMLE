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
#' The package exposes general, composable building blocks rather than
#' paper-specific drivers. Specific published methods (Balzer et al. 2021;
#' Nakato et al. 2026) are documented in the vignettes as compositions of these
#' blocks (and kept as internal reference implementations for verification).
#'
#' @section Main functions:
#' Stage 1 (cluster-specific endpoints):
#' * [stage1_mean()] — one measurement-adjusted cluster mean.
#' * [stage1_ratio()] — ratio of two cluster means.
#' * [stage1_apply()] — general engine: map any per-cluster estimator over clusters.
#'
#' Stage 2 and helpers:
#' * [stage2_tmle()] — cluster-level TMLE with Adaptive Pre-specification.
#' * [add_weights()] — analysis weights; [cand_adj()] — APS candidate sets.
#' * [as_stage2_list()], [flatten_stage2()] — list-based result helpers.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom data.table data.table as.data.table copy is.data.table setDT
#'   setnames := .SD .N .I
#' @importFrom stats glm predict plogis qlogis qt qnorm pt pnorm var model.matrix
#'   step binomial gaussian
NULL
