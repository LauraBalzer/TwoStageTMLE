#' Estimate a measurement-adjusted mean from censored data (internal)
#'
#' Thin wrapper around [ltmle::ltmle()] that estimates \eqn{E[Y]} from a data
#' frame whose columns are `[optional baseline covariates W..., C, Y]`, treating
#' `C` (a censoring node built with [ltmle::BinaryToCensoring()]) as the
#' missingness mechanism for `Y`. This is the atomic Stage-1 building block used
#' by [stage1_mean()] and [stage1_ratio()].
#'
#' Defaults match `ltmle`'s own defaults so that callers which omit an argument
#' get identical numbers to a bare `ltmle()` call; callers pass `sl_library`,
#' `det_q_fun`, and/or `variance_method` explicitly when they need non-default
#' behaviour (e.g. `variance_method = "ic"` and a deterministic-Q closure for
#' the HIV viral-suppression endpoint).
#'
#' @param data A `data.frame` whose last two columns are `C` (censoring node)
#'   and `Y` (outcome); any columns before `C` are treated as baseline
#'   covariates by `ltmle`.
#' @param sl_library SuperLearner library passed to [ltmle::ltmle()]. Default
#'   `"glm"` (plain main-terms GLM; no SuperLearner cross-validation).
#' @param det_q_fun Optional `deterministic.Q.function` passed to `ltmle`.
#' @param variance_method One of `"tmle"` (default), `"ic"`, `"iptw"`.
#'
#' @return The scalar TMLE point estimate (`fit$estimates["tmle"]`, which equals
#'   `summary(fit)$treatment$estimate`).
#'
#' @keywords internal
fit_censored_mean <- function(data, sl_library = "glm", det_q_fun = NULL,
                              variance_method = "tmle") {
  fit <- suppressMessages(suppressWarnings(
    ltmle::ltmle(
      data                     = data,
      Anodes                   = NULL,
      abar                     = NULL,
      Cnodes                   = "C",
      Ynodes                   = "Y",
      SL.library               = sl_library,
      deterministic.Q.function = det_q_fun,
      variance.method          = variance_method,
      estimate.time            = FALSE
    )
  ))
  unname(fit$estimates["tmle"])
}
