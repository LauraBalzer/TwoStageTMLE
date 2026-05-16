#' Wald-type confidence intervals and hypothesis test
#'
#' Returns point estimate, two-sided 100*(1 - sig_level)% confidence interval,
#' p-value, and coverage / rejection indicators given a point estimate on a
#' (possibly transformed) scale together with its standard error. For
#' `goal = "aRR"` and `goal = "OR"` the input `psi_hat` and `se` are assumed
#' to be on the log scale; the returned estimate and confidence interval are
#' back-transformed to the original scale via `exp()`.
#'
#' Uses Student's t-distribution when `df <= 40`, normal otherwise — matching
#' the heuristic in `from_analysis/stage2.R`.
#'
#' @param goal One of `"RD"`, `"aRR"`, `"OR"`.
#' @param psi True parameter value if known (e.g., in simulation). Used to
#'   compute bias and coverage. Default `NA`.
#' @param psi_hat Point estimate. For `aRR`/`OR` this must be on the log scale.
#' @param se Standard error of `psi_hat` on the (possibly transformed) scale.
#' @param df Degrees of freedom for the t-distribution. If `df > 40` the normal
#'   approximation is used.
#' @param sig_level Significance level. Default `0.05`.
#' @param one_sided Logical. If `TRUE`, a one-sided p-value is reported.
#' @param alt_smaller Logical. If `one_sided = TRUE`, `TRUE` means the
#'   alternative is that the intervention *decreases* the mean outcome.
#'
#' @return A one-row `data.table` with columns `est`, `ci_lo`, `ci_hi`, `se`,
#'   `pval`, `bias`, `cover`, `reject`.
#'
#' @keywords internal
get_inference <- function(goal = "RD", psi = NA_real_, psi_hat, se,
                          df = 99, sig_level = 0.05,
                          one_sided = FALSE, alt_smaller = NULL) {
  goal <- match.arg(goal, c("RD", "aRR", "OR"))
  if (one_sided && is.null(alt_smaller)) {
    stop("`one_sided = TRUE` requires `alt_smaller` to be TRUE or FALSE.")
  }

  tstat <- psi_hat / se

  if (df > 40) {
    cutoff <- stats::qnorm(sig_level / 2, lower.tail = FALSE)
    pval <- if (one_sided) {
      stats::pnorm(tstat, lower.tail = alt_smaller)
    } else {
      2 * stats::pnorm(abs(tstat), lower.tail = FALSE)
    }
  } else {
    cutoff <- stats::qt(sig_level / 2, df = df, lower.tail = FALSE)
    pval <- if (one_sided) {
      stats::pt(tstat, df = df, lower.tail = alt_smaller)
    } else {
      2 * stats::pt(abs(tstat), df = df, lower.tail = FALSE)
    }
  }

  ci_lo <- psi_hat - cutoff * se
  ci_hi <- psi_hat + cutoff * se

  if (goal != "RD") {
    psi_hat <- exp(psi_hat)
    ci_lo <- exp(ci_lo)
    ci_hi <- exp(ci_hi)
  }

  bias <- psi_hat - psi
  cover <- ci_lo <= psi & psi <= ci_hi
  reject <- as.numeric(pval < sig_level)

  data.table::data.table(
    est = psi_hat, ci_lo = ci_lo, ci_hi = ci_hi, se = se,
    pval = pval, bias = bias, cover = cover, reject = reject
  )
}
