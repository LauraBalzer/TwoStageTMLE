# Shared test helpers.

# Path to the bundled reference fixtures.
fixture_path <- function(name) {
  testthat::test_path("fixtures", name)
}

# Cluster-level input ready for stage2_tmle(), from the Tier A complex_eff
# deep dive (locked output from TwoStageTMLE-main code).
get_cluster_input <- function() {
  deep <- readRDS(fixture_path("tierA_deep_complex_eff_rep1.rds"))
  d <- data.table::as.data.table(deep$data.C)
  d[, alpha := 1]
  d[, nIndv_Y := nIndv]
  data.table::setnames(d, "tmle", "Y")
  list(data = d, deep = deep)
}

# Old-fixture column-name normaliser so parity comparisons are symmetric.
extract_old <- function(x) {
  if (is.list(x) && "est.df" %in% names(x)) x <- x$est.df
  est_val <- if ("psi.hat" %in% names(x)) x$psi.hat else x$est
  list(
    est = est_val, se = x$se,
    ci_lo = x$CI.lo, ci_hi = x$CI.hi,
    pval = x$pval,
    cover = if (is.null(x$cover)) NA else as.logical(x$cover),
    reject = if (is.null(x$reject)) NA else as.logical(x$reject)
  )
}

extract_new <- function(x) {
  if (is.list(x) && "est" %in% names(x) && data.table::is.data.table(x$est)) {
    x <- x$est
  }
  list(
    est = x$est, se = x$se,
    ci_lo = x$ci_lo, ci_hi = x$ci_hi,
    pval = x$pval,
    cover = as.logical(x$cover),
    reject = as.logical(x$reject)
  )
}

expect_parity <- function(old, new, tol = 1e-8) {
  o <- extract_old(old); n <- extract_new(new)
  testthat::expect_equal(n$est,    o$est,    tolerance = tol)
  testthat::expect_equal(n$se,     o$se,     tolerance = tol)
  testthat::expect_equal(n$ci_lo,  o$ci_lo,  tolerance = tol)
  testthat::expect_equal(n$ci_hi,  o$ci_hi,  tolerance = tol)
  testthat::expect_equal(n$pval,   o$pval,   tolerance = tol)
  if (!is.na(o$cover))  testthat::expect_equal(n$cover,  o$cover)
  if (!is.na(o$reject)) testthat::expect_equal(n$reject, o$reject)
}

# Wrapper that suppresses the "non-integer #successes in a binomial glm"
# warnings that are expected when fitting a binomial GLM to a probability-
# valued outcome (Stage 1 returns cluster-level means in [0,1]).
quiet_stage2 <- function(...) suppressWarnings(stage2_tmle(...))
