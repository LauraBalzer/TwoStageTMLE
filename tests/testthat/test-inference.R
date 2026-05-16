# Tests for get_inference: CI / p-value / coverage / one-sided variant.

test_that("Two-sided RD with df > 40 uses normal", {
  # psi_hat = 0.05, se = 0.025 -> CI roughly [0.001, 0.099], excludes 0.
  inf <- TwoStageTMLE:::get_inference(
    goal = "RD", psi = 0, psi_hat = 0.05, se = 0.025, df = 100
  )
  z <- stats::qnorm(0.975)
  expect_equal(inf$ci_lo, 0.05 - z * 0.025)
  expect_equal(inf$ci_hi, 0.05 + z * 0.025)
  expect_equal(inf$pval, 2 * stats::pnorm(2, lower.tail = FALSE), tolerance = 1e-8)
  expect_false(inf$cover)
})

test_that("Two-sided RD with df <= 40 uses t", {
  inf <- TwoStageTMLE:::get_inference(
    goal = "RD", psi = 0, psi_hat = 0.05, se = 0.025, df = 28
  )
  t_crit <- stats::qt(0.975, df = 28)
  expect_equal(inf$ci_lo, 0.05 - t_crit * 0.025)
  expect_equal(inf$ci_hi, 0.05 + t_crit * 0.025)
})

test_that("aRR / OR are back-transformed via exp()", {
  inf_rr <- TwoStageTMLE:::get_inference(
    goal = "aRR", psi = 1, psi_hat = log(1.2), se = 0.1, df = 28
  )
  expect_equal(inf_rr$est, 1.2)
  t_crit <- stats::qt(0.975, df = 28)
  expect_equal(inf_rr$ci_lo, exp(log(1.2) - t_crit * 0.1))
  expect_equal(inf_rr$ci_hi, exp(log(1.2) + t_crit * 0.1))
})

test_that("One-sided test requires alt_smaller", {
  expect_error(
    TwoStageTMLE:::get_inference(
      goal = "RD", psi = 0, psi_hat = -0.1, se = 0.04, df = 28,
      one_sided = TRUE
    ),
    "alt_smaller"
  )
})

test_that("One-sided alt_smaller = TRUE matches pt() lower-tail", {
  inf <- TwoStageTMLE:::get_inference(
    goal = "RD", psi = 0, psi_hat = -0.1, se = 0.04, df = 28,
    one_sided = TRUE, alt_smaller = TRUE
  )
  expect_equal(inf$pval, stats::pt(-2.5, df = 28, lower.tail = TRUE),
               tolerance = 1e-8)
})

test_that("Coverage and rejection indicators are correct", {
  inf <- TwoStageTMLE:::get_inference(
    goal = "RD", psi = 0, psi_hat = 0.0001, se = 0.05, df = 28
  )
  expect_true(inf$cover)
  expect_equal(inf$reject, 0)

  inf_reject <- TwoStageTMLE:::get_inference(
    goal = "RD", psi = 0, psi_hat = 0.5, se = 0.05, df = 28
  )
  expect_false(inf_reject$cover)
  expect_equal(inf_reject$reject, 1)
})
