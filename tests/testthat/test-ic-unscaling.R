# Focused tests for the IC-unscaling bug fix.
#
# The bug in from_analysis/stage2.R:249-250 added `scale_value_min` to
# the IC components. Because var(X + c) = var(X), this doesn't affect
# any user-facing variance / SE / CI / p-value — but it does affect:
#   (a) returned IC vectors (DY1, DY0, DY) when accessed directly, and
#   (b) the verbose "Solve EIF: mean(DY)" diagnostic.
#
# The package applies the fix: residuals are multiplied by (max - min)
# only, with no additive shift.

# Helper: minimal synthetic Stage 2 input with fully controllable Y.
synth_clust <- function(Y, J = 6) {
  stopifnot(length(Y) == J, J %% 2 == 0)
  data.table::data.table(
    id = seq_len(J),
    pair = rep(seq_len(J / 2), each = 2),
    nIndv = rep(50L, J),
    nIndv_Y = rep(50L, J),
    U = 1L,
    A = rep(c(0, 1), J / 2),
    X1.c = stats::runif(J, -1, 1),
    X2.c = stats::runif(J, -1, 1),
    Y = Y,
    alpha = 1
  )
}

test_that("Identity rescaling: outcome already in [0,1] is untouched", {
  rescale <- TwoStageTMLE:::scale_outcome(c(0.1, 0.4, 0.9))
  expect_equal(rescale$scale_min, 0)
  expect_equal(rescale$scale_max, 1)
  expect_equal(rescale$y, c(0.1, 0.4, 0.9))
})

test_that("Scaling shifts when min(Y) < 0 or max(Y) > 1", {
  rescale_neg <- TwoStageTMLE:::scale_outcome(c(-0.5, 0.2, 0.8))
  expect_equal(rescale_neg$scale_min, -0.5)
  expect_equal(rescale_neg$scale_max, 1)
  expect_equal(rescale_neg$y, c(0, 0.7 / 1.5, 1.3 / 1.5))

  rescale_pos <- TwoStageTMLE:::scale_outcome(c(0, 1.2, 1.8))
  expect_equal(rescale_pos$scale_min, 0)
  expect_equal(rescale_pos$scale_max, 1.8)
})

test_that("IC unscaling is multiplication only (residual semantics)", {
  # Construct vdata with controlled residual contributions.
  # If DY1_scaled = c(-0.1, 0.05, ...), then DY1_unscaled = DY1_scaled * (max-min).
  vdata <- data.table::data.table(
    id = 1:6,
    pair = rep(1:3, each = 2),
    alpha = 1,
    H.1W = c(2, 0, 2, 0, 2, 0),
    H.0W = c(0, 2, 0, 2, 0, 2),
    Y = c(0.6, 0.4, 0.7, 0.3, 0.5, 0.5),
    Qbar1W.star = c(0.5, 0.5, 0.6, 0.4, 0.5, 0.5),
    Qbar0W.star = c(0.5, 0.5, 0.6, 0.4, 0.5, 0.5)
  )

  out1 <- TwoStageTMLE:::get_ic_variance(
    goal = "RD", target = "clust", vdata = vdata,
    R1 = 0.55, R0 = 0.45,
    scale_min = 0, scale_max = 1, doing_cv = FALSE
  )
  out2 <- TwoStageTMLE:::get_ic_variance(
    goal = "RD", target = "clust", vdata = vdata,
    R1 = 0.55, R0 = 0.45,
    scale_min = -0.5, scale_max = 1.5, doing_cv = FALSE
  )

  # The IC components for run 2 should be the run-1 ICs multiplied by
  # (1.5 - (-0.5)) / (1 - 0) = 2 — with NO additive shift.
  expect_equal(out2$DY1, out1$DY1 * 2)
  expect_equal(out2$DY0, out1$DY0 * 2)
  expect_equal(out2$DY,  out1$DY  * 2)

  # Variances scale by the square of the multiplier (no shift involved).
  expect_equal(out2$var_R1, out1$var_R1 * 4)
  expect_equal(out2$var_R0, out1$var_R0 * 4)
  expect_equal(out2$var_break, out1$var_break * 4)
})
