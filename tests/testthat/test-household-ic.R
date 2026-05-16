# Tests for the household_ic reimplementation.
# Verifies behaviour matches Josh's transcription of ltmle:::HouseholdIC.

test_that("household_ic: NULL id is a pass-through", {
  m <- matrix(rnorm(20), ncol = 2)
  expect_identical(TwoStageTMLE:::household_ic(m, id = NULL), m)
})

test_that("household_ic: sums within cluster and scales by n_clusters/n_rows", {
  # 4 records in 2 clusters; sums within cluster, then * 2/4 = 1/2.
  ic <- matrix(c(1, 2, 3, 4), ncol = 1)
  ids <- c(1, 1, 2, 2)
  out <- TwoStageTMLE:::household_ic(ic, ids)
  expect_equal(unname(as.vector(out)), c(3 * 0.5, 7 * 0.5))
})

test_that("household_ic: errors when id length mismatches", {
  ic <- matrix(1:6, ncol = 2)
  expect_error(
    TwoStageTMLE:::household_ic(ic, id = c(1, 2)),
    "must equal the number of rows"
  )
})

test_that("household_ic: coerces a vector to a 1-column matrix", {
  ic <- c(1, 2, 3, 4)
  ids <- c(1, 1, 2, 2)
  out <- TwoStageTMLE:::household_ic(ic, ids)
  expect_true(is.matrix(out))
  expect_equal(ncol(out), 1)
  expect_equal(unname(as.vector(out)), c(3 * 0.5, 7 * 0.5))
})
