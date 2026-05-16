# Tests for input validation messages.

test_that("Stage 2 errors on missing required columns", {
  setup <- get_cluster_input()
  d <- setup$data
  d_no_alpha <- data.table::copy(d); d_no_alpha[, alpha := NULL]
  expect_error(stage2_tmle(data = d_no_alpha, goal = "RD"),
               "missing required columns.*alpha")
})

test_that("Stage 2 errors when U is not all 1", {
  setup <- get_cluster_input()
  d <- setup$data
  d[1, U := 2]
  expect_error(stage2_tmle(data = d, goal = "RD"),
               "U.*must equal 1")
})

test_that("Stage 2 errors on non-binary A", {
  setup <- get_cluster_input()
  d <- setup$data
  d[1, A := 2]
  expect_error(stage2_tmle(data = d, goal = "RD"),
               "A.*binary")
})

test_that("Stage 2 errors on break_match = FALSE without pair column", {
  setup <- get_cluster_input()
  d <- setup$data
  d[, pair := NULL]
  expect_error(stage2_tmle(data = d, goal = "RD", break_match = FALSE),
               "pair.*column")
})

test_that("Stage 1 errors on missing covariate column", {
  expect_error(stage1_tmle(
    data = data.table::data.table(id = 1:5, A = 0, Delta = 1, Y = 0),
    ind_cov = c("nonexistent"), clust_cov = character(0)
  ), "missing required columns")
})
