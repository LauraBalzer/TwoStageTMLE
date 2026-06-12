# Tests for the list-based result helpers: as_stage2_list() + flatten_stage2().

make_clust <- function(J = 30, seed = 5L) {
  set.seed(seed)
  d <- data.frame(
    id = 1:J, A = rep(c(1L, 0L), each = J / 2), U = 1,
    W = rbinom(J, 1, 0.5), nIndv = sample(40:100, J, replace = TRUE)
  )
  d$Y <- plogis(-0.4 + 0.5 * d$A + 0.5 * d$W + rnorm(J, 0, 0.4)) * 0.6 + 0.1
  d$alpha <- 1
  d
}

test_that("as_stage2_list builds the leaf shape", {
  skip_if_not_installed("ltmle")
  d <- make_clust()
  rd  <- stage2_tmle(d, goal = "RD",  target = "clust")
  arr <- stage2_tmle(d, goal = "aRR", target = "clust")

  leaf <- as_stage2_list(RD = rd, aRR = arr)

  expect_true(all(c("arms", "RD", "aRR") %in% names(leaf)))
  expect_named(leaf$arms, c("txt", "con"))
  expect_named(leaf$arms$txt, c("est", "ci_lo", "ci_hi", "se"))
  # q_adj is carried PER GOAL (APS is run per goal)
  expect_named(leaf$RD, c("est", "ci_lo", "ci_hi", "se", "pval", "q_adj"))
  expect_named(leaf$aRR, c("est", "ci_lo", "ci_hi", "se", "pval", "q_adj"))
  # leaf values agree with the source data.tables
  expect_equal(unname(leaf$RD["est"]),   rd$est)
  expect_equal(unname(leaf$aRR["pval"]), arr$pval)
  expect_equal(unname(leaf$RD["q_adj"]),  as.numeric(rd$q_adj))
  expect_equal(unname(leaf$aRR["q_adj"]), as.numeric(arr$q_adj))
  expect_equal(unname(leaf$arms$txt["est"]), rd$txt_est)
})

test_that("flatten_stage2 collapses a nested list to one row per leaf", {
  skip_if_not_installed("ltmle")
  d <- make_clust()
  rd  <- stage2_tmle(d, goal = "RD",  target = "clust")
  arr <- stage2_tmle(d, goal = "aRR", target = "clust")
  leaf <- as_stage2_list(RD = rd, aRR = arr)

  nested <- list(ALL = list(tmle = leaf, unadj = leaf))
  flat <- flatten_stage2(nested, label_names = c("subgroup", "estimator"))

  expect_equal(nrow(flat), 2L)
  expect_true(all(c("subgroup", "estimator", "txt_est",
                    "RD_est", "RD_pval", "RD_q_adj", "aRR_est", "aRR_q_adj") %in% names(flat)))
  expect_equal(flat$subgroup, c("ALL", "ALL"))
  expect_equal(flat$estimator, c("tmle", "unadj"))
  expect_equal(flat$RD_est[1], rd$est)
  expect_equal(flat$aRR_q_adj[1], as.numeric(arr$q_adj))
})
