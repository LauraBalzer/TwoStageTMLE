# Tests for stage1_mediated() — the Nakato et al. (2026) ratio-based
# Stage 1 estimator for the multilevel-mediation-missing data problem.
#
# Strategy:
#  (1) Generate a small CRT with a DGP that mirrors
#      Nakato/multi-level-mediation-missingness-main/generate_data.R.
#  (2) Run stage1_mediated().
#  (3) Run an inline reference reimplementation of Nakato's getYc() on
#      the same data and check the `tmle` and `unadj` columns match
#      to machine precision (both implementations call the same
#      ltmle() with identical arguments, so they must agree).
#  (4) Smoke-test the end-to-end driver two_stage_tmle_mediated().

# -- inline DGP (Nakato 2026 'main' DGP, condensed) ---------------------

gen_mediated_cluster <- function(N, j) {
  UE1 <- rep(stats::runif(1, -1, 1), N)
  UE2 <- rep(stats::runif(1, -1, 1), N)
  E1 <- rep(stats::rnorm(1, UE1[1], 0.5), N)
  E2 <- rep(stats::rnorm(1, UE2[1], 0.5), N)
  W1 <- stats::rnorm(N, 0, 0.5)
  W2 <- stats::rnorm(N, 0, 0.5)
  W3 <- stats::rbinom(N, 1, 0.5)

  pY1_star <- stats::plogis(0.5 + W1 + W2 - W3 + 0.25 * E1 + 0.25 * E2)
  Y1_star  <- as.numeric(stats::runif(N) < pY1_star)

  data.table::data.table(
    id = j, n = N, E1 = E1, E2 = E2,
    W1 = W1, W2 = W2, W3 = W3,
    UE1 = UE1, UE2 = UE2, Y1_star = Y1_star
  )
}

gen_mediated_crt <- function(J = 12, N_mean = 120, N_sd = 5, seed = 1L) {
  set.seed(seed)
  N_j <- round(stats::rnorm(J, N_mean, N_sd))
  dat <- data.table::rbindlist(
    lapply(seq_len(J), function(j) gen_mediated_cluster(N_j[j], j))
  )
  A_per_cluster <- sample(c(rep(1L, J %/% 2), rep(0L, J - J %/% 2)))
  dat[, A := A_per_cluster[id]]

  # Screening probability depends on A (the intervention boosts screening)
  pM <- stats::plogis(-0.5 + 0.6 * dat$A + 0.5 * dat$W1 + 0.4 * dat$W2 -
                        0.4 * dat$A * dat$W3 + 0.4 * (1 - dat$A) * dat$W3 +
                        0.1 * dat$UE1 + 0.1 * dat$UE2)
  dat[, M := as.integer(stats::runif(.N) < pM)]
  dat[, Y1 := Y1_star * M]

  # Outcome: only happens among the screened-and-eligible
  pY2 <- stats::plogis(0.2 + 0.1 * dat$A + dat$W1 + 0.5 * dat$W2 +
                        2 * dat$W3 + 0.2 * dat$E1 + 0.2 * dat$E2)
  Y2  <- as.integer(stats::runif(nrow(dat)) < pY2)
  dat[, Y2 := ifelse(Y1_star == 1 & M == 1, Y2, 0L)]

  dat[, c("UE1", "UE2", "Y1_star", "n") := NULL]
  dat[]
}

# -- inline reference reimplementation of Nakato's getYc() --------------

ref_getYc <- function(O, indv_cov, clust_cov,
                     sl_library = c("SL.mean", "SL.glm", "SL.gam")) {
  O <- as.data.frame(O)
  clusters <- sort(unique(O$id))
  rows <- lapply(clusters, function(j) {
    Oc <- O[O$id == j, , drop = FALSE]
    numerator <- mean(Oc$Y2)
    den_unadj <- mean(Oc[Oc$M == 1, "Y1"])
    den_fit <- suppressMessages(suppressWarnings(
      ltmle::ltmle(
        data = Oc[, c(indv_cov, "M", "Y1")],
        Anodes = "M", Ynodes = "Y1", abar = 1,
        SL.library = sl_library, estimate.time = FALSE, stratify = TRUE
      )
    ))
    den_tmle <- unname(den_fit$estimates["tmle"])
    data.frame(
      id      = j,
      unadj   = numerator / den_unadj,
      tmle    = numerator / den_tmle
    )
  })
  do.call(rbind, rows)
}

# -- tests --------------------------------------------------------------

test_that("stage1_mediated() returns the expected shape", {
  skip_if_not_installed("ltmle")
  skip_if_not_installed("gam")

  dat <- gen_mediated_crt(J = 8, N_mean = 80, seed = 42L)
  s1 <- stage1_mediated(
    data = dat,
    ind_cov  = c("W1", "W2", "W3"),
    clust_cov = c("E1", "E2"),
    sl_library = c("SL.mean", "SL.glm")
  )

  expect_s3_class(s1, "data.table")
  expect_equal(nrow(s1), 8L)
  expect_true(all(c("id", "A", "E1", "E2", "nIndv", "U",
                    "screened", "eligible", "unadj", "tmle") %in% names(s1)))
  expect_true(all(s1$U == 1L))
  expect_true(all(is.finite(s1$tmle)))
  # ratio of probabilities is bounded by 1 in expectation; finite-sample
  # excursions slightly above 1 are allowed but not far above.
  expect_true(all(s1$tmle >= 0 & s1$tmle <= 1.5))
})

test_that("stage1_mediated() matches Nakato's getYc() to machine precision", {
  skip_if_not_installed("ltmle")
  skip_if_not_installed("gam")

  dat <- gen_mediated_crt(J = 6, N_mean = 80, seed = 99L)
  ind_cov  <- c("W1", "W2", "W3")
  clust_cov <- c("E1", "E2")
  lib <- c("SL.mean", "SL.glm")

  # SuperLearner consumes RNG state for CV-fold assignment, so reset the
  # seed identically before each call to make the parity comparison exact.
  set.seed(20260101L)
  s1  <- stage1_mediated(data = dat, ind_cov = ind_cov,
                         clust_cov = clust_cov, sl_library = lib)
  set.seed(20260101L)
  ref <- ref_getYc(dat, indv_cov = ind_cov, clust_cov = clust_cov,
                   sl_library = lib)

  # Compare cluster-by-cluster after aligning on id
  s1 <- s1[order(id)]
  ref <- ref[order(ref$id), ]
  expect_equal(s1$id, ref$id)
  expect_equal(s1$unadj, ref$unadj, tolerance = 1e-12)
  expect_equal(s1$tmle,  ref$tmle,  tolerance = 1e-12)
})

test_that("stage1_mediated() input validation catches missing columns", {
  bad <- data.table::data.table(id = 1:5, A = 0L, M = 1L, Y1 = 0L)  # no Y2
  expect_error(
    stage1_mediated(bad, ind_cov = character(0), clust_cov = character(0)),
    "missing required columns.*Y2"
  )
})

test_that("stage1_mediated() input validation catches non-binary columns", {
  dat <- gen_mediated_crt(J = 4, N_mean = 30, seed = 7L)
  dat$M[1L] <- 2L
  expect_error(
    stage1_mediated(dat, ind_cov = c("W1"), clust_cov = c("E1")),
    "M.*binary"
  )
})

test_that("two_stage_tmle_mediated() runs end-to-end", {
  skip_if_not_installed("ltmle")
  skip_if_not_installed("gam")

  dat <- gen_mediated_crt(J = 12, N_mean = 100, seed = 2024L)
  cand <- cand_adj(c("E1", "E2"))

  out <- suppressWarnings(two_stage_tmle_mediated(
    data = dat,
    ind_cov  = c("W1", "W2", "W3"),
    clust_cov = c("E1", "E2"),
    sl_library = c("SL.mean", "SL.glm"),
    goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    remove_pscore = TRUE
  ))

  expect_named(out, c("est", "stage1"))
  expect_s3_class(out$est, "data.table")
  expect_equal(nrow(out$est), 1L)
  expect_true(all(c("est", "ci_lo", "ci_hi", "se", "pval") %in% names(out$est)))
  expect_true(is.finite(out$est$est))
  expect_true(is.finite(out$est$se) && out$est$se > 0)
  expect_true(out$est$ci_lo <= out$est$est && out$est$est <= out$est$ci_hi)

  # Stage 1 cluster table is intact and has all four candidate endpoints.
  expect_equal(nrow(out$stage1), 12L)
  expect_true(all(c("screened", "eligible", "unadj", "tmle") %in% names(out$stage1)))
})
