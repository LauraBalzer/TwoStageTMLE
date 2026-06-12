# Tests for the layered Stage-1 API: stage1_apply / stage1_mean / stage1_ratio
# (and the internal fit_censored_mean).
#
# Strategy mirrors test-stage1-mediated.R: generate a small CRT with a censored
# outcome, then check the new functions against an inline reference that calls
# ltmle() with identical arguments -> they must agree to machine precision.

gen_censored_crt <- function(J = 8, N_mean = 90, seed = 1L) {
  set.seed(seed)
  rows <- lapply(seq_len(J), function(j) {
    N  <- round(stats::rnorm(1, N_mean, 5))
    W1 <- stats::rnorm(N)
    W2 <- stats::rbinom(N, 1, 0.5)
    A  <- as.integer(j %% 2 == 0)
    # outcome (eligibility) measured only among ~75% -> NA = censored
    elig_full <- stats::rbinom(N, 1, stats::plogis(-0.2 + 0.7 * W1 + 0.3 * W2))
    measured  <- stats::runif(N) < 0.78
    elig <- elig_full; elig[!measured] <- NA_integer_
    # a fully-observed binary "started" outcome for the ratio numerator
    started <- stats::rbinom(N, 1, 0.35 + 0.05 * A)
    data.table::data.table(id = j, A = A, W1 = W1, W2 = W2,
                           eligible = elig, started = started, country = j %% 2)
  })
  data.table::rbindlist(rows)
}

# inline reference: unadjusted + adjusted censored mean for one column
ref_censored_pair <- function(oc, col, adj, lib) {
  y <- oc[[col]]
  base <- data.frame(C = ltmle::BinaryToCensoring(is.censored = is.na(y)), Y = y)
  f0 <- suppressMessages(suppressWarnings(
    ltmle::ltmle(base, Anodes = NULL, abar = NULL, Cnodes = "C", Ynodes = "Y",
                 SL.library = lib, estimate.time = FALSE)))
  adjf <- cbind(as.data.frame(oc[, adj, with = FALSE]), base)
  f1 <- suppressMessages(suppressWarnings(
    ltmle::ltmle(adjf, Anodes = NULL, abar = NULL, Cnodes = "C", Ynodes = "Y",
                 SL.library = lib, estimate.time = FALSE)))
  list(unadj = unname(f0$estimates["tmle"]), tmle = unname(f1$estimates["tmle"]))
}

test_that("stage1_apply assembles the canonical cluster table", {
  skip_if_not_installed("ltmle")
  dat <- gen_censored_crt(J = 6, N_mean = 60, seed = 11L)
  s1 <- stage1_apply(dat, function(oc) c(m = mean(oc$started)),
                     clust_cov = "country")
  expect_s3_class(s1, "data.table")
  expect_equal(nrow(s1), 6L)
  expect_true(all(c("id", "nIndv", "U", "A", "country", "m") %in% names(s1)))
  expect_true(all(s1$U == 1L))
})

test_that("stage1_mean matches an inline ltmle reference (machine precision)", {
  skip_if_not_installed("ltmle")
  dat <- gen_censored_crt(J = 6, N_mean = 80, seed = 21L)
  lib <- "glm"; adj <- c("W1", "W2")

  s1_adj   <- stage1_mean(dat, y_col = "eligible", adj_covars = adj,  sl_library = lib)[order(id)]
  s1_unadj <- stage1_mean(dat, y_col = "eligible", adj_covars = NULL, sl_library = lib)[order(id)]

  ref <- t(vapply(sort(unique(dat$id)), function(j) {
    r <- ref_censored_pair(dat[id == j], "eligible", adj, lib)
    c(unadj = r$unadj, tmle = r$tmle)
  }, numeric(2)))

  expect_equal(s1_adj$est,   ref[, "tmle"],  tolerance = 1e-12)   # adjusted
  expect_equal(s1_unadj$est, ref[, "unadj"], tolerance = 1e-12)   # unadjusted = adj_covars NULL
})

test_that("stage1_ratio (plain-mean numerator / censored denominator) matches reference", {
  skip_if_not_installed("ltmle")
  dat <- gen_censored_crt(J = 6, N_mean = 80, seed = 31L)
  lib <- "glm"; w1 <- c("W1", "W2")

  # opal coverage/yield shape: numerator `started` is fully observed (-> plain
  # mean via the no-NA short-circuit), denominator `eligible` is censored.
  s1_adj   <- stage1_ratio(dat, num_col = "started", den_col = "eligible",
                           den_adj = w1, sl_library = lib)[order(id)]
  s1_unadj <- stage1_ratio(dat, num_col = "started", den_col = "eligible",
                           den_adj = NULL, sl_library = lib)[order(id)]

  ref <- t(vapply(sort(unique(dat$id)), function(j) {
    oc  <- dat[id == j]
    num <- mean(oc$started)
    den <- ref_censored_pair(oc, "eligible", w1, lib)
    c(unadj = num / den$unadj, tmle = num / den$tmle)
  }, numeric(2)))

  expect_equal(s1_adj$est,   ref[, "tmle"],  tolerance = 1e-12)   # adjusted denominator
  expect_equal(s1_unadj$est, ref[, "unadj"], tolerance = 1e-12)   # unadjusted denominator
})

test_that("stage1_mean accepts a function-valued adj_covars (per-cluster set)", {
  skip_if_not_installed("ltmle")
  dat <- gen_censored_crt(J = 6, N_mean = 80, seed = 41L)
  lib <- "glm"
  # shrink the adjustment set for clusters with < 78 rows
  adj_fun <- function(oc) if (nrow(oc) < 78) "W1" else c("W1", "W2")

  s1 <- stage1_mean(dat, y_col = "eligible", adj_covars = adj_fun, sl_library = lib)
  s1 <- s1[order(id)]

  ref <- t(vapply(sort(unique(dat$id)), function(j) {
    oc  <- dat[id == j]
    adj <- if (nrow(oc) < 78) "W1" else c("W1", "W2")
    r   <- ref_censored_pair(oc, "eligible", adj, lib)
    c(unadj = r$unadj, tmle = r$tmle)
  }, numeric(2)))

  expect_equal(s1$est, ref[, "tmle"], tolerance = 1e-12)
  # at least one cluster on each side of the size threshold, so the function path matters
  sizes <- dat[, .N, by = id]$N
  expect_true(any(sizes < 78) && any(sizes >= 78))
})
