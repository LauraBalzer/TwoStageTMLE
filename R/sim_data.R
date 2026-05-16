# Simulation data-generating process for vignettes and tests.
# Reproduces the DGPs of Balzer et al. (2021) § 5 (main / complex) and
# § 9.6 (supplementary / simple). These functions are internal — exported
# users should not rely on them — but they are deterministic given the
# seed, so they form the basis of the package's reproducibility tests.

#' Pair-match clusters and randomise treatment within pairs
#'
#' Uses [nbpMatching::nonbimatch()] on a univariate matching covariate.
#' Returns the pair id and the cluster-level treatment assignment.
#'
#' @param match_data Numeric vector of length J (one matching covariate per
#'   cluster).
#' @return A `data.frame` with columns `pairs` and `A_c`.
#' @keywords internal
match_me <- function(match_data) {
  if (!requireNamespace("nbpMatching", quietly = TRUE)) {
    stopf("Simulation DGP requires the `nbpMatching` package (Suggests).")
  }
  dist <- nbpMatching::distancematrix(
    nbpMatching::gendistance(data.frame(matchData = match_data))
  )
  matches <- nbpMatching::nonbimatch(dist)
  grpA <- as.numeric(matches$halves[, "Group1.Row"])
  grpB <- as.numeric(matches$halves[, "Group2.Row"])

  J <- length(match_data)
  n_pairs <- J / 2
  pairs <- A_c <- rep(NA_real_, J)
  pairs[grpA] <- seq_len(n_pairs)
  pairs[grpB] <- seq_len(n_pairs)
  A <- stats::rbinom(n_pairs, 1, 0.5)
  A_c[grpA] <- A
  A_c[grpB] <- 1 - A
  data.frame(pairs = pairs, A_c = A_c)
}


#' Aggregate individual-level covariates to cluster level
#' @keywords internal
do_aggregation <- function(O) {
  X_c <- stats::aggregate(O[, c("X1", "X2")], by = list(O$id), FUN = mean)
  n_C <- stats::aggregate(O[, "U"], by = list(O$id), FUN = sum)$x
  X1.c <- rep(X_c[, "X1"], times = n_C)
  X2.c <- rep(X_c[, "X2"], times = n_C)
  nIndv <- rep(n_C, times = n_C)
  cbind(O, nIndv = nIndv, X1.c = X1.c, X2.c = X2.c)
}


#' Generate baseline covariates and cluster-level randomisation for the main DGP
#' @keywords internal
gen_baseline_complex <- function(n, J, N) {
  id <- rep(seq_len(J), times = n)
  U1_cc <- stats::runif(J, -1, 1)
  U2_cc <- stats::runif(J, -1, 1)
  U3_cc <- stats::rnorm(J, 0, 1)
  U1.c <- rep(U1_cc, times = n)
  U2.c <- rep(U2_cc, times = n)
  U3.c <- rep(U3_cc, times = n)
  X1 <- stats::rnorm(N, U1.c, 0.5)
  X2 <- stats::rnorm(N, U2.c, 0.5)
  X1.c <- rep(stats::aggregate(X1, by = list(id), FUN = mean)$x, times = n)
  X2.c <- rep(stats::aggregate(X2, by = list(id), FUN = mean)$x, times = n)
  match_out <- match_me(U3_cc)
  A.c <- rep(match_out$A_c, times = n)
  pair <- rep(match_out$pairs, times = n)
  data.frame(id = id, U1.c = U1.c, U2.c = U2.c, U3.c = U3.c,
             X1.c = X1.c, X2.c = X2.c, X1 = X1, X2 = X2,
             A.c = A.c, pair = pair)
}


#' Main simulation DGP (Balzer 2021 § 5) — baseline + post-baseline missingness
#'
#' @param effect Logical. If `TRUE`, the intervention has an effect on Y.
#' @param n Integer vector of candidate cluster sizes (sampled with
#'   replacement to assign one to each cluster).
#' @param J Number of clusters (must be even for pair matching).
#' @param verbose Logical.
#'
#' @return A `data.frame` at the individual level with the observed
#'   variables plus counterfactual outcomes `Y.1`, `Y.0`. Pre-aggregated
#'   cluster covariates `X1.c`, `X2.c` are joined on.
#'
#' @keywords internal
gen_data_complex <- function(effect, n, J, verbose = FALSE) {
  n <- sample(rep(n, times = J), J)
  N <- sum(n)
  df <- gen_baseline_complex(n = n, J = J, N = N)

  X1 <- df$X1; X2 <- df$X2; X1.c <- df$X1.c; X2.c <- df$X2.c; U.c <- df$U3.c
  A_ind <- list(1, 0, df$A.c)

  M <- vector("list", 3)
  U.M <- stats::runif(N)
  for (i in seq_along(A_ind)) {
    p <- stats::plogis(
      -1 + 2 * A_ind[[i]] + 1 * (X1 + X2) +
        0.2 * (1 - A_ind[[i]]) * (X1.c + X2.c) + 0.25 * U.c
    )
    M[[i]] <- as.numeric(U.M < p)
  }

  beta_A <- if (effect) -2.5 else 0
  beta_M <- if (effect)  4   else 0
  Y <- vector("list", 3)
  U.Y <- stats::runif(N)
  for (i in seq_along(A_ind)) {
    p <- stats::plogis(
      1 + beta_A * A_ind[[i]] + beta_M * M[[i]] +
        0.5 * (X1 + X2) + 0.2 * (X1.c + X2.c) + 0.25 * U.c
    )
    Y[[i]] <- as.numeric(U.Y < p)
  }

  U.Delta <- stats::runif(N)
  txt <- stats::plogis( 3 - 3 * M[[3]] - 0.5 * (X1 + X2))
  con <- stats::plogis(-2 + 3 * M[[3]] + 0.5 * (X1 + X2))
  pscore <- A_ind[[3]] * txt + (1 - A_ind[[3]]) * con
  Delta <- as.numeric(U.Delta < pscore)

  O <- data.frame(
    id = df$id, pair = df$pair, U = 1, X1 = X1, X2 = X2,
    A = A_ind[[3]], M = M[[3]],
    Delta = Delta,
    Y = Y[[3]], Y.1 = Y[[1]], Y.0 = Y[[2]]
  )
  do_aggregation(O)
}


#' Supplementary simulation DGP (Balzer 2021 § 9.6) — baseline-only missingness
#'
#' @inheritParams gen_data_complex
#' @keywords internal
gen_data_simple <- function(effect, n, J, verbose = FALSE) {
  n <- sample(rep(n, times = J), J)
  N <- sum(n)
  id <- rep(seq_len(J), times = n)

  U1.c_short <- stats::runif(J, 1.75, 2.25)
  U2.c_short <- stats::rnorm(J, 0, 1)
  U3.c_short <- stats::rnorm(J, 0, 1)

  match_out <- match_me(U3.c_short)
  A.c <- match_out$A_c
  pair <- rep(match_out$pairs, times = n)
  A_ind <- list(1, 0, rep(A.c, times = n))

  U1.c <- rep(U1.c_short, times = n)
  U2.c <- rep(U2.c_short, times = n)
  U3.c <- rep(U3.c_short, times = n)

  X1 <- stats::rnorm(N, U1.c, 1)
  X2 <- stats::rnorm(N, U2.c, 1)
  X1.c <- rep(stats::aggregate(X1, by = list(id), FUN = mean)$x, times = n)
  X2.c <- rep(stats::aggregate(X2, by = list(id), FUN = mean)$x, times = n)

  if (effect) { beta_A <- 0.15; beta_AX <- 0.15 } else { beta_A <- 0; beta_AX <- 0 }
  U.Y <- stats::runif(N)
  Y <- vector("list", 3)
  for (i in seq_along(A_ind)) {
    p <- stats::plogis(
      -4 + beta_A * A_ind[[i]] + 0.3 * X1.c + 0.3 * X2.c +
        0.4 * X1 + 0.2 * X2 + 0.5 * X1.c * X1 +
        beta_AX * X1 * A_ind[[i]] + 0.3 * U3.c
    )
    Y[[i]] <- as.numeric(U.Y < p)
  }

  U.Delta <- stats::runif(N)
  pscore <- stats::plogis(
    4 - 0.25 * A_ind[[3]] - 0.5 * X1.c - 0.1 * X2.c -
      0.1 * X2 - 0.75 * X1 - 0.75 * X1 * A_ind[[3]]
  )
  Delta <- as.numeric(U.Delta < pscore)

  O <- data.frame(
    id = id, pair = pair, U = 1, X1 = X1, X2 = X2,
    A = A_ind[[3]],
    Delta = Delta,
    Y = Y[[3]], Y.1 = Y[[1]], Y.0 = Y[[2]]
  )
  do_aggregation(O)
}


#' Compute true treatment-specific means and effect measures from full counterfactual data
#' @keywords internal
get_truth <- function(X_all) {
  truth_C <- stats::aggregate(X_all[, c("U", "Y.1", "Y.0")],
                              by = list(X_all$id), FUN = sum)
  Y1.c <- truth_C$Y.1 / truth_C$U
  Y0.c <- truth_C$Y.0 / truth_C$U
  Risk1.C <- mean(Y1.c)
  Risk0.C <- mean(Y0.c)
  RD <- Risk1.C - Risk0.C
  aRR <- Risk1.C / Risk0.C
  km <- get_variability(X0 = Y0.c, X1 = Y1.c,
                        pairs = stats::aggregate(X_all$pair,
                                                 by = list(X_all$id),
                                                 FUN = mean)[, 2])
  data.frame(Risk1.C = Risk1.C, Risk0.C = Risk0.C,
             RD = RD, aRR = aRR,
             k.con = km$k.con, k.txt = km$k.txt, km = km$km)
}


#' Population truth via a large counterfactual sample
#' @keywords internal
get_truth_pop <- function(do_complex, effect, n, J) {
  X_all <- if (do_complex) gen_data_complex(effect, n, J)
           else gen_data_simple(effect, n, J)
  get_truth(X_all)
}


#' @keywords internal
get_variability <- function(X1, X0, pairs) {
  k.con <- stats::sd(X0) / mean(X0)
  k.txt <- stats::sd(X1) / mean(X1)
  temp <- unique(pairs)
  n_pairs <- length(temp)
  varm <- pi_ <- rep(NA_real_, n_pairs)
  for (i in seq_len(n_pairs)) {
    varm[i] <- stats::var(X0[pairs == temp[i]])
    pi_[i]  <- mean(X0[pairs == temp[i]])
  }
  km <- mean(sqrt(varm) / pi_)
  data.frame(k.con = k.con, k.txt = k.txt, km = km)
}


#' Observed data with missing outcomes masked
#' @keywords internal
get_obs <- function(O, exclude) {
  O$Y <- O$Y * O$Delta
  if (exclude) {
    O <- O[O$Delta == 1, ]
    O$Y.1 <- NULL; O$Y.0 <- NULL; O$X1.c <- NULL; O$X2.c <- NULL
    O <- do_aggregation(O)
  } else {
    O$Y.1 <- NULL; O$Y.0 <- NULL
  }
  O
}
