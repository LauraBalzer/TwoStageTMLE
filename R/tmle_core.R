# Core TMLE machinery for Stage 2: initial outcome regression, propensity
# score, fluctuation parameter, targeting step, and the do_tmle() wrapper.
#
# These functions assume the input data has already been normalised:
#   - column 'Y' for the (cluster-level or individual-level) outcome
#   - column 'A' for the treatment indicator
#   - column 'alpha' for analysis weights
#   - column 'U' = 1 as a dummy adjustment variable for the unadjusted estimator
#   - column 'id' for cluster identity
#
# Outputs are returned with augmented columns rather than mutated in place;
# data.table copies are used at entry to keep semantics clear.


#' Initial estimation of E(Y | A, W)
#'
#' Fits a regression of `Y` on `A` and the chosen adjustment variables and
#' returns predicted Qbar(A,W), Qbar(1,W), Qbar(0,W). The form of the
#' regression is controlled by `q_form`:
#'
#' * `"glm"` — main-terms logistic GLM (Stage 2 outcome is scaled to `[0,1]`)
#' * `"step"` — `stats::step()` on the main-terms GLM
#' * `"step.interaction"` — `step()` with full pairwise interactions
#'
#' The `"lasso"` / `"mars"` / `"mars.corP"` algorithms from `from_analysis/`
#' are deliberately omitted from the package because they relied on
#' SuperLearner-internal screeners (`screen.glmnet`, `screen.corP`) and on
#' SEARCH-trial-specific covariate forcing. Users who want those can fit a
#' Super Learner externally and pass in the resulting `q_out`.
#'
#' @param train Data with `Y`, `A`, `alpha`, and the columns named in `q_adj`.
#' @param q_adj Character vector of adjustment-variable column names, or `'U'`.
#' @param q_form One of `"glm"`, `"step"`, `"step.interaction"`. Default `"glm"`.
#' @param q_out Pre-fit regression object, e.g., from a previous training fold.
#'   If supplied, no new fit is performed.
#' @param verbose Logical.
#'
#' @return A list with `q_adj`, `q_form`, `glm_out` (the fit), and a
#'   `data.table` `train` augmented with columns `QbarAW`, `Qbar1W`, `Qbar0W`.
#' @keywords internal
do_init_qbar <- function(train, q_adj, q_form = "glm", q_out = NULL,
                         verbose = FALSE) {
  q_form <- match.arg(q_form, c("glm", "step", "step.interaction"))
  train <- as_dt_copy(train)
  if (is.null(q_adj)) q_adj <- "U"

  train_temp <- as.data.frame(train[, c(q_adj, "A", "Y"), with = FALSE])
  X1 <- X0 <- train_temp
  X1$A <- 1; X0$A <- 0

  if (is.null(q_out)) {
    q_out <- suppressWarnings(stats::glm(
      Y ~ ., family = "binomial", data = train_temp, weights = train$alpha
    ))
    if (q_form == "step") {
      q_out <- stats::step(q_out, direction = "both", trace = 0, k = 2)
    } else if (q_form == "step.interaction") {
      q_out <- stats::step(q_out, scope = Y ~ .^2, direction = "both",
                           trace = 0, k = 2)
    }
    vmsg(verbose, "do_init_qbar: %s fit complete", q_form)
  }

  QbarAW <- stats::predict(q_out, newdata = train_temp, type = "response")
  Qbar1W <- stats::predict(q_out, newdata = X1, type = "response")
  Qbar0W <- stats::predict(q_out, newdata = X0, type = "response")

  train[, `:=`(QbarAW = QbarAW, Qbar1W = Qbar1W, Qbar0W = Qbar0W)]
  list(q_adj = q_adj, q_form = q_form, glm_out = q_out, train = train)
}


#' Propensity-score estimation and clever covariate construction
#'
#' Fits P(A = 1 | W) using the adjustment variables in `g_adj` and the
#' algorithm `g_form`, bounds the predicted pscore to `[0.025, 0.975]`, and
#' adds the clever covariates `H.1W = I(A=1)/g(1|W)`, `H.0W = I(A=0)/g(0|W)`,
#' `H.AW = H.1W - H.0W` to `train`.
#'
#' The bounding is a safety net; for a randomised trial with the correct
#' design g(1|W) = 0.5 and the bound never bites.
#'
#' @inheritParams do_init_qbar
#' @param g_adj Character vector of adjustment variables for the pscore.
#' @param g_form One of `"glm"`, `"step"`, `"step.interaction"`.
#' @param p_out Pre-fit pscore object. If supplied, no new fit is performed.
#'
#' @return A list with `g_adj`, `g_form`, `p_out`, and `train` augmented with
#'   columns `pscore`, `H.1W`, `H.0W`, `H.AW`.
#' @keywords internal
get_clever_cov <- function(train, g_adj, g_form = "glm", p_out = NULL,
                           verbose = FALSE) {
  g_form <- match.arg(g_form, c("glm", "step", "step.interaction"))
  train <- as_dt_copy(train)
  if (is.null(g_adj)) g_adj <- "U"

  train_temp <- as.data.frame(train[, c(g_adj, "A"), with = FALSE])

  if (is.null(p_out)) {
    p_out <- suppressWarnings(stats::glm(
      A ~ ., family = "binomial", data = train_temp, weights = train$alpha
    ))
    if (g_form == "step") {
      p_out <- stats::step(p_out, direction = "both", trace = 0, k = 2)
    } else if (g_form == "step.interaction") {
      p_out <- stats::step(p_out, scope = A ~ .^2, direction = "both",
                           trace = 0, k = 2)
    }
    vmsg(verbose, "get_clever_cov: %s pscore fit complete", g_form)
  }

  pscore <- stats::predict(p_out, newdata = train_temp, type = "response")
  pscore <- pmin(pmax(pscore, 0.025), 0.975)
  A.train <- train$A
  H.1W <- A.train / pscore
  H.0W <- (1 - A.train) / (1 - pscore)
  H.AW <- H.1W - H.0W

  train[, `:=`(pscore = pscore, H.1W = H.1W, H.0W = H.0W, H.AW = H.AW)]
  list(g_adj = g_adj, g_form = g_form, p_out = p_out, train = train)
}


#' Fluctuation coefficient (TMLE targeting step)
#'
#' One-dimensional clever covariate is used for `goal = "RD"`; a two-
#' dimensional clever covariate is used for `"aRR"` and `"OR"` so that
#' relative effects can be targeted via the Delta method.
#'
#' Skips the update when within-arm variance is below `1e-4` (per the
#' `from_analysis/` code; see `questions_for_laura.md` Q1).
#'
#' @param train Data augmented by `do_init_qbar()` and `get_clever_cov()`.
#' @param goal One of `"RD"`, `"aRR"`, `"OR"`.
#' @param verbose Logical.
#' @return A named numeric vector of fluctuation coefficients.
#' @keywords internal
get_epsilon <- function(train, goal, verbose = FALSE) {
  goal <- match.arg(goal, c("RD", "aRR", "OR"))
  A <- train$A
  Y <- train$Y

  skip <- (stats::var(Y[A == 1]) < 1e-4) | (stats::var(Y[A == 0]) < 1e-4)

  if (goal == "RD") {
    if (!skip) {
      fit <- suppressWarnings(stats::glm(
        Y ~ -1 + offset(stats::qlogis(train$QbarAW)) + train$H.AW,
        family = "binomial", weights = train$alpha
      ))
      eps <- fit$coefficients
    } else {
      eps <- 0
    }
    names(eps) <- "H.AW"
  } else {
    if (!skip) {
      fit <- suppressWarnings(stats::glm(
        Y ~ -1 + offset(stats::qlogis(train$QbarAW)) + train$H.0W + train$H.1W,
        family = "binomial", weights = train$alpha
      ))
      eps <- fit$coefficients
    } else {
      eps <- c(0, 0)
    }
    names(eps) <- c("H.0W", "H.1W")
  }
  vmsg(verbose, "get_epsilon: eps = %s", paste(round(eps, 6), collapse = ", "))
  eps
}


#' Targeted predictions Qbar*(A,W), Qbar*(1,W), Qbar*(0,W)
#'
#' Updates the initial Qbar predictions using the fluctuation coefficient
#' `eps`. For `goal = "RD"`, the one-dimensional update is applied; otherwise
#' the two-dimensional update.
#'
#' @inheritParams get_epsilon
#' @param eps Output of [get_epsilon()].
#' @return `train` augmented with columns `QbarAW.star`, `Qbar1W.star`, `Qbar0W.star`.
#' @keywords internal
do_targeting <- function(train, eps, goal) {
  goal <- match.arg(goal, c("RD", "aRR", "OR"))
  train <- as_dt_copy(train)

  g1W <- train$pscore
  g0W <- 1 - g1W

  if (goal == "RD") {
    QbarAW_star <- stats::plogis(stats::qlogis(train$QbarAW) + eps * train$H.AW)
    Qbar1W_star <- stats::plogis(stats::qlogis(train$Qbar1W) + eps / g1W)
    Qbar0W_star <- stats::plogis(stats::qlogis(train$Qbar0W) - eps / g0W)
  } else {
    QbarAW_star <- stats::plogis(
      stats::qlogis(train$QbarAW) + eps["H.0W"] * train$H.0W + eps["H.1W"] * train$H.1W
    )
    Qbar1W_star <- stats::plogis(stats::qlogis(train$Qbar1W) + eps["H.1W"] / g1W)
    Qbar0W_star <- stats::plogis(stats::qlogis(train$Qbar0W) + eps["H.0W"] / g0W)
  }
  train[, `:=`(QbarAW.star = QbarAW_star,
               Qbar1W.star = Qbar1W_star,
               Qbar0W.star = Qbar0W_star)]
  train
}


#' Master Stage 2 TMLE point estimator
#'
#' Runs the four-step TMLE: (1) initial Qbar, (2) clever covariate, (3)
#' fluctuation, (4) targeted predictions and treatment-specific means. Output
#' is unscaled back to the original outcome scale using `scale_min` /
#' `scale_max`.
#'
#' @param train Data already including `Y` (scaled to `[0,1]`), `A`, `alpha`,
#'   `id`, `U`, and any adjustment columns.
#' @param goal,target,q_adj,q_form,g_adj,g_form,q_out,p_out,scale_min,scale_max,doing_cv,verbose
#'   See the corresponding arguments of [stage2_tmle()].
#' @return A list with the full fit metadata plus `R1`, `R0`, `var_R1`,
#'   `var_R0`, `var_pair`, `var_break`, and the augmented training data.
#' @keywords internal
do_tmle <- function(train, goal, target = "clust",
                    q_adj, q_form = "glm",
                    g_adj = NULL, g_form = "glm",
                    q_out = NULL, p_out = NULL,
                    scale_min = 0, scale_max = 1,
                    doing_cv = FALSE, verbose = FALSE) {
  goal <- match.arg(goal, c("RD", "aRR", "OR"))
  target <- match.arg(target, c("clust", "indv", "pre.post"))

  q <- do_init_qbar(train = train, q_adj = q_adj, q_form = q_form,
                    q_out = q_out, verbose = verbose)
  train <- q$train

  g <- get_clever_cov(train = train, g_adj = g_adj, g_form = g_form,
                      p_out = p_out, verbose = verbose)
  train <- g$train

  eps <- get_epsilon(train = train, goal = goal, verbose = verbose)
  train <- do_targeting(train = train, eps = eps, goal = goal)

  J <- length(unique(train$id))
  if (nrow(train) > J && target == "clust") {
    if (!doing_cv) vmsg(verbose, "data=indv; target=clust")
    R1 <- mean(
      stats::aggregate(data.frame(x = train$alpha * train$Qbar1W.star),
                       by = list(train$id), FUN = sum)$x
    )
    R0 <- mean(
      stats::aggregate(data.frame(x = train$alpha * train$Qbar0W.star),
                       by = list(train$id), FUN = sum)$x
    )
  } else {
    R1 <- mean(train$alpha * train$Qbar1W.star)
    R0 <- mean(train$alpha * train$Qbar0W.star)
  }

  R1 <- R1 * (scale_max - scale_min) + scale_min
  R0 <- R0 * (scale_max - scale_min) + scale_min

  variance_out <- get_ic_variance(
    goal = goal, target = target, vdata = train,
    R1 = R1, R0 = R0,
    scale_min = scale_min, scale_max = scale_max,
    doing_cv = doing_cv, verbose = verbose
  )

  list(
    train = train,
    q_adj = q$q_adj, q_form = q$q_form, q_out = q$glm_out,
    g_adj = g$g_adj, g_form = g$g_form, p_out = g$p_out,
    eps = eps, R1 = R1, R0 = R0,
    var_R1 = variance_out$var_R1, var_R0 = variance_out$var_R0,
    var_pair = variance_out$var_pair, var_break = variance_out$var_break
  )
}


#' Direct unadjusted treatment-specific means (no TMLE update)
#'
#' Convenience path equivalent to `do_tmle()` with `q_adj = g_adj = "U"` but
#' more numerically stable with rare outcomes. Does not incorporate weights;
#' should not be used when `target = "indv"` with cluster-level data.
#'
#' @keywords internal
do_unadjusted <- function(train, goal, scale_min = 0, scale_max = 1,
                          verbose = FALSE) {
  goal <- match.arg(goal, c("RD", "aRR", "OR"))
  train <- as_dt_copy(train)
  N <- nrow(train)

  Qbar1W <- R1 <- mean(train[A == 1, Y])
  Qbar0W <- R0 <- mean(train[A == 0, Y])
  QbarAW <- ifelse(train$A == 1, R1, R0)

  pscore <- mean(train$A)
  H.1W <- train$A / pscore
  H.0W <- (1 - train$A) / (1 - pscore)
  H.AW <- H.1W - H.0W

  train[, `:=`(
    QbarAW = QbarAW, Qbar1W = Qbar1W, Qbar0W = Qbar0W,
    pscore = pscore, H.1W = H.1W, H.0W = H.0W, H.AW = H.AW,
    QbarAW.star = QbarAW, Qbar1W.star = Qbar1W, Qbar0W.star = Qbar0W
  )]

  R1 <- R1 * (scale_max - scale_min) + scale_min
  R0 <- R0 * (scale_max - scale_min) + scale_min

  variance_out <- get_ic_variance(
    goal = goal, target = "indv", vdata = train,
    R1 = R1, R0 = R0,
    scale_min = scale_min, scale_max = scale_max,
    doing_cv = FALSE, verbose = verbose
  )

  list(
    train = train,
    q_adj = "U", q_form = NA_character_, q_out = NULL,
    g_adj = "U", g_form = NA_character_, p_out = NULL,
    eps = 0, R1 = R1, R0 = R0,
    var_R1 = variance_out$var_R1, var_R0 = variance_out$var_R0,
    var_pair = variance_out$var_pair, var_break = variance_out$var_break
  )
}
