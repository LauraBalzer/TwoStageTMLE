#' Influence-curve-based variance estimator for Stage 2 TMLE
#'
#' Computes the influence-curve components `DY1`, `DY0`, and the goal-scale
#' IC `DY`, then returns variance estimates on the original outcome scale.
#' For pair-matched designs, also returns the variance treating each pair as
#' the independent unit.
#'
#' **Bug fix:** the previous implementation in `from_analysis/stage2.R:249-250`
#' added `scale_min` when unscaling the IC components. That is wrong because
#' influence-curve components are residuals (`alpha * H * (Y - Q*)`); shifting
#' the residual by a constant is incorrect. The package multiplies by
#' `(scale_max - scale_min)` only. The bug was variance-invariant (because
#' `var(X + c) = var(X)`), so the change does not affect downstream SE / CI /
#' p-value. See `questions_for_laura.md` Q2.
#'
#' @param goal One of `"RD"`, `"aRR"`, `"OR"`.
#' @param target One of `"clust"`, `"indv"`, `"pre.post"`.
#' @param vdata Training data augmented with `alpha`, `H.1W`, `H.0W`, `Y`,
#'   `Qbar1W.star`, `Qbar0W.star`, `id`, and (for pair-matched) `pair`.
#' @param R1,R0 Treatment-specific means on the original outcome scale.
#' @param sample_effect Logical. `TRUE` (default) uses the sample-effect IC
#'   (no centering term); `FALSE` adds the centering term for the population
#'   effect.
#' @param scale_min,scale_max Outcome-rescaling parameters.
#' @param doing_cv Logical. Suppresses some print/verbose output when called
#'   from inside cross-validation.
#' @param verbose Logical.
#'
#' @return A list with components `R1`, `R0`, `DY1`, `var_R1`, `DY0`,
#'   `var_R0`, `DY`, `var_break`, `DY_paired`, `var_pair`.
#'
#' @keywords internal
get_ic_variance <- function(goal, target = "clust", vdata,
                            R1 = NA_real_, R0 = NA_real_,
                            sample_effect = TRUE,
                            scale_min = 0, scale_max = 1,
                            doing_cv = FALSE, verbose = FALSE) {
  goal <- match.arg(goal, c("RD", "aRR", "OR"))
  target <- match.arg(target, c("clust", "indv", "pre.post"))

  J <- length(unique(vdata$id))

  if (sample_effect) {
    DY1 <- vdata$alpha * vdata$H.1W * (vdata$Y - vdata$Qbar1W.star)
    DY0 <- vdata$alpha * vdata$H.0W * (vdata$Y - vdata$Qbar0W.star)
  } else {
    DY1 <- vdata$alpha * (
      vdata$H.1W * (vdata$Y - vdata$Qbar1W.star) + vdata$Qbar1W.star - R1
    )
    DY0 <- vdata$alpha * (
      vdata$H.0W * (vdata$Y - vdata$Qbar0W.star) + vdata$Qbar0W.star - R0
    )
  }

  # Unscale: multiplication only. (See bug-fix note above.)
  DY1 <- DY1 * (scale_max - scale_min)
  DY0 <- DY0 * (scale_max - scale_min)

  if (length(DY1) > J) {
    if (target == "clust") {
      if (!doing_cv) vmsg(verbose, "data=indv; target=clust")
      DY1 <- stats::aggregate(DY1, by = list(vdata$id), FUN = sum)[, -1]
      DY0 <- stats::aggregate(DY0, by = list(vdata$id), FUN = sum)[, -1]
    } else if (target == "indv") {
      if (!doing_cv) vmsg(verbose, "data=indv; target=indv")
      DY1 <- c(household_ic(as.matrix(DY1), id = vdata$id))
      DY0 <- c(household_ic(as.matrix(DY0), id = vdata$id))
    } else if (target == "pre.post") {
      if (!doing_cv) vmsg(verbose, "PRE POST")
      DY1 <- DY1[DY1 != 0]
      DY0 <- DY0[DY0 != 0]
    }
  }

  DY <- switch(goal,
    RD  = DY1 - DY0,
    aRR = (1 / R1) * DY1 - (1 / R0) * DY0,
    OR  = (1 / R1) * DY1 + (1 / (1 - R1)) * DY1
        - (1 / (1 - R0)) * DY0 - (1 / R0) * DY0
  )

  if (!doing_cv) vmsg(verbose, "Solve EIF: %g", mean(DY))

  var_R1 <- stats::var(DY1) / J
  var_R0 <- stats::var(DY0) / J
  var_break <- stats::var(DY) / J

  if ("pair" %in% names(vdata)) {
    pairC <- stats::aggregate(vdata, by = list(vdata$id), FUN = mean)$pair
    pairs <- unique(pairC)
    n_pairs <- length(pairs)
    DY_paired <- numeric(n_pairs)
    for (i in seq_len(n_pairs)) {
      these <- pairC %in% pairs[i]
      DY_paired[i] <- 0.5 * sum(DY[these])
    }
    var_pair <- stats::var(DY_paired) / n_pairs
  } else {
    DY_paired <- NA_real_
    var_pair <- NA_real_
  }

  list(R1 = R1, R0 = R0,
       DY1 = DY1, var_R1 = var_R1,
       DY0 = DY0, var_R0 = var_R0,
       DY = DY, var_break = var_break,
       DY_paired = DY_paired, var_pair = var_pair)
}
