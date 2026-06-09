# TwoStageTMLE

R implementation of the Two-Stage Targeted Minimum Loss-Based Estimator
(TMLE) for cluster randomized trials, following Balzer *et al.* (2021).

* **Stage 1** — within each cluster, an individual-level TMLE adjusts for
  differential outcome measurement on baseline `W` and post-baseline `M`
  covariates.
* **Stage 2** — a cluster-level TMLE compares treated to control clusters
  with Adaptive Pre-specification choosing the adjustment variable that
  maximises empirical efficiency, while preserving Type-I error control.

Supports risk difference, arithmetic risk ratio, and odds-ratio targets;
pair-matched and unmatched designs; cluster-level and individual-level
target effects (Benitez *et al.* 2023 weighting).

## Installation

```r
# install.packages("remotes")
remotes::install_github("LauraBalzer/TwoStageTMLE")
```

## Quick start

```r
library(TwoStageTMLE)
library(data.table)

# Simulate one CRT (internal helper used by the package tests)
set.seed(1)
X_all <- TwoStageTMLE:::gen_data_complex(effect = TRUE, n = c(100, 150, 200), J = 30)
O <- TwoStageTMLE:::get_obs(X_all, exclude = FALSE)

# Run the full Two-Stage TMLE. `two_stage_tmle()` is an internal reference
# implementation of the Balzer (2021) driver; the public API is the general
# building blocks (stage1_mean/stage1_apply + add_weights + stage2_tmle), shown
# composed in vignette("two-stage-tmle").
out <- TwoStageTMLE:::two_stage_tmle(
  data = O,
  ind_cov = c("X1", "X2", "M"),
  clust_cov = c("X1.c", "X2.c"),
  goal = "RD", target = "clust", break_match = FALSE,
  do_data_adapt = TRUE,
  cand_q_adj = cand_adj(c("X1.c", "X2.c"))$cand_q_adj,
  cand_q_form = cand_adj(c("X1.c", "X2.c"))$cand_q_form,
  cand_g_adj = cand_adj(c("X1.c", "X2.c"))$cand_g_adj,
  cand_g_form = cand_adj(c("X1.c", "X2.c"))$cand_g_form,
  remove_pscore = TRUE
)
out$est
```

For details, see `vignette("two-stage-tmle")` and
`vignette("reproduce-paper")`.

## Main functions

The package exposes general, composable building blocks; specific published
methods (Balzer 2021, Nakato 2026) are documented in the vignettes as
compositions of these (and kept as internal reference implementations).

* `stage1_mean()` — one measurement-adjusted cluster mean
* `stage1_ratio()` — ratio of two cluster means
* `stage1_apply()` — general engine: map any per-cluster estimator over clusters
* `stage2_tmle()` — cluster-level TMLE with Adaptive Pre-specification
* `add_weights()` — compute analysis weights per Benitez *et al.* (2023)
* `cand_adj()` — build candidate adjustment specifications for APS
* `as_stage2_list()`, `flatten_stage2()` — list-based result helpers

## References

* Balzer L.B. *et al.* (2021). Two-Stage TMLE to Reduce Bias and Improve
  Efficiency in Cluster Randomized Trials.
  <https://arxiv.org/abs/2106.15737>
* Balzer L.B., van der Laan M., Petersen M., SEARCH Consortium (2016).
  Adaptive pre-specification in randomized trials with and without
  pair-matching. *Statistics in Medicine* 35(10):4528-4545.
* Benitez A. *et al.* (2023). Comparative methods for the analysis of
  cluster randomized trials. arXiv:2110.09633.

## License

GPL (>= 3).
