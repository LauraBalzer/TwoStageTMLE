# Bitwise parity with the locked Tier B unit fixtures.
# Each fixture exercises one code path in quiet_stage2() and was generated
# from from_analysis/ with the documented bug fixes monkey-patched in.

test_that("Tier B fixtures: 13/13 bitwise parity", {
  setup <- get_cluster_input()
  clust_in <- setup$data
  deep <- setup$deep
  ref <- readRDS(fixture_path("tierB_unit_fixtures.rds"))

  cand <- cand_adj(c("X1.c", "X2.c"))

  # B1: RD, break match, APS, target = clust (paper-default path)
  expect_parity(ref$B1_RD_break_aps_clust, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$RD
  ))

  # B2: aRR, keep match, APS
  expect_parity(ref$B2_aRR_keep_aps, quiet_stage2(
    data = clust_in, goal = "aRR", target = "clust", break_match = FALSE,
    do_data_adapt = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$aRR
  ))

  # B3: OR, break match, APS
  expect_parity(ref$B3_OR_break_aps, quiet_stage2(
    data = clust_in, goal = "OR", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = NA
  ))

  # B4: non-adaptive, explicit q_adj / g_adj
  expect_parity(ref$B4_RD_break_fixed, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = FALSE, q_adj = "X1.c", g_adj = "X2.c",
    psi = deep$truth$RD
  ))

  # B5: unadjusted shortcut
  expect_parity(ref$B5_RD_break_unadj, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    unadjusted = TRUE, psi = deep$truth$RD
  ))

  # B6: one-sided test, alt_smaller = TRUE
  expect_parity(ref$B6_RD_break_one_sided, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE, one_sided = TRUE, alt_smaller = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$RD
  ))

  # B7: do_cv_variance — produces extra columns; reuse the basic est columns
  expect_parity(ref$B7_RD_break_cvvar, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE, do_cv_variance = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$RD
  ))

  # B8: return_ic = TRUE — comparison via the est sub-element
  expect_parity(ref$B8_RD_break_returnIC, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE, return_ic = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$RD
  ))

  # B9: remove_pscore = TRUE
  expect_parity(ref$B9_RD_break_remove_pscore, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE, remove_pscore = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$RD
  ))

  # B10: Qform = "step"
  expect_parity(ref$B10_RD_step, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE,
    cand_q_adj = list(c("X1.c", "X2.c")), cand_q_form = "step",
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$RD
  ))

  # B11: Qform = "step.interaction"
  expect_parity(ref$B11_RD_step_interaction, quiet_stage2(
    data = clust_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE,
    cand_q_adj = list(c("X1.c", "X2.c")), cand_q_form = "step.interaction",
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form,
    psi = deep$truth$RD
  ))

  # B12: V-fold CV (60 indpt units, V = 5)
  big_in <- rbind(clust_in, clust_in)
  big_in[, id := seq_len(.N)]
  big_in[, pair := rep(seq_len(.N / 2), each = 2)]
  expect_parity(ref$B12_RD_break_vfold, quiet_stage2(
    data = big_in, goal = "RD", target = "clust", break_match = TRUE,
    do_data_adapt = TRUE, V = 5,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form
  ))

  # B13: target = "indv" with cluster-level data
  clust_indv <- data.table::copy(clust_in)
  J <- nrow(clust_indv); N_tot <- sum(clust_indv$nIndv)
  clust_indv[, alpha := J / N_tot * nIndv]
  expect_parity(ref$B13_RD_target_indv, quiet_stage2(
    data = clust_indv, goal = "RD", target = "indv", break_match = TRUE,
    do_data_adapt = TRUE,
    cand_q_adj = cand$cand_q_adj, cand_q_form = cand$cand_q_form,
    cand_g_adj = cand$cand_g_adj, cand_g_form = cand$cand_g_form
  ))
})
