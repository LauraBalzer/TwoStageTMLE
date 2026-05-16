# Bitwise parity with the Tier A deep-dive (TwoStageTMLE-main code).
# The 2021 paper-companion code always removes the Q-chosen variable from
# the pscore candidate set (see questions_for_laura.md Q14), so this test
# uses remove_pscore = TRUE and a candidate set without the unadjusted 'U'.

test_that("Tier A deep-dive: 6/6 bitwise parity for RD/aRR/OR x break/keep", {
  setup <- get_cluster_input()
  clust_in <- setup$data
  deep <- setup$deep
  clust_cov <- deep$clust.cov

  cand_no_u <- list(
    cand_q_adj = as.list(clust_cov),
    cand_q_form = rep("glm", length(clust_cov)),
    cand_g_adj = as.list(clust_cov),
    cand_g_form = rep("glm", length(clust_cov))
  )

  for (goal in c("RD", "aRR", "OR")) {
    for (bm in c(TRUE, FALSE)) {
      key <- sprintf("%s__%s", goal, ifelse(bm, "break_match", "keep_match"))
      psi_true <- switch(goal,
                         RD = deep$truth$RD,
                         aRR = deep$truth$aRR,
                         OR = NA)
      old <- deep$stage2[[key]]
      new <- quiet_stage2(
        data = clust_in, goal = goal, target = "clust", break_match = bm,
        do_data_adapt = TRUE, remove_pscore = TRUE,
        cand_q_adj = cand_no_u$cand_q_adj, cand_q_form = cand_no_u$cand_q_form,
        cand_g_adj = cand_no_u$cand_g_adj, cand_g_form = cand_no_u$cand_g_form,
        psi = psi_true
      )
      expect_parity(old, new)
    }
  }
})
