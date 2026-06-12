# List-based result helpers. stage2_tmle() returns a one-row data.table (kept
# for parity with the original Stage2); these convert one-or-more such results
# into the plain nested list the applied pipelines emit, and flatten lists back
# to a data.frame for report tables / golden-master comparison.


#' Convert Stage-2 estimate row(s) into a tidy per-estimator list ("leaf")
#'
#' Turns the one-row `data.table` returned by [stage2_tmle()] — or several such
#' rows, one per goal (e.g. `"RD"` and `"aRR"`) — into a plain named list:
#' arm means under `arms$txt` / `arms$con`, one effect block per goal (named by
#' the goal), and the selected adjustment indices under `adjustment`.
#'
#' @param ... One or more one-row results from [stage2_tmle()]. Names are used
#'   as the goal labels; if unnamed, each result's own `goal` is used when
#'   present, else `g1`, `g2`, ...
#'
#' @return A plain `list` with `arms` (`txt`, `con` numeric vectors) and one
#'   numeric vector per goal (`est`, `ci_lo`, `ci_hi`, `se`, `pval`, `q_adj`),
#'   where `q_adj` is the APS-selected outcome-adjustment index **for that goal**
#'   (APS is run per goal, so it can differ between, say, RD and aRR). Arm means
#'   are taken from the first result (they do not depend on the goal).
#' @export
as_stage2_list <- function(...) {
  ests <- list(...)
  if (!length(ests)) stop("supply at least one stage2_tmle() result")
  ests <- lapply(ests, as.data.frame)
  goal_names <- names(ests)
  if (is.null(goal_names)) goal_names <- rep("", length(ests))
  for (i in seq_along(ests)) {
    if (!nzchar(goal_names[i])) {
      g <- ests[[i]][["goal"]]
      goal_names[i] <- if (!is.null(g)) as.character(g)[1] else paste0("g", i)
    }
  }

  first <- ests[[1]]
  out <- list(
    arms = list(
      txt = c(est = first$txt_est, ci_lo = first$txt_ci_lo,
              ci_hi = first$txt_ci_hi, se = first$txt_se),
      con = c(est = first$con_est, ci_lo = first$con_ci_lo,
              ci_hi = first$con_ci_hi, se = first$con_se)
    )
  )
  for (i in seq_along(ests)) {
    e <- ests[[i]]
    # q_adj is per goal: APS is run separately for each goal, so the selected
    # adjustment index can differ between e.g. RD and aRR.
    out[[goal_names[i]]] <- c(est = e$est, ci_lo = e$ci_lo, ci_hi = e$ci_hi,
                              se = e$se, pval = e$pval, q_adj = e$q_adj)
  }
  out
}


#' Flatten a nested Stage-2 result list into a data.frame
#'
#' Recursively walks a nested list of leaves produced by [as_stage2_list()]
#' (e.g. endpoint -> subgroup -> estimator -> leaf) and returns one row per leaf,
#' with the list path as leading label columns and the leaf's arm means, per-goal
#' effects, and adjustment indices as the remaining columns. Useful for report
#' tables and for comparing against golden masters.
#'
#' @param x A nested list whose leaves are [as_stage2_list()] outputs (a leaf is
#'   recognised by having an `arms` element).
#' @param label_names Optional character vector naming the path levels (e.g.
#'   `c("endpoint","subgroup","estimator")`); used as the label column names.
#'
#' @return A `data.frame`, one row per leaf.
#' @export
flatten_stage2 <- function(x, label_names = NULL) {
  rows <- list()
  walk <- function(node, path) {
    if (is.list(node) && !is.null(node[["arms"]])) {       # a leaf
      leaf <- node
      flat <- c(
        txt_est = leaf$arms$txt[["est"]], txt_ci_lo = leaf$arms$txt[["ci_lo"]],
        txt_ci_hi = leaf$arms$txt[["ci_hi"]],
        con_est = leaf$arms$con[["est"]], con_ci_lo = leaf$arms$con[["ci_lo"]],
        con_ci_hi = leaf$arms$con[["ci_hi"]]
      )
      goals <- setdiff(names(leaf), "arms")
      for (g in goals) {
        v <- leaf[[g]]
        names(v) <- paste0(g, "_", names(v))   # e.g. RD_est, RD_q_adj, aRR_est, ...
        flat <- c(flat, v)                      # numeric: arm means + per-goal effects + q_adj
      }
      labs <- as.list(path)
      names(labs) <- paste0("label", seq_along(path))
      rows[[length(rows) + 1L]] <<- c(labs, as.list(flat))
      return(invisible())
    }
    if (is.list(node)) {
      nms <- names(node); if (is.null(nms)) nms <- as.character(seq_along(node))
      for (k in seq_along(node)) walk(node[[k]], c(path, nms[k]))
    }
  }
  walk(x, character(0))
  if (!length(rows)) return(data.frame())

  df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  depth <- sum(grepl("^label[0-9]+$", names(df)))   # leading path/label columns
  if (!is.null(label_names) && length(label_names) == depth) {
    names(df)[seq_len(depth)] <- label_names
  }
  df
}
