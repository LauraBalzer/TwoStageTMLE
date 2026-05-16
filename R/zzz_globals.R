# Quiet R CMD check NOTEs about undefined globals from data.table NSE.
# All names here are column names referenced inside data.table expressions.
utils::globalVariables(c(
  # data.table operators (silenced via :=, .SD already imported, but kept here
  # for older R-CMD-check pedantry)
  ".", ".SD", ".N", ".I",
  # column names used in data.table NSE
  "id", "pair", "A", "Y", "U", "alpha",
  "nIndv", "nIndv_Y",
  "Delta",
  "QbarAW", "Qbar1W", "Qbar0W",
  "Qbar1W.star", "Qbar0W.star", "QbarAW.star",
  "H.1W", "H.0W", "H.AW", "pscore",
  "indpt_unit",
  # columns produced by get_inference() and other helpers used inline
  "ci_lo", "ci_hi", "se", "est", "tmle"
))
