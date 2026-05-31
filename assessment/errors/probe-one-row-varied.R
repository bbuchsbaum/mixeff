# Probe variant: one-row but with a non-constant y (to isolate "too few rows" vs
# "constant response" in the mixeff message).  The real structural problem is
# n=1 — not enough data to fit any model regardless of y variation.

library(lme4)
library(mixeff)

# One row, non-constant y impossible (scalar), but use y != 0 to avoid any
# "constant" short-circuit that might mask the true "too few rows" check.
# Actually with n=1 there is only one y value — so it IS constant by definition.
# Let's try two rows, same group (1 level), varied y — isolates the grouping issue.

two_rows_one_group <- data.frame(
  y       = c(1.0, 2.5),
  x       = c(0.1, 0.9),
  subject = factor(c("A", "A"))  # only 1 level
)

cat("========== lme4::lmer — 2 rows, 1 group ==========\n")
msg_lme4 <- tryCatch(
  lme4::lmer(y ~ x + (1 | subject), data = two_rows_one_group),
  error   = function(e) paste0("[ERROR] ", conditionMessage(e)),
  warning = function(w) paste0("[WARNING] ", conditionMessage(w))
)
cat(msg_lme4, "\n\n")

cat("========== mixeff::lmm — 2 rows, 1 group ==========\n")
msg_mixeff <- tryCatch(
  mixeff::lmm(y ~ x + (1 | subject), data = two_rows_one_group,
              control = mixeff::mm_control(verbose = -1)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("[ERROR class=", cls, "]\n", conditionMessage(e))
  },
  warning = function(w) paste0("[WARNING] ", conditionMessage(w))
)
cat(msg_mixeff, "\n\n")

# Also try pure one-row with non-zero varied — impossible but confirms message source
cat("========== mixeff::lmm — truly 1 row ==========\n")
one_row_varied <- data.frame(y = 3.7, x = 0.5, subject = factor("A"))
msg_one_row <- tryCatch(
  mixeff::lmm(y ~ x + (1 | subject), data = one_row_varied,
              control = mixeff::mm_control(verbose = -1)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("[ERROR class=", cls, "]\n", conditionMessage(e))
  },
  warning = function(w) paste0("[WARNING] ", conditionMessage(w))
)
cat(msg_one_row, "\n\n")

cat("Done.\n")
