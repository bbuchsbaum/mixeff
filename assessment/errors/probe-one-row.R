# Probe: single-row data frame — error message quality comparison
# Tests both lmm()/glmm() (mixeff) vs lmer()/glmer() (lme4)

library(lme4)
library(mixeff)

# ---- 1. Single-row data frame ------------------------------------------------
one_row <- data.frame(
  y       = 1.5,
  x       = 0.3,
  subject = factor("A")
)

cat("========== lme4::lmer — one-row ==========\n")
lme4_lmer_msg <- tryCatch(
  lme4::lmer(y ~ x + (1 | subject), data = one_row),
  error   = function(e) paste0("[ERROR] ", conditionMessage(e)),
  warning = function(w) paste0("[WARNING] ", conditionMessage(w))
)
cat(lme4_lmer_msg, "\n\n")

cat("========== mixeff::lmm — one-row ==========\n")
mixeff_lmm_msg <- tryCatch(
  mixeff::lmm(y ~ x + (1 | subject), data = one_row,
              control = mixeff::mm_control(verbose = -1)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("[ERROR class=", cls, "]\n", conditionMessage(e))
  },
  warning = function(w) paste0("[WARNING] ", conditionMessage(w))
)
cat(mixeff_lmm_msg, "\n\n")

# ---- 2. glmm / glmer counterpart --------------------------------------------
cat("========== lme4::glmer — one-row (binomial) ==========\n")
one_row_bin <- data.frame(
  y       = 1L,
  x       = 0.3,
  subject = factor("A")
)
lme4_glmer_msg <- tryCatch(
  lme4::glmer(y ~ x + (1 | subject), data = one_row_bin, family = binomial),
  error   = function(e) paste0("[ERROR] ", conditionMessage(e)),
  warning = function(w) paste0("[WARNING] ", conditionMessage(w))
)
cat(lme4_glmer_msg, "\n\n")

cat("========== mixeff::glmm — one-row (binomial) ==========\n")
mixeff_glmm_msg <- tryCatch(
  mixeff::glmm(y ~ x + (1 | subject), data = one_row_bin, family = binomial,
               control = mixeff::mm_control(verbose = -1)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("[ERROR class=", cls, "]\n", conditionMessage(e))
  },
  warning = function(w) paste0("[WARNING] ", conditionMessage(w))
)
cat(mixeff_glmm_msg, "\n\n")

cat("Done.\n")
