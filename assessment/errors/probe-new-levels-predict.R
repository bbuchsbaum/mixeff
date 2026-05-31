## Probe: new-levels-predict
## Trigger predict(newdata) where newdata contains a grouping level not seen
## during training, without allow.new.levels = TRUE.
## Compare lme4::lmer + predict vs mixeff::lmm + predict.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

set.seed(42)

# ---- Training data: subjects 1-10, 4 obs each ----------------------------
n_subj  <- 10L
n_obs   <- 4L
subject <- factor(rep(seq_len(n_subj), each = n_obs))
x       <- rnorm(n_subj * n_obs)
y       <- 1.5 + 0.7 * x + rnorm(n_subj * n_obs, sd = 0.5) +
           rep(rnorm(n_subj, sd = 0.8), each = n_obs)
train   <- data.frame(y = y, x = x, subject = subject)

# ---- New data: one row from subject "99" (never seen during training) ----
# Note: lme4 and mixeff both only need predictor columns in newdata for
# predict(); the response `y` is not required. We omit it here to focus the
# probe on the new-level signal rather than a data-validation red herring.
newdata_new_level <- data.frame(
  x       = 0.5,
  subject = factor("99", levels = "99")   # completely new level
)

# ============================================================
# 1. lme4::lmer
# ============================================================
cat("=== lme4 ===\n")
fit_lme4 <- lmer(y ~ x + (1 | subject), data = train, REML = TRUE)

lme4_err <- tryCatch(
  predict(fit_lme4, newdata = newdata_new_level, allow.new.levels = FALSE),
  error   = function(e) e,
  warning = function(w) w
)

cat("lme4 condition class: ", paste(class(lme4_err), collapse = ", "), "\n")
cat("lme4 message:\n  ", conditionMessage(lme4_err), "\n\n")

# ============================================================
# 2. mixeff::lmm
# ============================================================
cat("=== mixeff ===\n")
fit_mm <- lmm(y ~ x + (1 | subject), data = train,
              control = mm_control(verbose = -1L))

mm_err <- tryCatch(
  predict(fit_mm, newdata = newdata_new_level, allow.new.levels = FALSE),
  error   = function(e) e,
  warning = function(w) w
)

cat("mixeff condition class: ", paste(class(mm_err), collapse = ", "), "\n")
cat("mixeff message:\n  ", conditionMessage(mm_err), "\n\n")

# ============================================================
# 3. Verify allow.new.levels = TRUE works (should succeed silently)
# ============================================================
cat("=== allow.new.levels = TRUE (should succeed) ===\n")

lme4_ok <- tryCatch(
  predict(fit_lme4, newdata = newdata_new_level, allow.new.levels = TRUE),
  error = function(e) e
)
cat("lme4  allow.new.levels=TRUE: ", if (is.numeric(lme4_ok)) round(lme4_ok, 4) else conditionMessage(lme4_ok), "\n")

mm_ok <- tryCatch(
  predict(fit_mm, newdata = newdata_new_level, allow.new.levels = TRUE),
  error = function(e) e
)
cat("mixeff allow.new.levels=TRUE:", if (is.numeric(mm_ok)) round(mm_ok, 4) else conditionMessage(mm_ok), "\n\n")

# ============================================================
# 4. Summary comparison
# ============================================================
cat("=== SUMMARY ===\n")
cat("Scenario: predict on newdata with an unseen grouping level\n")
cat("  (no allow.new.levels = TRUE supplied)\n\n")
cat("lme4 error class   :", paste(class(lme4_err), collapse = ", "), "\n")
cat("lme4 message       :", conditionMessage(lme4_err), "\n\n")
cat("mixeff error class :", paste(class(mm_err), collapse = ", "), "\n")
cat("mixeff message     :", conditionMessage(mm_err), "\n")
