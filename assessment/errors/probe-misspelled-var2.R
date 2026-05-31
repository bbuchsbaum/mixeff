## Probe addendum: glmm with a supported family + misspelled variable
library(lme4)
library(mixeff)

set.seed(42)
df <- data.frame(
  y = rbinom(40, 1, 0.5),
  x = rnorm(40),
  g = factor(rep(seq_len(10), each = 4))
)

cat("=== lme4::glmer (binomial) — misspelled fixed predictor ===\n")
lme4_glmer_bin <- tryCatch(
  lme4::glmer(y ~ xx + (1 | g), data = df, family = binomial()),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(lme4_glmer_bin), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(lme4_glmer_bin), "\n\n")

cat("=== mixeff::glmm (binomial) — misspelled fixed predictor ===\n")
mixeff_glmm_bin <- tryCatch(
  mixeff::glmm(y ~ xx + (1 | g), data = df, family = binomial(),
               control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mixeff_glmm_bin), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(mixeff_glmm_bin), "\n\n")

cat("=== mixeff::glmm (binomial) — misspelled grouping variable ===\n")
mixeff_glmm_grp <- tryCatch(
  mixeff::glmm(y ~ x + (1 | gg), data = df, family = binomial(),
               control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mixeff_glmm_grp), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(mixeff_glmm_grp), "\n\n")
