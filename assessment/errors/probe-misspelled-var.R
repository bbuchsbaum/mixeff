## Probe: misspelled-var — formula references a column not in data
## Compares lme4::lmer / lme4::glmer vs mixeff::lmm / mixeff::glmm

library(lme4)
library(mixeff)

# ---- shared data: y, x, g exist; 'xx' does NOT ----
set.seed(42)
df <- data.frame(
  y = rnorm(40),
  x = rnorm(40),
  g = factor(rep(seq_len(10), each = 4))
)

cat("=== lme4::lmer — misspelled fixed-effect predictor ===\n")
lme4_lmer_err <- tryCatch(
  lme4::lmer(y ~ xx + (1 | g), data = df),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(lme4_lmer_err), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(lme4_lmer_err), "\n\n")

cat("=== lme4::glmer — misspelled fixed-effect predictor ===\n")
lme4_glmer_err <- tryCatch(
  lme4::glmer(y ~ xx + (1 | g), data = df, family = gaussian()),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(lme4_glmer_err), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(lme4_glmer_err), "\n\n")

cat("=== mixeff::lmm — misspelled fixed-effect predictor ===\n")
mixeff_lmm_err <- tryCatch(
  mixeff::lmm(y ~ xx + (1 | g), data = df, control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mixeff_lmm_err), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(mixeff_lmm_err), "\n\n")

cat("=== mixeff::glmm — misspelled fixed-effect predictor ===\n")
mixeff_glmm_err <- tryCatch(
  mixeff::glmm(y ~ xx + (1 | g), data = df, family = gaussian(),
               control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mixeff_glmm_err), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(mixeff_glmm_err), "\n\n")

cat("=== lme4::lmer — misspelled grouping variable ===\n")
lme4_lmer_grp <- tryCatch(
  lme4::lmer(y ~ x + (1 | gg), data = df),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(lme4_lmer_grp), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(lme4_lmer_grp), "\n\n")

cat("=== mixeff::lmm — misspelled grouping variable ===\n")
mixeff_lmm_grp <- tryCatch(
  mixeff::lmm(y ~ x + (1 | gg), data = df, control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mixeff_lmm_grp), collapse = ", "), "\n")
cat("message:\n")
cat(conditionMessage(mixeff_lmm_grp), "\n\n")
