## Probe: na-in-predictor — NAs in a predictor with default na.action
## Compares lme4::lmer / lme4::glmer vs mixeff::lmm / mixeff::glmm
## Also checks compile_model() and the lmm() + glmm() entry points.

library(lme4)
library(mixeff)

# ---- shared data: y, x, g exist; x has 3 NAs ----
set.seed(42)
n <- 40
df <- data.frame(
  y = rnorm(n),
  x = rnorm(n),
  g = factor(rep(seq_len(10), each = 4))
)
df$x[c(2, 15, 30)] <- NA   # inject NAs into fixed-effect predictor

# Also a version with NA in the response
df_na_y <- df
df_na_y$y[c(5, 20)] <- NA
df_na_y$x <- rnorm(n)       # clean x

# And a version with NA in the grouping variable
df_na_g <- df
df_na_g$x <- rnorm(n)       # clean x
df_na_g$g[c(7, 22)] <- NA

cat("============================================================\n")
cat("SCENARIO A: NA in fixed-effect predictor (x has 3 NAs)\n")
cat("============================================================\n\n")

cat("--- lme4::lmer ---\n")
lme4_lmer_err <- tryCatch(
  withCallingHandlers(
    lme4::lmer(y ~ x + (1 | g), data = df),
    warning = function(w) { message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") }
  ),
  error = function(e) e
)
if (inherits(lme4_lmer_err, "condition")) {
  cat("class:", paste(class(lme4_lmer_err), collapse = ", "), "\n")
  cat("message:\n  ", conditionMessage(lme4_lmer_err), "\n\n")
} else {
  cat("(succeeded — lme4 silently dropped NA rows via na.omit)\n")
  cat("nobs:", lme4::getME(lme4_lmer_err, "n"), "(expected", n - 3, "after NA drop)\n\n")
}

cat("--- lme4::glmer (binomial on 0/1 y) ---\n")
df_bin <- df
df_bin$y <- as.integer(df$y > 0)
lme4_glmer_err <- tryCatch(
  withCallingHandlers(
    lme4::glmer(y ~ x + (1 | g), data = df_bin, family = binomial()),
    warning = function(w) { message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") }
  ),
  error = function(e) e
)
if (inherits(lme4_glmer_err, "condition")) {
  cat("class:", paste(class(lme4_glmer_err), collapse = ", "), "\n")
  cat("message:\n  ", conditionMessage(lme4_glmer_err), "\n\n")
} else {
  cat("(succeeded — lme4 silently dropped NA rows via na.omit)\n")
  cat("nobs:", lme4::getME(lme4_glmer_err, "n"), "(expected", n - 3, "after NA drop)\n\n")
}

cat("--- mixeff::compile_model (pre-fit, NA in x) ---\n")
mx_compile_err <- tryCatch(
  mixeff::compile_model(y ~ x + (1 | g), data = df),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mx_compile_err), collapse = ", "), "\n")
cat("message:\n  ", conditionMessage(mx_compile_err), "\n\n")

cat("--- mixeff::lmm (NA in x) ---\n")
mx_lmm_err <- tryCatch(
  mixeff::lmm(y ~ x + (1 | g), data = df, control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mx_lmm_err), collapse = ", "), "\n")
cat("message:\n  ", conditionMessage(mx_lmm_err), "\n\n")

cat("--- mixeff::glmm (binomial, NA in x) ---\n")
mx_glmm_err <- tryCatch(
  mixeff::glmm(y ~ x + (1 | g), data = df_bin, family = binomial(),
               control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mx_glmm_err), collapse = ", "), "\n")
cat("message:\n  ", conditionMessage(mx_glmm_err), "\n\n")

cat("============================================================\n")
cat("SCENARIO B: NA in response variable (y has 2 NAs)\n")
cat("============================================================\n\n")

cat("--- lme4::lmer (NA in y) ---\n")
lme4_lmer_na_y <- tryCatch(
  withCallingHandlers(
    lme4::lmer(y ~ x + (1 | g), data = df_na_y),
    warning = function(w) { message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") }
  ),
  error = function(e) e
)
cat("class:", paste(class(lme4_lmer_na_y), collapse = ", "), "\n")
if (inherits(lme4_lmer_na_y, "condition")) {
  cat("message:\n  ", conditionMessage(lme4_lmer_na_y), "\n\n")
} else {
  cat("(succeeded — lme4 silently dropped NA rows)\n")
  cat("nobs:", lme4::getME(lme4_lmer_na_y, "n"), "(expected", n - 2, "after NA drop)\n\n")
}

cat("--- mixeff::lmm (NA in y) ---\n")
mx_lmm_na_y <- tryCatch(
  mixeff::lmm(y ~ x + (1 | g), data = df_na_y, control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mx_lmm_na_y), collapse = ", "), "\n")
if (inherits(mx_lmm_na_y, "condition")) {
  cat("message:\n  ", conditionMessage(mx_lmm_na_y), "\n\n")
} else {
  cat("(succeeded — mixeff did NOT catch NA in y)\n\n")
}

cat("============================================================\n")
cat("SCENARIO C: NA in grouping variable (g has 2 NAs)\n")
cat("============================================================\n\n")

cat("--- lme4::lmer (NA in g) ---\n")
lme4_lmer_na_g <- tryCatch(
  withCallingHandlers(
    lme4::lmer(y ~ x + (1 | g), data = df_na_g),
    warning = function(w) { message("WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") }
  ),
  error = function(e) e
)
cat("class:", paste(class(lme4_lmer_na_g), collapse = ", "), "\n")
if (inherits(lme4_lmer_na_g, "condition")) {
  cat("message:\n  ", conditionMessage(lme4_lmer_na_g), "\n\n")
} else {
  cat("(succeeded — lme4 silently dropped NA rows)\n")
  cat("nobs:", lme4::getME(lme4_lmer_na_g, "n"), "(expected", n - 2, "after NA drop)\n\n")
}

cat("--- mixeff::lmm (NA in g) ---\n")
mx_lmm_na_g <- tryCatch(
  mixeff::lmm(y ~ x + (1 | g), data = df_na_g, control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mx_lmm_na_g), collapse = ", "), "\n")
if (inherits(mx_lmm_na_g, "condition")) {
  cat("message:\n  ", conditionMessage(mx_lmm_na_g), "\n\n")
} else {
  cat("(succeeded — mixeff did NOT catch NA in g)\n\n")
}

cat("============================================================\n")
cat("SCENARIO D: glmm() custom na.action refusal\n")
cat("============================================================\n\n")

cat("--- mixeff::glmm with custom na.action=na.pass ---\n")
mx_glmm_na_action <- tryCatch(
  mixeff::glmm(y ~ x + (1 | g), data = df_bin, family = binomial(),
               na.action = na.pass,
               control = mm_control(verbose = -1)),
  error   = function(e) e,
  warning = function(w) w
)
cat("class:", paste(class(mx_glmm_na_action), collapse = ", "), "\n")
cat("message:\n  ", conditionMessage(mx_glmm_na_action), "\n\n")

cat("Done.\n")
