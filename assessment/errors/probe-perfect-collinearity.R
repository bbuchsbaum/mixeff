## probe-perfect-collinearity.R
## Compare lme4 vs mixeff error messages when two fixed-effect predictors
## are perfectly collinear (x2 = 2 * x1).

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

set.seed(42)
n_subj <- 10
n_obs  <- 4

df <- data.frame(
  subject = factor(rep(seq_len(n_subj), each = n_obs)),
  x1      = rnorm(n_subj * n_obs),
  y       = rnorm(n_subj * n_obs)
)
df$x2 <- 2 * df$x1   # perfect collinearity

cat("=== DATA SUMMARY ===\n")
cat(sprintf("nrow=%d, cor(x1,x2)=%.6f\n", nrow(df), cor(df$x1, df$x2)))
cat("x2 = 2 * x1 exactly (perfect linear dependence)\n\n")

# ---- lme4 lmer ---------------------------------------------------------------
cat("=== lme4::lmer — perfect collinearity ===\n")
lmer_result <- tryCatch(
  {
    fit <- lme4::lmer(y ~ x1 + x2 + (1 | subject), data = df, REML = TRUE)
    list(fit = fit, warning = NULL, error = NULL)
  },
  warning = function(w) list(fit = NULL, warning = conditionMessage(w), error = NULL),
  error   = function(e) list(fit = NULL, warning = NULL, error = conditionMessage(e))
)

# lmer often warns but fits (rank-deficient); capture all
lmer_fit_result <- withCallingHandlers(
  tryCatch(
    lme4::lmer(y ~ x1 + x2 + (1 | subject), data = df, REML = TRUE),
    error = function(e) { cat("lme4 ERROR:", conditionMessage(e), "\n"); NULL }
  ),
  warning = function(w) {
    cat("lme4 WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("lme4 MESSAGE:", conditionMessage(m), "\n")
    invokeRestart("muffleMessage")
  }
)

if (!is.null(lmer_fit_result)) {
  cat("lme4 fit succeeded (rank-deficient?)\n")
  coef_vals <- lme4::fixef(lmer_fit_result)
  cat("fixef:", paste(names(coef_vals), round(coef_vals, 4), sep = "=", collapse = ", "), "\n")
  cat("isSingular:", lme4::isSingular(lmer_fit_result), "\n")
  cat("class:", paste(class(lmer_fit_result), collapse = ", "), "\n")
}

cat("\n")

# ---- lme4 glmer (Gaussian) ---------------------------------------------------
cat("=== lme4::glmer — perfect collinearity ===\n")
glmer_fit_result <- withCallingHandlers(
  tryCatch(
    lme4::glmer(y ~ x1 + x2 + (1 | subject), data = df, family = gaussian()),
    error = function(e) { cat("lme4::glmer ERROR:", conditionMessage(e), "\n"); NULL }
  ),
  warning = function(w) {
    cat("lme4::glmer WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("lme4::glmer MESSAGE:", conditionMessage(m), "\n")
    invokeRestart("muffleMessage")
  }
)

if (!is.null(glmer_fit_result)) {
  cat("lme4::glmer fit succeeded\n")
  coef_vals2 <- lme4::fixef(glmer_fit_result)
  cat("fixef:", paste(names(coef_vals2), round(coef_vals2, 4), sep = "=", collapse = ", "), "\n")
}

cat("\n")

# ---- mixeff lmm --------------------------------------------------------------
cat("=== mixeff::lmm — perfect collinearity ===\n")
mm_lmm_result <- withCallingHandlers(
  tryCatch(
    mixeff::lmm(y ~ x1 + x2 + (1 | subject), data = df,
                control = mixeff::mm_control(verbose = -1)),
    error = function(e) {
      cat("mixeff::lmm ERROR class:", paste(class(e), collapse = ", "), "\n")
      cat("mixeff::lmm ERROR message:", conditionMessage(e), "\n")
      NULL
    }
  ),
  warning = function(w) {
    cat("mixeff::lmm WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("mixeff::lmm MESSAGE:", conditionMessage(m), "\n")
    invokeRestart("muffleMessage")
  }
)

if (!is.null(mm_lmm_result)) {
  cat("mixeff::lmm fit succeeded\n")
  coef_mm <- mixeff::fixef(mm_lmm_result)
  cat("fixef:", paste(names(coef_mm), round(coef_mm, 4), sep = "=", collapse = ", "), "\n")
}

cat("\n")

# ---- mixeff glmm (Binomial) --------------------------------------------------
# Build a binary outcome for the glmm probe
set.seed(42)
df_bin <- df
df_bin$ybin <- rbinom(nrow(df_bin), 1, 0.5)

cat("=== lme4::glmer — perfect collinearity (binomial/logit) ===\n")
glmer_bin_result <- withCallingHandlers(
  tryCatch(
    lme4::glmer(ybin ~ x1 + x2 + (1 | subject), data = df_bin, family = binomial()),
    error = function(e) { cat("lme4::glmer(binomial) ERROR:", conditionMessage(e), "\n"); NULL }
  ),
  warning = function(w) {
    cat("lme4::glmer(binomial) WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("lme4::glmer(binomial) MESSAGE:", conditionMessage(m), "\n")
    invokeRestart("muffleMessage")
  }
)
if (!is.null(glmer_bin_result)) {
  coef_glmer_bin <- lme4::fixef(glmer_bin_result)
  cat("lme4::glmer(binomial) fit succeeded\n")
  cat("fixef:", paste(names(coef_glmer_bin), round(coef_glmer_bin, 4), sep = "=", collapse = ", "), "\n")
}

cat("\n=== mixeff::glmm — perfect collinearity (binomial/logit) ===\n")
mm_glmm_result <- withCallingHandlers(
  tryCatch(
    mixeff::glmm(ybin ~ x1 + x2 + (1 | subject), data = df_bin, family = binomial(),
                 control = mixeff::mm_control(verbose = -1)),
    error = function(e) {
      cat("mixeff::glmm ERROR class:", paste(class(e), collapse = ", "), "\n")
      cat("mixeff::glmm ERROR message:", conditionMessage(e), "\n")
      NULL
    }
  ),
  warning = function(w) {
    cat("mixeff::glmm WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("mixeff::glmm MESSAGE:", conditionMessage(m), "\n")
    invokeRestart("muffleMessage")
  }
)

if (!is.null(mm_glmm_result)) {
  cat("mixeff::glmm fit succeeded\n")
  coef_glmm <- mixeff::fixef(mm_glmm_result)
  cat("fixef:", paste(names(coef_glmm), round(coef_glmm, 4), sep = "=", collapse = ", "), "\n")
}

cat("\n")

# ---- Also probe integer collinearity (x2 = x1, identical) -------------------
cat("=== BONUS: identical predictors (x2 = x1) ===\n")
df2 <- df
df2$x2 <- df2$x1

cat("-- lme4::lmer identical predictors --\n")
withCallingHandlers(
  tryCatch(
    lme4::lmer(y ~ x1 + x2 + (1 | subject), data = df2, REML = TRUE),
    error = function(e) { cat("lme4 ERROR:", conditionMessage(e), "\n"); NULL }
  ),
  warning = function(w) {
    cat("lme4 WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }
)

cat("-- mixeff::lmm identical predictors --\n")
withCallingHandlers(
  tryCatch(
    mixeff::lmm(y ~ x1 + x2 + (1 | subject), data = df2,
                control = mixeff::mm_control(verbose = -1)),
    error = function(e) {
      cat("mixeff::lmm ERROR class:", paste(class(e), collapse = ", "), "\n")
      cat("mixeff::lmm ERROR message:", conditionMessage(e), "\n")
      NULL
    }
  ),
  warning = function(w) {
    cat("mixeff::lmm WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }
)

cat("\nDone.\n")
