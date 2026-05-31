## Probe: contrast/test_effect referencing a nonexistent coefficient name
## Scenario: contrast-nonexistent-coef
## Tests: lme4/lmer + lme4/glmer vs. mixeff lmm() + mixeff glmm()
## Also tests test_effect() for nonexistent term

library(lme4)
library(mixeff)

# ── Shared data ────────────────────────────────────────────────────────────────
set.seed(42)
n_subj <- 20
n_obs  <- 4
N      <- n_subj * n_obs

dat <- data.frame(
  y       = rnorm(N, mean = 10, sd = 2),
  x       = rnorm(N),
  subject = factor(rep(seq_len(n_subj), each = n_obs))
)
dat$y_bin <- as.integer(dat$y > 10)

cat("=== DATA SUMMARY ===\n")
cat(sprintf("N=%d, subjects=%d, obs/subject=%d\n\n", N, n_subj, n_obs))

# ── 1. lme4::lmer — contrast via linearHypothesis (not built-in contrast verb) ─
# lme4 has no first-class contrast() with named coefficient lookup.
# The canonical post-fit verbs are fixef() and multcomp::glht().
# We demonstrate what lme4 produces when a user tries the natural approaches.

cat("=== lme4::lmer — accessing nonexistent fixed-effect term ===\n")

fit_lme4 <- lmer(y ~ x + (1 | subject), data = dat)
cat("Fixed-effect names in lme4 fit:\n")
print(names(lme4::fixef(fit_lme4)))

# Approach 1: subscript fixef() with a wrong name
cat("\n--- lme4 fixef()['NONEXISTENT'] ---\n")
err_lme4_fixef <- tryCatch({
  val <- lme4::fixef(fit_lme4)["NONEXISTENT"]
  cat(sprintf("Result (no error): %s\n", format(val)))
  val
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", conditionMessage(e)))
  e
}, warning = function(w) {
  cat(sprintf("WARNING: %s\n", conditionMessage(w)))
  w
})

# Approach 2: manual contrast vector with wrong length (simulates "coefficient by name")
cat("\n--- lme4 manual contrast with wrong-dimension L (via multcomp-style) ---\n")
cat("(lme4 has no built-in named-coef contrast verb; NA subscript is silent)\n")

# ── 2. lme4::glmer — same scenario for GLMM ──────────────────────────────────
cat("\n=== lme4::glmer — accessing nonexistent fixed-effect term ===\n")

fit_glmer <- glmer(y_bin ~ x + (1 | subject), data = dat, family = binomial())
cat("Fixed-effect names in glmer fit:\n")
print(names(lme4::fixef(fit_glmer)))

cat("\n--- lme4 glmer fixef()['NONEXISTENT'] ---\n")
err_glmer_fixef <- tryCatch({
  val <- lme4::fixef(fit_glmer)["NONEXISTENT"]
  cat(sprintf("Result (no error): %s\n", format(val)))
  val
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", conditionMessage(e)))
  e
})

# ── 3. mixeff::lmm — contrast() with wrong-dimension L ───────────────────────
cat("\n=== mixeff::lmm — contrast() with incorrect column count ===\n")
fit_mm <- lmm(y ~ x + (1 | subject), data = dat)
cat("Fixed-effect names in mixeff lmm fit:\n")
print(names(fixef(fit_mm)))
cat(sprintf("Number of fixed effects: %d\n", length(fixef(fit_mm))))

# Wrong number of columns in L (too few)
cat("\n--- mixeff contrast() with L having wrong ncol (too few) ---\n")
err_mm_contrast_dim <- tryCatch({
  L_bad_dim <- matrix(c(1, 0), nrow = 1)  # 2 cols instead of 2 — actually this
  # is fine for 2-coef model; let's try 3 cols for a 2-coef model
  L_bad_dim <- matrix(c(1, 0, 0), nrow = 1)  # 3 cols for 2-coef model
  result <- contrast(fit_mm, L_bad_dim)
  cat("Result (no error): returned without error\n")
  print(result)
  result
}, error = function(e) {
  cat(sprintf("ERROR class: %s\n", paste(class(e), collapse=", ")))
  cat(sprintf("ERROR message: %s\n", conditionMessage(e)))
  e
})

# ── 4. mixeff::lmm — test_effect() with nonexistent term name ─────────────────
cat("\n=== mixeff::lmm — test_effect() with nonexistent term ===\n")
cat(sprintf("Available fixed-effect terms: %s\n",
            paste(mixeff:::mm_fixed_effect_terms(fit_mm), collapse=", ")))

cat("\n--- mixeff test_effect(fit, 'NONEXISTENT_TERM') ---\n")
err_mm_test_effect <- tryCatch({
  result <- test_effect(fit_mm, "NONEXISTENT_TERM")
  cat("Result (no error): returned without error\n")
  print(result)
  result
}, error = function(e) {
  cat(sprintf("ERROR class: %s\n", paste(class(e), collapse=", ")))
  cat(sprintf("ERROR message: %s\n", conditionMessage(e)))
  e
})

# ── 5. mixeff::lmm — confint() with nonexistent parm name ────────────────────
cat("\n=== mixeff::lmm — confint() with nonexistent parm name ===\n")
err_mm_confint <- tryCatch({
  result <- confint(fit_mm, parm = "NONEXISTENT_COEF", method = "wald")
  cat("Result (no error): returned without error\n")
  print(result)
  result
}, error = function(e) {
  cat(sprintf("ERROR class: %s\n", paste(class(e), collapse=", ")))
  cat(sprintf("ERROR message: %s\n", conditionMessage(e)))
  e
})

# ── 6. mixeff::glmm — summary() with post-fit verb using bad coef name ────────
cat("\n=== mixeff::glmm — summary() coefficients table (normal) ===\n")
fit_glmm <- tryCatch({
  glmm(y_bin ~ x + (1 | subject), data = dat, family = binomial())
}, error = function(e) {
  cat(sprintf("glmm() fit ERROR: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(fit_glmm)) {
  cat("Fixed-effect names in mixeff glmm fit:\n")
  # glmm summary coefficients
  s <- summary(fit_glmm, tests = "coefficients")
  print(s)

  # glmm has no contrast() method — verify what happens
  cat("\n--- mixeff contrast() dispatched to glmm (no method registered) ---\n")
  err_glmm_contrast <- tryCatch({
    L_ok <- matrix(c(0, 1), nrow = 1)
    result <- contrast(fit_glmm, L_ok)
    cat("Result (no error): returned without error\n")
    print(result)
    result
  }, error = function(e) {
    cat(sprintf("ERROR class: %s\n", paste(class(e), collapse=", ")))
    cat(sprintf("ERROR message: %s\n", conditionMessage(e)))
    e
  })

  # Try a bad-name scenario via glmm — if there is a method for glmm contrast
  cat("\n--- mixeff contrast() on glmm with wrong-dim L ---\n")
  err_glmm_contrast_dim <- tryCatch({
    L_bad <- matrix(c(1, 0, 0), nrow = 1)  # 3 cols for 2-coef glmm
    result <- contrast(fit_glmm, L_bad)
    cat("Result (no error): returned without error\n")
    print(result)
    result
  }, error = function(e) {
    cat(sprintf("ERROR class: %s\n", paste(class(e), collapse=", ")))
    cat(sprintf("ERROR message: %s\n", conditionMessage(e)))
    e
  })
}

# ── 7. mixeff::lmm — contrast() with a named vector (named approach) ──────────
# The contrast() function takes a numeric L matrix, not names.
# But what if someone passes a named vector where a name doesn't match?
# This is the closest analog: passing the right-size L but with a misleading
# rowname that doesn't match any coefficient.
# The real "nonexistent coef name" scenario is test_effect() — checked in #4.
# Let's also check what happens when someone tries to build L by name:
cat("\n=== mixeff::lmm — contrast() with named L, wrong rowname (semantic test) ===\n")
cat("(The L matrix is size-checked by column count, not by name)\n")
err_mm_named_L <- tryCatch({
  L_named <- matrix(c(0, 1), nrow = 1)
  rownames(L_named) <- "NONEXISTENT_COEF_LABEL"
  colnames(L_named) <- c("(Intercept)", "x")
  result <- contrast(fit_mm, L_named)
  cat("Result (no error): returned, contrast is valid numerically\n")
  cat(sprintf("Contrast rowname used: %s\n", rownames(result$L)))
  print(result$table[, c("contrast", "estimate", "p_value", "status"), drop=FALSE])
  result
}, error = function(e) {
  cat(sprintf("ERROR class: %s\n", paste(class(e), collapse=", ")))
  cat(sprintf("ERROR message: %s\n", conditionMessage(e)))
  e
})

cat("\n=== DONE ===\n")
