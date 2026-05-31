## Probe: singular-fit scenario
## Compares lme4 and mixeff messages when a model fits singularly.
## A random-slope model on data with zero slope variation is a reliable trigger.

library(lme4)
library(mixeff)

## ---- helper to capture messages/warnings/errors verbatim -------------------
capture_all <- function(expr_thunk) {
  msgs   <- character()
  warns  <- character()
  err    <- NULL
  result <- withCallingHandlers(
    tryCatch(expr_thunk(), error = function(e) { err <<- e; NULL }),
    message = function(m) { msgs  <<- c(msgs,  conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  list(result = result, messages = msgs, warnings = warns, error = err)
}

## ---- Dataset that reliably induces a singular fit --------------------------
## Scenario A: random slope with zero between-subject variation in slope
## All subjects have the exact same x sequence → slope variance = 0
set.seed(42)
n_subj <- 10
n_obs  <- 4
subject <- factor(rep(seq_len(n_subj), each = n_obs))
x       <- rep(seq_len(n_obs), times = n_subj)          # identical within each subject
y       <- 1 + 0.5 * x + rnorm(n_subj * n_obs, sd = 0.5)

df_singular <- data.frame(y = y, x = x, subject = subject)

cat("======================================================\n")
cat("SCENARIO A: random slope on constant x — zero slope variance\n")
cat("======================================================\n\n")

## --- lme4 -------------------------------------------------------------------
cat("--- lme4::lmer ---\n")
cap_lme4 <- capture_all(function() {
  lme4::lmer(y ~ x + (1 + x | subject), data = df_singular, REML = TRUE)
})
cat("Warnings:\n")
if (length(cap_lme4$warnings)) {
  for (w in cap_lme4$warnings) cat("  [WARNING]", w, "\n")
} else {
  cat("  (none)\n")
}
cat("Messages:\n")
if (length(cap_lme4$messages)) {
  for (m in cap_lme4$messages) cat("  [MESSAGE]", m, "\n")
} else {
  cat("  (none)\n")
}
if (!is.null(cap_lme4$error)) {
  cat("Error:\n  [ERROR]", conditionMessage(cap_lme4$error), "\n")
}
if (!is.null(cap_lme4$result)) {
  cat("isSingular:", isSingular(cap_lme4$result), "\n")
  cat("VarCorr summary:\n")
  print(lme4::VarCorr(cap_lme4$result))
}

cat("\n--- lme4::lmer + lmerTest::as.lmerModLmerTest (post-fit) ---\n")
cap_lmertest <- capture_all(function() {
  fit <- lme4::lmer(y ~ x + (1 + x | subject), data = df_singular, REML = TRUE)
  lmerTest::as.lmerModLmerTest(fit)
})
if (length(cap_lmertest$warnings)) for (w in cap_lmertest$warnings) cat("  [WARNING]", w, "\n")

cat("\n--- mixeff::lmm ---\n")
cap_mixeff <- capture_all(function() {
  mixeff::lmm(y ~ x + (1 + x | subject), data = df_singular,
              control = mm_control(verbose = -1L))
})
cat("Warnings:\n")
if (length(cap_mixeff$warnings)) {
  for (w in cap_mixeff$warnings) cat("  [WARNING]", w, "\n")
} else {
  cat("  (none)\n")
}
cat("Messages:\n")
if (length(cap_mixeff$messages)) {
  for (m in cap_mixeff$messages) cat("  [MESSAGE]", m, "\n")
} else {
  cat("  (none)\n")
}
if (!is.null(cap_mixeff$error)) {
  cat("Error:\n  [ERROR]", conditionMessage(cap_mixeff$error), "\n")
  if (!is.null(attr(cap_mixeff$error, "reason_code"))) {
    cat("  reason_code:", attr(cap_mixeff$error, "reason_code"), "\n")
  }
  cat("  class:", paste(class(cap_mixeff$error), collapse = ", "), "\n")
}
if (!is.null(cap_mixeff$result)) {
  fit_mm <- cap_mixeff$result
  cat("fit_status:", fit_mm$fit_status, "\n")
  cat("is_singular:", is_singular(fit_mm), "\n")
  cat("\nprint(fit):\n")
  print(fit_mm)
  cat("\nVarCorr:\n")
  print(VarCorr(fit_mm))
}

## ---- Scenario B: intercept-only random effect with insufficient groups -----
## Very few groups → variance collapses to boundary = 0
cat("\n\n======================================================\n")
cat("SCENARIO B: very few groups → intercept variance at boundary\n")
cat("======================================================\n\n")

set.seed(7)
df_few <- data.frame(
  y       = rnorm(6),
  x       = 1:6,
  subject = factor(c(1,1,1,2,2,2))
)

cat("--- lme4::lmer (2 groups, random intercept) ---\n")
cap_lme4_b <- capture_all(function() {
  lme4::lmer(y ~ x + (1 | subject), data = df_few, REML = TRUE)
})
if (length(cap_lme4_b$warnings)) for (w in cap_lme4_b$warnings) cat("  [WARNING]", w, "\n")
if (!is.null(cap_lme4_b$error))   cat("  [ERROR]",   conditionMessage(cap_lme4_b$error), "\n")
if (!is.null(cap_lme4_b$result)) {
  cat("isSingular:", isSingular(cap_lme4_b$result), "\n")
  print(lme4::VarCorr(cap_lme4_b$result))
}

cat("\n--- mixeff::lmm (2 groups, random intercept) ---\n")
cap_mixeff_b <- capture_all(function() {
  mixeff::lmm(y ~ x + (1 | subject), data = df_few,
              control = mm_control(verbose = -1L))
})
if (length(cap_mixeff_b$warnings)) for (w in cap_mixeff_b$warnings) cat("  [WARNING]", w, "\n")
if (length(cap_mixeff_b$messages)) for (m in cap_mixeff_b$messages) cat("  [MESSAGE]", m, "\n")
if (!is.null(cap_mixeff_b$error)) {
  cat("  [ERROR]", conditionMessage(cap_mixeff_b$error), "\n")
  cat("  class:", paste(class(cap_mixeff_b$error), collapse = ", "), "\n")
}
if (!is.null(cap_mixeff_b$result)) {
  fit_b <- cap_mixeff_b$result
  cat("fit_status:", fit_b$fit_status, "\n")
  cat("is_singular:", is_singular(fit_b), "\n")
  cat("\nprint(fit):\n")
  print(fit_b)
  cat("\nVarCorr:\n")
  print(VarCorr(fit_b))
}

## ---- Scenario C: overparameterised random structure -------------------------
## (1 + x + z | subject) with 2 observations per subject → guaranteed singular
cat("\n\n======================================================\n")
cat("SCENARIO C: overparameterised random structure (rank > obs per group)\n")
cat("======================================================\n\n")

set.seed(99)
df_over <- data.frame(
  y       = rnorm(20),
  x       = rnorm(20),
  z       = rnorm(20),
  subject = factor(rep(1:10, each = 2))
)

cat("--- lme4::lmer (random slope + slope2 with 2 obs/group) ---\n")
cap_lme4_c <- capture_all(function() {
  lme4::lmer(y ~ x + z + (1 + x + z | subject), data = df_over, REML = TRUE)
})
if (length(cap_lme4_c$warnings)) for (w in cap_lme4_c$warnings) cat("  [WARNING]", w, "\n")
if (!is.null(cap_lme4_c$error))   cat("  [ERROR]",   conditionMessage(cap_lme4_c$error), "\n")
if (!is.null(cap_lme4_c$result)) {
  cat("isSingular:", isSingular(cap_lme4_c$result), "\n")
  print(lme4::VarCorr(cap_lme4_c$result))
}

cat("\n--- mixeff::lmm (random slope + slope2 with 2 obs/group) ---\n")
cap_mixeff_c <- capture_all(function() {
  mixeff::lmm(y ~ x + z + (1 + x + z | subject), data = df_over,
              control = mm_control(verbose = -1L))
})
if (length(cap_mixeff_c$warnings)) for (w in cap_mixeff_c$warnings) cat("  [WARNING]", w, "\n")
if (length(cap_mixeff_c$messages)) for (m in cap_mixeff_c$messages) cat("  [MESSAGE]", m, "\n")
if (!is.null(cap_mixeff_c$error)) {
  cat("  [ERROR]", conditionMessage(cap_mixeff_c$error), "\n")
  cat("  class:", paste(class(cap_mixeff_c$error), collapse = ", "), "\n")
}
if (!is.null(cap_mixeff_c$result)) {
  fit_c <- cap_mixeff_c$result
  cat("fit_status:", fit_c$fit_status, "\n")
  cat("is_singular:", is_singular(fit_c), "\n")
  cat("\nprint(fit):\n")
  print(fit_c)
  cat("\nVarCorr:\n")
  print(VarCorr(fit_c))
}

cat("\n\nDone.\n")
