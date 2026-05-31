## Probe: slope-no-variation
##
## Scenario: a random slope is requested for a predictor that does not vary
## within any group.  The slope is structurally unidentifiable from the data.
##
## We build a dataset where `treatment` is a BETWEEN-subjects factor: every
## subject sees only ONE value of treatment, so (1 + treatment | subject) has
## no within-subject variation to estimate the slope from.
##
## We test four surfaces:
##   1. lme4::lmer()
##   2. lme4::glmer()
##   3. mixeff::lmm()
##   4. mixeff::glmm()

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

## ---- data setup ------------------------------------------------------------
set.seed(42)
n_subjects <- 20
# Between-subjects: each subject is in exactly ONE treatment condition.
subject    <- factor(rep(seq_len(n_subjects), each = 4))
treatment  <- factor(rep(c("A", "B"), each = n_subjects * 2))
y_cont     <- rnorm(n_subjects * 4)
y_bin      <- rbinom(n_subjects * 4, 1, 0.5)
df <- data.frame(subject = subject, treatment = treatment,
                 y_cont = y_cont, y_bin = y_bin)

# Verify: treatment does NOT vary within subject
within_var <- tapply(as.integer(df$treatment), df$subject, function(x) length(unique(x)))
stopifnot(all(within_var == 1L))
cat("Confirmed: treatment does not vary within any subject (all within_var == 1)\n\n")

## ---- helper ----------------------------------------------------------------
capture_msg <- function(expr_call) {
  err <- warn <- NULL
  val <- tryCatch(
    withCallingHandlers(
      eval(expr_call),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }
  )
  list(
    result  = val,
    error   = err,
    warning = warn
  )
}

## ===========================================================================
## 1. lme4::lmer  — random slope for between-subjects factor
## ===========================================================================
cat("=== lme4::lmer(y_cont ~ 1 + (1 + treatment | subject)) ===\n")
r_lmer <- capture_msg(quote(
  lmer(y_cont ~ 1 + (1 + treatment | subject), data = df,
       control = lmerControl(optimizer = "bobyqa"))
))
cat("  error  : ", if (is.null(r_lmer$error))   "(none)"    else r_lmer$error,   "\n")
cat("  warning: ", if (is.null(r_lmer$warning))  "(none)"    else paste(r_lmer$warning, collapse = " | "), "\n")
cat("  result : ", if (is.null(r_lmer$result))   "(NULL)"    else class(r_lmer$result)[1], "\n\n")

## ===========================================================================
## 2. lme4::glmer  — random slope for between-subjects factor (binary)
## ===========================================================================
cat("=== lme4::glmer(y_bin ~ 1 + (1 + treatment | subject), family=binomial) ===\n")
r_glmer <- capture_msg(quote(
  glmer(y_bin ~ 1 + (1 + treatment | subject), data = df,
        family = binomial,
        control = glmerControl(optimizer = "bobyqa"))
))
cat("  error  : ", if (is.null(r_glmer$error))   "(none)"    else r_glmer$error,   "\n")
cat("  warning: ", if (is.null(r_glmer$warning))  "(none)"    else paste(r_glmer$warning, collapse = " | "), "\n")
cat("  result : ", if (is.null(r_glmer$result))   "(NULL)"    else class(r_glmer$result)[1], "\n\n")

## ===========================================================================
## 3. mixeff::lmm  — random slope for between-subjects factor
## ===========================================================================
cat("=== mixeff::lmm(y_cont ~ 1 + (1 + treatment | subject)) ===\n")
r_lmm <- capture_msg(quote(
  lmm(y_cont ~ 1 + (1 + treatment | subject), data = df,
      control = mm_control(verbose = -1L))
))
cat("  error  : ", if (is.null(r_lmm$error))   "(none)"    else r_lmm$error,   "\n")
cat("  warning: ", if (is.null(r_lmm$warning))  "(none)"    else paste(r_lmm$warning, collapse = " | "), "\n")
cat("  result : ", if (is.null(r_lmm$result))   "(NULL)"    else class(r_lmm$result)[1], "\n")
if (!is.null(r_lmm$result)) {
  cat("  diagnostics:\n")
  print(diagnostics(r_lmm$result))
  cat("  fit_status:", fit_status(r_lmm$result), "\n")
}
cat("\n")

## ===========================================================================
## 4. mixeff::glmm  — random slope for between-subjects factor (binary)
## ===========================================================================
cat("=== mixeff::glmm(y_bin ~ 1 + (1 + treatment | subject), family='bernoulli') ===\n")
r_glmm <- capture_msg(quote(
  glmm(y_bin ~ 1 + (1 + treatment | subject), data = df,
       family = "bernoulli",
       control = mm_control(verbose = -1L))
))
cat("  error  : ", if (is.null(r_glmm$error))   "(none)"    else r_glmm$error,   "\n")
cat("  warning: ", if (is.null(r_glmm$warning))  "(none)"    else paste(r_glmm$warning, collapse = " | "), "\n")
cat("  result : ", if (is.null(r_glmm$result))   "(NULL)"    else class(r_glmm$result)[1], "\n")
if (!is.null(r_glmm$result)) {
  cat("  diagnostics:\n")
  print(diagnostics(r_glmm$result))
  cat("  fit_status:", fit_status(r_glmm$result), "\n")
}
cat("\n")

## ===========================================================================
## 5. Also test numeric predictor that is constant within group (pure numeric)
##    to ensure the audit also catches numeric slopes with no within-group SD
## ===========================================================================
cat("=== Numeric constant-within-group slope ===\n")
df2 <- df
# Give each subject the same numeric score (= their subject id, so between-subj)
df2$score <- as.numeric(df2$subject)
cat("-- mixeff::lmm with numeric between-subjects predictor --\n")
r_lmm2 <- capture_msg(quote(
  lmm(y_cont ~ 1 + (1 + score | subject), data = df2,
      control = mm_control(verbose = -1L))
))
cat("  error  : ", if (is.null(r_lmm2$error))   "(none)"    else r_lmm2$error,   "\n")
cat("  warning: ", if (is.null(r_lmm2$warning))  "(none)"    else paste(r_lmm2$warning, collapse = " | "), "\n")
cat("  result : ", if (is.null(r_lmm2$result))   "(NULL)"    else class(r_lmm2$result)[1], "\n")
if (!is.null(r_lmm2$result)) {
  cat("  diagnostics:\n")
  print(diagnostics(r_lmm2$result))
}
cat("\n")

cat("-- lme4::lmer with numeric between-subjects predictor --\n")
r_lmer2 <- capture_msg(quote(
  lmer(y_cont ~ 1 + (1 + score | subject), data = df2,
       control = lmerControl(optimizer = "bobyqa"))
))
cat("  error  : ", if (is.null(r_lmer2$error))   "(none)"    else r_lmer2$error,   "\n")
cat("  warning: ", if (is.null(r_lmer2$warning))  "(none)"    else paste(r_lmer2$warning, collapse = " | "), "\n")
cat("  result : ", if (is.null(r_lmer2$result))   "(NULL)"    else class(r_lmer2$result)[1], "\n\n")

cat("Done.\n")
