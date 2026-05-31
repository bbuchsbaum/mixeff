## probe-mismatched-weights.R
## Scenario: weights vector has wrong length (too short, too long, zero-length)
## Compare lme4::lmer / lme4::glmer vs mixeff::lmm / mixeff::glmm

library(lme4)
library(mixeff)

## ------------------------------------------------------------------
## Minimal reproducible data
## ------------------------------------------------------------------
set.seed(42)
n <- 40
df <- data.frame(
  y       = rnorm(n),
  x       = rnorm(n),
  subject = factor(rep(seq_len(10), each = 4)),
  count   = rpois(n, lambda = 3)
)

w_correct   <- runif(n, 0.5, 1.5)   # length 40  -- correct
w_short     <- w_correct[1:20]       # length 20  -- too short
w_long      <- c(w_correct, w_correct) # length 80 -- too long
w_zero      <- numeric(0)            # length 0   -- empty

capture_msg <- function(expr) {
  tryCatch(
    withCallingHandlers(
      { expr; "--- NO ERROR / NO WARNING ---" },
      warning = function(w) {
        invokeRestart("muffleWarning")
      }
    ),
    error   = function(e) paste0("[ERROR] ", conditionMessage(e)),
    warning = function(w) paste0("[WARNING] ", conditionMessage(w))
  )
}

cat("========================================\n")
cat("SCENARIO: mismatched-weights (wrong length)\n")
cat("========================================\n\n")

## ------------------------------------------------------------------
## 1. lme4::lmer — too-short weights
## ------------------------------------------------------------------
cat("--- lme4::lmer  (weights too short: 20, need 40) ---\n")
msg <- capture_msg(lmer(y ~ x + (1 | subject), data = df, weights = w_short))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 2. lme4::lmer — too-long weights
## ------------------------------------------------------------------
cat("--- lme4::lmer  (weights too long: 80, need 40) ---\n")
msg <- capture_msg(lmer(y ~ x + (1 | subject), data = df, weights = w_long))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 3. lme4::lmer — zero-length weights
## ------------------------------------------------------------------
cat("--- lme4::lmer  (weights zero-length) ---\n")
msg <- capture_msg(lmer(y ~ x + (1 | subject), data = df, weights = w_zero))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 4. lme4::glmer — too-short weights  (binomial)
## ------------------------------------------------------------------
df$ybin <- as.integer(df$y > 0)
cat("--- lme4::glmer (weights too short: 20, need 40) ---\n")
msg <- capture_msg(glmer(ybin ~ x + (1 | subject), data = df,
                         family = binomial(), weights = w_short))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 5. mixeff::lmm — too-short weights
## ------------------------------------------------------------------
cat("--- mixeff::lmm (weights too short: 20, need 40) ---\n")
msg <- capture_msg(lmm(y ~ x + (1 | subject), data = df,
                       weights = w_short, control = mm_control(verbose = -1)))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 6. mixeff::lmm — too-long weights
## ------------------------------------------------------------------
cat("--- mixeff::lmm (weights too long: 80, need 40) ---\n")
msg <- capture_msg(lmm(y ~ x + (1 | subject), data = df,
                       weights = w_long, control = mm_control(verbose = -1)))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 7. mixeff::lmm — zero-length weights
## ------------------------------------------------------------------
cat("--- mixeff::lmm (weights zero-length) ---\n")
msg <- capture_msg(lmm(y ~ x + (1 | subject), data = df,
                       weights = w_zero, control = mm_control(verbose = -1)))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 8. mixeff::glmm — weights (reserved, any non-NULL should error)
## ------------------------------------------------------------------
cat("--- mixeff::glmm (weights non-NULL — reserved param) ---\n")
msg <- capture_msg(glmm(ybin ~ x + (1 | subject), data = df,
                        family = binomial(), weights = w_short,
                        control = mm_control(verbose = -1)))
cat(msg, "\n\n")

## ------------------------------------------------------------------
## 9. Condition class inspection for mixeff errors
## ------------------------------------------------------------------
cat("--- mixeff::lmm condition class (too-short weights) ---\n")
cls <- tryCatch(
  lmm(y ~ x + (1 | subject), data = df, weights = w_short,
      control = mm_control(verbose = -1)),
  error = function(e) class(e)
)
cat(paste(cls, collapse = ", "), "\n\n")

cat("--- mixeff::glmm condition class (weights non-NULL) ---\n")
cls2 <- tryCatch(
  glmm(ybin ~ x + (1 | subject), data = df,
       family = binomial(), weights = w_short,
       control = mm_control(verbose = -1)),
  error = function(e) class(e)
)
cat(paste(cls2, collapse = ", "), "\n\n")

cat("Done.\n")
