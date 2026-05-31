## probe-malformed-formula.R
## Probe how lme4 and mixeff handle syntactically broken / malformed formulas.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

set.seed(42)
df <- data.frame(
  y       = rnorm(60),
  x       = rnorm(60),
  subject = factor(rep(1:10, each = 6))
)

# Binary outcome for glmm
df$yb <- as.integer(df$y > 0)

capture_result <- function(expr) {
  tryCatch(
    withCallingHandlers(
      { val <- expr; list(ok = TRUE, value = val) },
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) invokeRestart("muffleMessage")
    ),
    error = function(e) {
      list(
        ok      = FALSE,
        message = conditionMessage(e),
        class   = class(e)
      )
    }
  )
}

print_result <- function(label, r) {
  if (r$ok) {
    cat(label, "=> (no error — returned value of class:", paste(class(r$value), collapse=", "), ")\n")
  } else {
    cat(label, "\n  message:", r$message, "\n  class:", paste(r$class, collapse=", "), "\n")
  }
}

separator <- function(title) {
  cat("\n", strrep("=", 70), "\n")
  cat("Scenario:", title, "\n")
  cat(strrep("=", 70), "\n\n")
}

# ---------------------------------------------------------------------------
# Scenario 1: Completely empty formula string
# ---------------------------------------------------------------------------
separator("Empty formula string")

print_result("lme4::lmer",  capture_result(lmer("", data = df)))
print_result("lme4::glmer", capture_result(glmer("", data = df, family = binomial())))
print_result("mixeff::lmm", capture_result(lmm("", data = df, control = mm_control(verbose = -1))))
print_result("mixeff::glmm",capture_result(glmm("", data = df, family = binomial(), control = mm_control(verbose = -1))))

# ---------------------------------------------------------------------------
# Scenario 2: Garbled string — double tilde
# ---------------------------------------------------------------------------
separator("Garbled string (double tilde: 'y ~~ x + (1 | subject)')")

garbled <- "y ~~ x + (1 | subject)"
print_result("lme4::lmer",  capture_result(lmer(garbled, data = df)))
print_result("lme4::glmer", capture_result(glmer(garbled, data = df, family = binomial())))
print_result("mixeff::lmm", capture_result(lmm(garbled, data = df, control = mm_control(verbose = -1))))
print_result("mixeff::glmm",capture_result(glmm(garbled, data = df, family = binomial(), control = mm_control(verbose = -1))))

# ---------------------------------------------------------------------------
# Scenario 3: Missing '|' in random term  (string form)
# ---------------------------------------------------------------------------
separator("Missing '|' in random term: 'y ~ x + (1  subject)'")

bad_re <- "y ~ x + (1  subject)"   # missing pipe
# Note: as.formula() would fail at parse time, so lmer must take character or
# the formula object. lmer() can take a character formula string directly.
print_result("lme4::lmer",  capture_result(lmer(bad_re, data = df)))
print_result("lme4::glmer", capture_result(glmer(bad_re, data = df, family = binomial())))
print_result("mixeff::lmm", capture_result(lmm(bad_re, data = df, control = mm_control(verbose = -1))))
print_result("mixeff::glmm",capture_result(glmm(bad_re, data = df, family = binomial(), control = mm_control(verbose = -1))))

# ---------------------------------------------------------------------------
# Scenario 4: Empty grouping factor '(1 | )'
# ---------------------------------------------------------------------------
separator("Empty grouping factor: 'y ~ x + (1 | )'")

empty_group <- "y ~ x + (1 | )"
print_result("lme4::lmer",  capture_result(lmer(empty_group, data = df)))
print_result("lme4::glmer", capture_result(glmer(empty_group, data = df, family = binomial())))
print_result("mixeff::lmm", capture_result(lmm(empty_group, data = df, control = mm_control(verbose = -1))))
print_result("mixeff::glmm",capture_result(glmm(empty_group, data = df, family = binomial(), control = mm_control(verbose = -1))))

# ---------------------------------------------------------------------------
# Scenario 5: Valid R formula but no random effects at all
# ---------------------------------------------------------------------------
separator("Valid R formula but no random effects: y ~ x")

no_re <- y ~ x
print_result("lme4::lmer",  capture_result(lmer(no_re, data = df)))
print_result("lme4::glmer", capture_result(glmer(yb ~ x, data = df, family = binomial())))
print_result("mixeff::lmm", capture_result(lmm(no_re, data = df, control = mm_control(verbose = -1))))
print_result("mixeff::glmm",capture_result(glmm(yb ~ x, data = df, family = binomial(), control = mm_control(verbose = -1))))

# ---------------------------------------------------------------------------
# Scenario 6: Non-formula, non-character input (numeric)
# ---------------------------------------------------------------------------
separator("Non-formula, non-character argument (numeric 42)")

print_result("lme4::lmer",  capture_result(lmer(42, data = df)))
print_result("lme4::glmer", capture_result(glmer(42, data = df, family = binomial())))
print_result("mixeff::lmm", capture_result(lmm(42, data = df, control = mm_control(verbose = -1))))
print_result("mixeff::glmm",capture_result(glmm(42, data = df, family = binomial(), control = mm_control(verbose = -1))))

# ---------------------------------------------------------------------------
# Scenario 7: Formula as NA
# ---------------------------------------------------------------------------
separator("Formula as NA")

print_result("lme4::lmer",  capture_result(lmer(NA, data = df)))
print_result("mixeff::lmm", capture_result(lmm(NA, data = df, control = mm_control(verbose = -1))))
print_result("mixeff::glmm",capture_result(glmm(NA, data = df, family = binomial(), control = mm_control(verbose = -1))))

cat("\n=== Done ===\n")
