library(lme4)
library(mixeff)

set.seed(42)
n <- 60
df <- data.frame(
  y_char = sample(c("low", "mid", "high"), n, replace = TRUE),
  y_bin  = sample(c("yes", "no"), n, replace = TRUE),
  x      = rnorm(n),
  subj   = factor(rep(1:10, each = 6))
)

# ---- mixeff::glmm with character response (multi-class, binomial family) ----
cat("===== mixeff::glmm multi-class character response + binomial =====\n")
r1 <- tryCatch(
  mixeff::glmm(y_char ~ x + (1 | subj), data = df, family = binomial(),
               control = mm_control(verbose = -1L)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("ERROR [class: ", cls, "]: ", conditionMessage(e))
  }
)
cat(r1, "\n\n")

# ---- mixeff::glmm with binary character response + binomial family ----
cat("===== mixeff::glmm binary character response + binomial =====\n")
r2 <- tryCatch(
  mixeff::glmm(y_bin ~ x + (1 | subj), data = df, family = binomial(),
               control = mm_control(verbose = -1L)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("ERROR [class: ", cls, "]: ", conditionMessage(e))
  }
)
cat(r2, "\n\n")

# ---- Collect full condition object for lmm to see all fields ----
cat("===== mixeff::lmm full condition fields =====\n")
cond <- tryCatch(
  mixeff::lmm(y_char ~ x + (1 | subj), data = df,
              control = mm_control(verbose = -1L)),
  error = function(e) e
)
cat("class:", paste(class(cond), collapse=", "), "\n")
cat("message:", conditionMessage(cond), "\n")
extra <- setdiff(names(cond), c("message", "call", "trace", "parent", "use_cli_format", "rlang"))
if (length(extra)) cat("extra fields:", paste(extra, collapse=", "), "\n")
cat("\n")

# ---- Same for glmm binary ----
cat("===== mixeff::glmm binary char full condition =====\n")
cond2 <- tryCatch(
  mixeff::glmm(y_bin ~ x + (1 | subj), data = df, family = binomial(),
               control = mm_control(verbose = -1L)),
  error = function(e) e
)
if (inherits(cond2, "condition")) {
  cat("class:", paste(class(cond2), collapse=", "), "\n")
  cat("message:", conditionMessage(cond2), "\n")
  extra2 <- setdiff(names(cond2), c("message", "call", "trace", "parent", "use_cli_format", "rlang"))
  if (length(extra2)) cat("extra fields:", paste(extra2, collapse=", "), "\n")
} else {
  cat("Succeeded (no error):", class(cond2), "\n")
}

cat("\nDone.\n")
