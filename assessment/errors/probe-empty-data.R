## probe-empty-data.R
## Probe error messages when fitting LMM/GLMM with a zero-row data frame.
## Tests both lme4/glmer and mixeff lmm()/glmm().

library(lme4)
library(mixeff)

# ---- helpers ----------------------------------------------------------------
capture_msg <- function(expr) {
  tryCatch(
    withCallingHandlers(
      { force(expr); list(type = "success", message = NA_character_, class = NA_character_) },
      warning = function(w) {
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(
        type    = "error",
        message = conditionMessage(e),
        class   = paste(class(e), collapse = ", ")
      )
    }
  )
}

# ---- build an empty data frame with the right structure ---------------------
df_full <- data.frame(
  y       = rnorm(40),
  x       = rep(0:3, 10),
  subject = factor(rep(seq_len(10), each = 4))
)
df_empty <- df_full[0L, ]   # zero rows, correct columns and types

cat("=== Data frame structure ===\n")
cat("nrow(df_empty):", nrow(df_empty), "\n")
cat("columns:", paste(names(df_empty), collapse=", "), "\n\n")

# ============================================================
# 1. lme4::lmer  — LMM
# ============================================================
cat("--- lme4::lmer (LMM, empty data) ---\n")
r_lmer <- capture_msg(
  lme4::lmer(y ~ x + (1 | subject), data = df_empty, REML = TRUE)
)
cat("type   :", r_lmer$type, "\n")
cat("class  :", r_lmer$class, "\n")
cat("message:", r_lmer$message, "\n\n")

# ============================================================
# 2. mixeff::lmm  — LMM
# ============================================================
cat("--- mixeff::lmm (LMM, empty data) ---\n")
r_lmm <- capture_msg(
  mixeff::lmm(y ~ x + (1 | subject), data = df_empty,
              control = mm_control(verbose = -1L))
)
cat("type   :", r_lmm$type, "\n")
cat("class  :", r_lmm$class, "\n")
cat("message:", r_lmm$message, "\n\n")

# ============================================================
# 3. lme4::glmer  — GLMM (binomial)
# ============================================================
df_full_bin <- data.frame(
  y       = rbinom(40, 1, 0.5),
  x       = rep(0:3, 10),
  subject = factor(rep(seq_len(10), each = 4))
)
df_empty_bin <- df_full_bin[0L, ]

cat("--- lme4::glmer (GLMM binomial, empty data) ---\n")
r_glmer <- capture_msg(
  lme4::glmer(y ~ x + (1 | subject), data = df_empty_bin, family = binomial())
)
cat("type   :", r_glmer$type, "\n")
cat("class  :", r_glmer$class, "\n")
cat("message:", r_glmer$message, "\n\n")

# ============================================================
# 4. mixeff::glmm  — GLMM (binomial)
# ============================================================
cat("--- mixeff::glmm (GLMM binomial, empty data) ---\n")
r_glmm <- capture_msg(
  mixeff::glmm(y ~ x + (1 | subject), data = df_empty_bin,
               family = binomial(),
               control = mm_control(verbose = -1L))
)
cat("type   :", r_glmm$type, "\n")
cat("class  :", r_glmm$class, "\n")
cat("message:", r_glmm$message, "\n\n")

# ============================================================
# Summary
# ============================================================
cat("=== SUMMARY ===\n")
cat(sprintf("lmer  : [%s] %s\n", r_lmer$type,  r_lmer$message))
cat(sprintf("lmm   : [%s] %s\n", r_lmm$type,   r_lmm$message))
cat(sprintf("glmer : [%s] %s\n", r_glmer$type,  r_glmer$message))
cat(sprintf("glmm  : [%s] %s\n", r_glmm$type,   r_glmm$message))
