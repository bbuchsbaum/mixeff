library(lme4)
library(mixeff)

# Build a small dataset where the response is a character vector
set.seed(42)
n <- 60
df <- data.frame(
  y_char = sample(c("low", "mid", "high"), n, replace = TRUE),  # character response
  x      = rnorm(n),
  subj   = factor(rep(1:10, each = 6))
)

cat("===== lme4::lmer with character response =====\n")
lme4_lmer_msg <- tryCatch(
  lme4::lmer(y_char ~ x + (1 | subj), data = df),
  error   = function(e) paste("ERROR:", conditionMessage(e)),
  warning = function(w) paste("WARNING:", conditionMessage(w))
)
cat(lme4_lmer_msg, "\n\n")

cat("===== lme4::glmer with character response =====\n")
lme4_glmer_msg <- tryCatch(
  lme4::glmer(y_char ~ x + (1 | subj), data = df, family = binomial),
  error   = function(e) paste("ERROR:", conditionMessage(e)),
  warning = function(w) paste("WARNING:", conditionMessage(w))
)
cat(lme4_glmer_msg, "\n\n")

cat("===== mixeff::lmm with character response =====\n")
mixeff_lmm_msg <- tryCatch(
  mixeff::lmm(y_char ~ x + (1 | subj), data = df,
              control = mm_control(verbose = -1L)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("ERROR [class: ", cls, "]: ", conditionMessage(e))
  },
  warning = function(w) paste("WARNING:", conditionMessage(w))
)
cat(mixeff_lmm_msg, "\n\n")

cat("===== mixeff::glmm with character response =====\n")
mixeff_glmm_msg <- tryCatch(
  mixeff::glmm(y_char ~ x + (1 | subj), data = df, family = "binomial",
               control = mm_control(verbose = -1L)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("ERROR [class: ", cls, "]: ", conditionMessage(e))
  },
  warning = function(w) paste("WARNING:", conditionMessage(w))
)
cat(mixeff_glmm_msg, "\n\n")

# Also test with binary character (yes/no) for glmm
cat("===== lme4::glmer binary character response =====\n")
df$y_bin <- sample(c("yes", "no"), n, replace = TRUE)
lme4_glmer_bin <- tryCatch(
  lme4::glmer(y_bin ~ x + (1 | subj), data = df, family = binomial),
  error   = function(e) paste("ERROR:", conditionMessage(e)),
  warning = function(w) paste("WARNING:", conditionMessage(w))
)
cat(lme4_glmer_bin, "\n\n")

cat("===== mixeff::glmm binary character response =====\n")
mixeff_glmm_bin <- tryCatch(
  mixeff::glmm(y_bin ~ x + (1 | subj), data = df, family = "binomial",
               control = mm_control(verbose = -1L)),
  error   = function(e) {
    cls <- paste(class(e), collapse = ", ")
    paste0("ERROR [class: ", cls, "]: ", conditionMessage(e))
  },
  warning = function(w) paste("WARNING:", conditionMessage(w))
)
cat(mixeff_glmm_bin, "\n\n")

cat("Done.\n")
