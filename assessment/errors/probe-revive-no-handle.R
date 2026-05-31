## probe-revive-no-handle.R
##
## Scenario: use a revived/serialized fit after the Rust handle is dropped.
## Tests both lme4/lmerTest and mixeff lmm()/glmm().
## Captures exact error/warning messages via tryCatch.

suppressPackageStartupMessages({
  library(lme4)
  library(mixeff)
})

## Explicit namespace aliases to avoid masking confusion
lme4_fixef  <- lme4::fixef
lme4_ranef  <- lme4::ranef
mm_fixef    <- mixeff::fixef
mm_ranef    <- mixeff::ranef

## ============================================================
## Shared data: simple random-intercept design
## ============================================================
set.seed(42L)
n_subjects <- 10L
n_per      <- 5L
subject    <- factor(rep(seq_len(n_subjects), each = n_per))
x          <- rep(seq_len(n_per) - 1L, n_subjects)
b0         <- rnorm(n_subjects, sd = 0.5)
y          <- 2 + 0.4 * x + b0[as.integer(subject)] + rnorm(length(x), sd = 0.3)
counts     <- rpois(length(x), lambda = exp(0.5 + 0.2 * x + b0[as.integer(subject)]))
df         <- data.frame(y = y, counts = counts, x = x, subject = subject)

## ============================================================
## Helper: save to RDS, reload, call revive() if applicable
## ============================================================
save_reload <- function(fit, use_revive = FALSE) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(fit, tf)
  out <- readRDS(tf)
  if (use_revive) out <- revive(out)
  out
}

## ============================================================
## lme4 / lmerTest — LMM
## ============================================================
cat("\n====== lme4 lmer: fit, saveRDS/readRDS, then post-fit verbs ======\n")

lme4_fit <- lmer(y ~ x + (1 | subject), data = df, REML = TRUE)
lme4_restored <- save_reload(lme4_fit, use_revive = FALSE)

## lme4 has no explicit revive — objects just survive serialization.
## Demonstrate that lme4 works without any special revival call:
cat("lme4 fixef after reload: ")
r_lme4_fixef <- tryCatch(lme4_fixef(lme4_restored), error = function(e) e)
if (inherits(r_lme4_fixef, "error")) {
  cat("ERROR:", conditionMessage(r_lme4_fixef), "\n")
} else {
  cat("OK —", paste(names(r_lme4_fixef), round(r_lme4_fixef, 4), sep = "=", collapse = ", "), "\n")
}

cat("lme4 predict after reload: ")
r_lme4_pred <- tryCatch(head(predict(lme4_restored), 3), error = function(e) e)
if (inherits(r_lme4_pred, "error")) {
  cat("ERROR:", conditionMessage(r_lme4_pred), "\n")
} else {
  cat("OK — first 3:", paste(round(r_lme4_pred, 4), collapse = ", "), "\n")
}

cat("lme4 ranef after reload: ")
r_lme4_ranef <- tryCatch(lme4_ranef(lme4_restored), error = function(e) e)
if (inherits(r_lme4_ranef, "error")) {
  cat("ERROR:", conditionMessage(r_lme4_ranef), "\n")
} else {
  cat("OK — n groups:", nrow(r_lme4_ranef$subject), "\n")
}

## ============================================================
## lme4 / lmerTest — GLMM (Poisson)
## ============================================================
cat("\n====== lme4 glmer (Poisson): fit, saveRDS/readRDS, then post-fit verbs ======\n")

lme4_glmm <- glmer(counts ~ x + (1 | subject), data = df, family = poisson())
lme4_glmm_restored <- save_reload(lme4_glmm, use_revive = FALSE)

cat("lme4 glmer fixef after reload: ")
r_g_fixef <- tryCatch(lme4_fixef(lme4_glmm_restored), error = function(e) e)
if (inherits(r_g_fixef, "error")) {
  cat("ERROR:", conditionMessage(r_g_fixef), "\n")
} else {
  cat("OK —", paste(names(r_g_fixef), round(r_g_fixef, 4), sep = "=", collapse = ", "), "\n")
}

cat("lme4 glmer predict after reload: ")
r_g_pred <- tryCatch(head(predict(lme4_glmm_restored, type = "response"), 3), error = function(e) e)
if (inherits(r_g_pred, "error")) {
  cat("ERROR:", conditionMessage(r_g_pred), "\n")
} else {
  cat("OK — first 3:", paste(round(r_g_pred, 4), collapse = ", "), "\n")
}

## ============================================================
## mixeff lmm — scenario 1: readRDS WITHOUT revive(), then use
## ============================================================
cat("\n====== mixeff lmm: readRDS WITHOUT revive() ======\n")

mm_fit <- lmm(y ~ x + (1 | subject), data = df, control = mm_control(verbose = -1))
mm_raw_restored <- save_reload(mm_fit, use_revive = FALSE)  # no revive

cat("mixeff fixef (no revive): ")
r_mm_fixef_norev <- tryCatch(fixef(mm_raw_restored), error = function(e) e)
if (inherits(r_mm_fixef_norev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_fixef_norev), "\n")
} else {
  cat("OK —", paste(names(r_mm_fixef_norev), round(r_mm_fixef_norev, 4), sep = "=", collapse = ", "), "\n")
}

cat("mixeff predict (no revive): ")
r_mm_pred_norev <- tryCatch(head(predict(mm_raw_restored), 3), error = function(e) e)
if (inherits(r_mm_pred_norev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_pred_norev), "\n")
} else {
  cat("OK — first 3:", paste(round(r_mm_pred_norev, 4), collapse = ", "), "\n")
}

cat("mixeff ranef (no revive): ")
r_mm_ranef_norev <- tryCatch(ranef(mm_raw_restored), error = function(e) e)
if (inherits(r_mm_ranef_norev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_ranef_norev), "\n")
} else {
  cat("OK — n groups:", nrow(r_mm_ranef_norev$subject), "\n")
}

cat("mixeff summary (no revive): ")
r_mm_summ_norev <- tryCatch(
  capture.output(summary(mm_raw_restored)),
  error = function(e) e,
  warning = function(w) w
)
if (inherits(r_mm_summ_norev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_summ_norev), "\n")
} else if (inherits(r_mm_summ_norev, "warning")) {
  cat("WARNING:", conditionMessage(r_mm_summ_norev), "\n")
} else {
  cat("OK\n")
}

cat("mixeff contrast (no revive): ")
L <- c("(Intercept)" = 0, x = 1)
r_mm_ct_norev <- tryCatch(contrast(mm_raw_restored, L), error = function(e) e)
if (inherits(r_mm_ct_norev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_ct_norev), "\n")
} else {
  cat("OK — p_value:", r_mm_ct_norev$table$p_value, "\n")
}

## ============================================================
## mixeff lmm — scenario 2: readRDS WITH revive(), then use
## ============================================================
cat("\n====== mixeff lmm: readRDS WITH revive() ======\n")

mm_revived <- save_reload(mm_fit, use_revive = TRUE)

cat("fit_handle_alive after revive: ", fit_handle_alive(mm_revived), "\n")

cat("mixeff fixef (revived): ")
r_mm_fixef_rev <- tryCatch(fixef(mm_revived), error = function(e) e)
if (inherits(r_mm_fixef_rev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_fixef_rev), "\n")
} else {
  cat("OK —", paste(names(r_mm_fixef_rev), round(r_mm_fixef_rev, 4), sep = "=", collapse = ", "), "\n")
}

cat("mixeff predict (revived): ")
r_mm_pred_rev <- tryCatch(head(predict(mm_revived), 3), error = function(e) e)
if (inherits(r_mm_pred_rev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_pred_rev), "\n")
} else {
  cat("OK — first 3:", paste(round(r_mm_pred_rev, 4), collapse = ", "), "\n")
}

cat("mixeff ranef (revived): ")
r_mm_ranef_rev <- tryCatch(ranef(mm_revived), error = function(e) e)
if (inherits(r_mm_ranef_rev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_ranef_rev), "\n")
} else {
  cat("OK — n groups:", nrow(r_mm_ranef_rev$subject), "\n")
}

cat("mixeff contrast (revived): ")
r_mm_ct_rev <- tryCatch(contrast(mm_revived, L), error = function(e) e)
if (inherits(r_mm_ct_rev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_ct_rev), "\n")
} else {
  cat("OK — p_value:", r_mm_ct_rev$table$p_value, "\n")
}

cat("mixeff random_blocks (revived): ")
r_mm_rb_rev <- tryCatch(random_blocks(mm_revived), error = function(e) e)
if (inherits(r_mm_rb_rev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_rb_rev), "\n")
} else {
  cat("OK — groups:", r_mm_rb_rev$table$group, "\n")
}

cat("mixeff inference_table (revived): ")
r_mm_inf_rev <- tryCatch(inference_table(mm_revived), error = function(e) e)
if (inherits(r_mm_inf_rev, "error")) {
  cat("ERROR:", conditionMessage(r_mm_inf_rev), "\n")
} else {
  cat("OK — methods:", paste(r_mm_inf_rev$table$method, collapse = ", "), "\n")
}

## ============================================================
## mixeff lmm — scenario 3: deliberately break revive by
## nulling artifact_json and artifact, then call revive()
## ============================================================
cat("\n====== mixeff lmm: revive() on broken fit (no artifact) ======\n")

broken_fit <- mm_fit
broken_fit$artifact <- NULL
broken_fit$fit$artifact_json <- NULL
broken_fit$rust_handle <- NULL

cat("revive() on broken fit: ")
r_broken_revive <- tryCatch(revive(broken_fit), error = function(e) e)
if (inherits(r_broken_revive, "error")) {
  cat("ERROR class:", class(r_broken_revive), "\n")
  cat("ERROR message:", conditionMessage(r_broken_revive), "\n")
} else {
  cat("OK (unexpected)\n")
}

## ============================================================
## mixeff lmm — scenario 4: revive() on a non-mm_fit object
## ============================================================
cat("\n====== mixeff revive() on non-mm_fit object ======\n")

cat("revive() on plain list: ")
r_bad_revive <- tryCatch(revive(list(a = 1)), error = function(e) e)
if (inherits(r_bad_revive, "error")) {
  cat("ERROR class:", class(r_bad_revive), "\n")
  cat("ERROR message:", conditionMessage(r_bad_revive), "\n")
} else {
  cat("OK (unexpected)\n")
}

cat("revive() on integer: ")
r_bad_revive2 <- tryCatch(revive(42L), error = function(e) e)
if (inherits(r_bad_revive2, "error")) {
  cat("ERROR class:", class(r_bad_revive2), "\n")
  cat("ERROR message:", conditionMessage(r_bad_revive2), "\n")
} else {
  cat("OK (unexpected)\n")
}

## ============================================================
## mixeff glmm — revive scenario
## ============================================================
cat("\n====== mixeff glmm: readRDS WITH revive() ======\n")

mm_glmm_fit <- tryCatch(
  glmm(counts ~ x + (1 | subject), data = df, family = poisson(),
       control = mm_control(verbose = -1)),
  error = function(e) e
)

if (inherits(mm_glmm_fit, "error")) {
  cat("glmm fit failed:", conditionMessage(mm_glmm_fit), "\n")
} else {
  mm_glmm_revived <- save_reload(mm_glmm_fit, use_revive = TRUE)
  cat("fit_handle_alive after glmm revive: ", fit_handle_alive(mm_glmm_revived), "\n")

  cat("mixeff glmm fixef (revived): ")
  r_gg_fixef <- tryCatch(fixef(mm_glmm_revived), error = function(e) e)
  if (inherits(r_gg_fixef, "error")) {
    cat("ERROR:", conditionMessage(r_gg_fixef), "\n")
  } else {
    cat("OK —", paste(names(r_gg_fixef), round(r_gg_fixef, 4), sep = "=", collapse = ", "), "\n")
  }

  cat("mixeff glmm ranef (revived): ")
  r_gg_ranef <- tryCatch(ranef(mm_glmm_revived), error = function(e) e)
  if (inherits(r_gg_ranef, "error")) {
    cat("ERROR:", conditionMessage(r_gg_ranef), "\n")
  } else {
    cat("OK — n groups:", nrow(r_gg_ranef$subject), "\n")
  }

  cat("mixeff glmm predict (revived): ")
  r_gg_pred <- tryCatch(head(predict(mm_glmm_revived, type = "response"), 3), error = function(e) e)
  if (inherits(r_gg_pred, "error")) {
    cat("ERROR:", conditionMessage(r_gg_pred), "\n")
  } else {
    cat("OK — first 3:", paste(round(r_gg_pred, 4), collapse = ", "), "\n")
  }
}

cat("\n====== probe complete ======\n")
