suppressPackageStartupMessages({ library(lme4); library(mixeff) })
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper to extract diagnostics as text
show_diag <- function(fit, label) {
  d <- tryCatch(diagnostics(fit), error=function(e) NULL)
  if (is.null(d)) { cat(sprintf("  [%s] diagnostics: NULL\n", label)); return() }
  cat(sprintf("  [%s] diagnostics class: %s\n", label, paste(class(d), collapse=", ")))
  # Try coercing to data.frame
  df_d <- tryCatch(as.data.frame(d), error=function(e) NULL)
  if (!is.null(df_d) && nrow(df_d) > 0) {
    cols <- intersect(c("code","severity","stage","message"), names(df_d))
    print(df_d[, cols], row.names=FALSE)
  } else if (is.list(d)) {
    # Try printing each entry
    for (item in d) {
      if (is.list(item) && !is.null(item$code))
        cat(sprintf("    %s [%s/%s]: %s\n", item$code, item$severity %||% "?", item$stage %||% "?", item$message %||% ""))
    }
  } else {
    print(d)
  }
}

# ---- Scenario 1: extreme scale mismatch ----
cat("=== SCENARIO 1: extreme scale mismatch ===\n")
set.seed(2024)
n_subj <- 15; n_obs <- 4
df <- data.frame(
  subject = factor(rep(seq_len(n_subj), each = n_obs)),
  x_small = rnorm(n_subj * n_obs, 0, 1),
  x_large = rnorm(n_subj * n_obs, 1e6, 1e3),
  y       = rnorm(n_subj * n_obs)
)

# lme4 tight budget
cat("\n--- lme4::lmer tight budget (maxfun=10) ---\n")
l4m1 <- character(); l4w1 <- character()
lf1 <- tryCatch(
  withCallingHandlers(
    lmer(y ~ x_small + x_large + (1 + x_large | subject), data=df,
         control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10),calc.derivs=TRUE)),
    message=function(m){l4m1<<-c(l4m1,conditionMessage(m));invokeRestart("muffleMessage")},
    warning=function(w){l4w1<<-c(l4w1,conditionMessage(w));invokeRestart("muffleWarning")}
  ),
  error=function(e) structure(list(msg=conditionMessage(e)),class="err_h")
)
if(inherits(lf1,"err_h")) cat(sprintf("  ERROR: %s\n",lf1$msg)) else cat(sprintf("  class: %s\n",class(lf1)[1]))
for(m in l4m1) cat(sprintf("  MSG : %s\n",trimws(m)))
for(w in l4w1) cat(sprintf("  WARN: %s\n",trimws(w)))

# lme4 default budget
cat("\n--- lme4::lmer DEFAULT budget ---\n")
l4m2 <- character(); l4w2 <- character()
lf2 <- tryCatch(
  withCallingHandlers(
    lmer(y ~ x_small + x_large + (1 + x_large | subject), data=df),
    message=function(m){l4m2<<-c(l4m2,conditionMessage(m));invokeRestart("muffleMessage")},
    warning=function(w){l4w2<<-c(l4w2,conditionMessage(w));invokeRestart("muffleWarning")}
  ),
  error=function(e) structure(list(msg=conditionMessage(e)),class="err_h")
)
if(inherits(lf2,"err_h")) cat(sprintf("  ERROR: %s\n",lf2$msg)) else cat(sprintf("  class: %s\n",class(lf2)[1]))
for(m in l4m2) cat(sprintf("  MSG : %s\n",trimws(m)))
for(w in l4w2) cat(sprintf("  WARN: %s\n",trimws(w)))

# mixeff lmm
cat("\n--- mixeff::lmm DEFAULT budget ---\n")
mm1w <- character()
mf1 <- tryCatch(
  withCallingHandlers(
    lmm(y ~ x_small + x_large + (1 + x_large | subject), data=df,
        control=mm_control(verbose=-1L)),
    message=function(m){cat(sprintf("  MSG: %s\n",trimws(conditionMessage(m))));invokeRestart("muffleMessage")},
    warning=function(w){mm1w<<-c(mm1w,conditionMessage(w));invokeRestart("muffleWarning")}
  ),
  error=function(e) structure(list(msg=conditionMessage(e),cls=class(e)),class="err_h")
)
for(w in mm1w) cat(sprintf("  WARN: %s\n",trimws(w)))
if(inherits(mf1,"err_h")) {
  cat(sprintf("  ERROR (%s): %s\n",paste(mf1$cls,collapse=","),mf1$msg))
} else {
  cat(sprintf("  fit_status: %s\n",fit_status(mf1)))
  # Show artifact-level diagnostics
  cat("  -- artifact$diagnostics --\n")
  for(item in mf1$artifact$diagnostics %||% list())
    cat(sprintf("    %s [%s/%s]: %s\n", item$code, item$severity, item$stage, item$message))
  # Show certificate diagnostics
  cat("  -- optimizer_certificate$diagnostics --\n")
  for(item in mf1$artifact$optimizer_certificate$diagnostics %||% list())
    cat(sprintf("    %s [%s/%s]: %s\n", item$code, item$severity, item$stage, item$message))
  show_diag(mf1, "lmm s1")
}

# ---- Scenario 2: near-zero within-group variance ----
cat("\n\n=== SCENARIO 2: near-zero within-group variance ===\n")
set.seed(999)
n2 <- 10; n_per2 <- 3
df2 <- data.frame(
  subject = factor(rep(seq_len(n2), each=n_per2)),
  x       = rep(rnorm(n2, sd=10), each=n_per2) + rnorm(n2*n_per2, sd=0.001),
  y       = rnorm(n2*n_per2)
)

cat("\n--- lme4::lmer s2 ---\n")
l4m3 <- character(); l4w3 <- character()
lf3 <- tryCatch(
  withCallingHandlers(
    lmer(y ~ x + (1 + x | subject), data=df2,
         control=lmerControl(optimizer="bobyqa")),
    message=function(m){l4m3<<-c(l4m3,conditionMessage(m));invokeRestart("muffleMessage")},
    warning=function(w){l4w3<<-c(l4w3,conditionMessage(w));invokeRestart("muffleWarning")}
  ),
  error=function(e) structure(list(msg=conditionMessage(e)),class="err_h")
)
if(inherits(lf3,"err_h")) cat(sprintf("  ERROR: %s\n",lf3$msg)) else cat(sprintf("  class: %s\n",class(lf3)[1]))
for(m in l4m3) cat(sprintf("  MSG : %s\n",trimws(m)))
for(w in l4w3) cat(sprintf("  WARN: %s\n",trimws(w)))

cat("\n--- mixeff::lmm s2 ---\n")
mm2w <- character()
mf2 <- tryCatch(
  withCallingHandlers(
    lmm(y ~ x + (1 + x | subject), data=df2,
        control=mm_control(verbose=-1L)),
    message=function(m){cat(sprintf("  MSG: %s\n",trimws(conditionMessage(m))));invokeRestart("muffleMessage")},
    warning=function(w){mm2w<<-c(mm2w,conditionMessage(w));invokeRestart("muffleWarning")}
  ),
  error=function(e) structure(list(msg=conditionMessage(e),cls=class(e)),class="err_h")
)
for(w in mm2w) cat(sprintf("  WARN: %s\n",trimws(w)))
if(inherits(mf2,"err_h")) {
  cat(sprintf("  ERROR (%s): %s\n",paste(mf2$cls,collapse=","),mf2$msg))
} else {
  cat(sprintf("  fit_status: %s\n",fit_status(mf2)))
  cat("  -- artifact$diagnostics --\n")
  for(item in mf2$artifact$diagnostics %||% list())
    cat(sprintf("    %s [%s/%s]: %s\n", item$code, item$severity, item$stage, item$message))
  cat("  -- optimizer_certificate$diagnostics --\n")
  for(item in mf2$artifact$optimizer_certificate$diagnostics %||% list())
    cat(sprintf("    %s [%s/%s]: %s\n", item$code, item$severity, item$stage, item$message))
  show_diag(mf2, "lmm s2")
}

# ---- Scenario 3: glmm convergence hard ----
cat("\n\n=== SCENARIO 3: glmm with extreme scale mismatch ===\n")
set.seed(2025)
n3 <- 20; n_per3 <- 5
df3 <- data.frame(
  subject = factor(rep(seq_len(n3), each=n_per3)),
  x_small = rnorm(n3*n_per3, 0, 1),
  x_large = rnorm(n3*n_per3, 1e6, 1e3),
  y       = rbinom(n3*n_per3, 1, 0.5)
)

cat("\n--- lme4::glmer s3 (tight budget) ---\n")
l4m4 <- character(); l4w4 <- character()
lf4 <- tryCatch(
  withCallingHandlers(
    glmer(y ~ x_small + x_large + (1 + x_large | subject), data=df3,
          family=binomial(),
          control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=20))),
    message=function(m){l4m4<<-c(l4m4,conditionMessage(m));invokeRestart("muffleMessage")},
    warning=function(w){l4w4<<-c(l4w4,conditionMessage(w));invokeRestart("muffleWarning")}
  ),
  error=function(e) structure(list(msg=conditionMessage(e)),class="err_h")
)
if(inherits(lf4,"err_h")) cat(sprintf("  ERROR: %s\n",lf4$msg)) else cat(sprintf("  class: %s\n",class(lf4)[1]))
for(m in l4m4) cat(sprintf("  MSG : %s\n",trimws(m)))
for(w in l4w4) cat(sprintf("  WARN: %s\n",trimws(w)))

cat("\n--- mixeff::glmm s3 ---\n")
mm3w <- character()
mf3 <- tryCatch(
  withCallingHandlers(
    glmm(y ~ x_small + x_large + (1 + x_large | subject), data=df3,
         family=binomial(), control=mm_control(verbose=-1L)),
    message=function(m){cat(sprintf("  MSG: %s\n",trimws(conditionMessage(m))));invokeRestart("muffleMessage")},
    warning=function(w){mm3w<<-c(mm3w,conditionMessage(w));invokeRestart("muffleWarning")}
  ),
  error=function(e) structure(list(msg=conditionMessage(e),cls=class(e)),class="err_h")
)
for(w in mm3w) cat(sprintf("  WARN: %s\n",trimws(w)))
if(inherits(mf3,"err_h")) {
  cat(sprintf("  ERROR (%s): %s\n",paste(mf3$cls,collapse=","),mf3$msg))
} else {
  cat(sprintf("  fit_status: %s\n",fit_status(mf3)))
  cat("  -- artifact$diagnostics --\n")
  for(item in mf3$artifact$diagnostics %||% list())
    cat(sprintf("    %s [%s/%s]: %s\n", item$code, item$severity, item$stage, item$message))
  cat("  -- optimizer_certificate$diagnostics --\n")
  for(item in mf3$artifact$optimizer_certificate$diagnostics %||% list())
    cat(sprintf("    %s [%s/%s]: %s\n", item$code, item$severity, item$stage, item$message))
  show_diag(mf3, "glmm s3")
}

cat("\n=== DONE ===\n")
