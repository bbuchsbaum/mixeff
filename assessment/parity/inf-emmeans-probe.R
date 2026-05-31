## Empirical parity probe: inf-emmeans
## Dataset: cake
## Formula: angle ~ recipe*temperature + (1|recipe:replicate)
## Focus: emmeans + pairs() vs lme4-emmeans
## Quantities: estimates, SE, df, p (for marginal means and pairwise contrasts)

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(mixeff)
})

cat("=== SESSION INFO ===\n")
cat("lme4 version:    ", as.character(packageVersion("lme4")), "\n")
cat("lmerTest version:", as.character(packageVersion("lmerTest")), "\n")
cat("emmeans version: ", as.character(packageVersion("emmeans")), "\n")
cat("mixeff version:  ", as.character(packageVersion("mixeff")), "\n\n")

## ── 1. Data ─────────────────────────────────────────────────────────────────
data(cake, package = "lme4")
cat("=== DATASET ===\n")
cat("nrow:", nrow(cake), " ncol:", ncol(cake), "\n")
cat("recipe levels:", paste(levels(cake$recipe), collapse = ", "), "\n")
cat("temperature levels:", paste(levels(cake$temperature), collapse = ", "), "\n\n")

## ── 2. Fit lme4/lmerTest ────────────────────────────────────────────────────
cat("=== lme4/lmerTest FIT ===\n")
t_lme4 <- system.time({
  fit_lme4 <- lmer(angle ~ recipe * temperature + (1 | recipe:replicate),
                   data = cake, REML = TRUE)
})
cat("lme4 wall-clock (seconds):", t_lme4["elapsed"], "\n\n")

## ── 3. Fit mixeff ────────────────────────────────────────────────────────────
cat("=== mixeff FIT ===\n")
t_mixeff <- system.time({
  fit_mm <- tryCatch(
    lmm(angle ~ recipe * temperature + (1 | recipe:replicate),
        data = cake, REML = TRUE,
        control = mm_control(verbose = -1L)),
    error = function(e) e
  )
})
cat("mixeff wall-clock (seconds):", t_mixeff["elapsed"], "\n\n")

if (inherits(fit_mm, "condition")) {
  cat("!!! mixeff FIT ERROR !!!\n")
  cat("class:", paste(class(fit_mm), collapse = ", "), "\n")
  cat("message:", conditionMessage(fit_mm), "\n")
  print(fit_mm)
  cat("\n=== OUTCOME: mixeff-error (fit failed) ===\n")
  quit(status = 0)
}
cat("mixeff fit_status:", fit_mm$fit_status, "\n\n")

## ── 4. emmeans on lme4 ───────────────────────────────────────────────────────
cat("=== emmeans on lme4 (recipe marginal means) ===\n")
emm_lme4_recipe <- tryCatch(
  emmeans(fit_lme4, ~ recipe),
  error = function(e) { cat("lme4 emmeans error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(emm_lme4_recipe)) {
  cat("-- emmeans(lme4, ~ recipe) --\n")
  print(summary(emm_lme4_recipe))
}

cat("\n=== emmeans on lme4 (temperature marginal means) ===\n")
emm_lme4_temp <- tryCatch(
  emmeans(fit_lme4, ~ temperature),
  error = function(e) { cat("lme4 emmeans error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(emm_lme4_temp)) {
  cat("-- emmeans(lme4, ~ temperature) --\n")
  print(summary(emm_lme4_temp))
}

cat("\n=== emmeans on lme4 (recipe * temperature interaction) ===\n")
emm_lme4_int <- tryCatch(
  emmeans(fit_lme4, ~ recipe * temperature),
  error = function(e) { cat("lme4 emmeans error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(emm_lme4_int)) {
  cat("-- emmeans(lme4, ~ recipe * temperature) --\n")
  print(summary(emm_lme4_int))
}

cat("\n=== pairs() on lme4 (recipe contrasts) ===\n")
pairs_lme4_recipe <- tryCatch(
  pairs(emm_lme4_recipe),
  error = function(e) { cat("lme4 pairs error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(pairs_lme4_recipe)) {
  cat("-- pairs(emmeans(lme4, ~ recipe)) --\n")
  print(summary(pairs_lme4_recipe))
}

cat("\n=== pairs() on lme4 (temperature contrasts) ===\n")
pairs_lme4_temp <- tryCatch(
  pairs(emm_lme4_temp),
  error = function(e) { cat("lme4 pairs error:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(pairs_lme4_temp)) {
  cat("-- pairs(emmeans(lme4, ~ temperature)) --\n")
  print(summary(pairs_lme4_temp))
}

## ── 5. emmeans on mixeff ──────────────────────────────────────────────────────
cat("\n=== emmeans on mixeff (recipe marginal means) ===\n")
emm_mm_recipe <- tryCatch(
  emmeans(fit_mm, ~ recipe),
  error = function(e) {
    cat("!!! mixeff emmeans error !!!\n")
    cat("class:", paste(class(e), collapse = ", "), "\n")
    cat("message:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(emm_mm_recipe)) {
  cat("-- emmeans(mixeff, ~ recipe) --\n")
  print(summary(emm_mm_recipe))
}

cat("\n=== emmeans on mixeff (temperature marginal means) ===\n")
emm_mm_temp <- tryCatch(
  emmeans(fit_mm, ~ temperature),
  error = function(e) {
    cat("!!! mixeff emmeans error !!!\n")
    cat("class:", paste(class(e), collapse = ", "), "\n")
    cat("message:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(emm_mm_temp)) {
  cat("-- emmeans(mixeff, ~ temperature) --\n")
  print(summary(emm_mm_temp))
}

cat("\n=== emmeans on mixeff (recipe * temperature interaction) ===\n")
emm_mm_int <- tryCatch(
  emmeans(fit_mm, ~ recipe * temperature),
  error = function(e) {
    cat("!!! mixeff emmeans error !!!\n")
    cat("class:", paste(class(e), collapse = ", "), "\n")
    cat("message:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(emm_mm_int)) {
  cat("-- emmeans(mixeff, ~ recipe * temperature) --\n")
  print(summary(emm_mm_int))
}

cat("\n=== pairs() on mixeff (recipe contrasts) ===\n")
pairs_mm_recipe <- tryCatch(
  if (!is.null(emm_mm_recipe)) pairs(emm_mm_recipe) else stop("no emmeans object"),
  error = function(e) {
    cat("!!! mixeff pairs error !!!\n")
    cat("class:", paste(class(e), collapse = ", "), "\n")
    cat("message:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(pairs_mm_recipe)) {
  cat("-- pairs(emmeans(mixeff, ~ recipe)) --\n")
  print(summary(pairs_mm_recipe))
}

cat("\n=== pairs() on mixeff (temperature contrasts) ===\n")
pairs_mm_temp <- tryCatch(
  if (!is.null(emm_mm_temp)) pairs(emm_mm_temp) else stop("no emmeans object"),
  error = function(e) {
    cat("!!! mixeff pairs error !!!\n")
    cat("class:", paste(class(e), collapse = ", "), "\n")
    cat("message:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(pairs_mm_temp)) {
  cat("-- pairs(emmeans(mixeff, ~ temperature)) --\n")
  print(summary(pairs_mm_temp))
}

## ── 6. Numerical comparison ───────────────────────────────────────────────────
cat("\n=== NUMERICAL COMPARISON ===\n")
tols <- list(estimate = 1e-4, se = 1e-4, df = 1.0, p = 1e-4)

compare_vec <- function(label, lme4_val, mm_val, tol, fmt = "%.6f") {
  lme4_val <- as.numeric(lme4_val)
  mm_val   <- as.numeric(mm_val)
  if (is.null(mm_val) || all(is.na(mm_val))) {
    cat(sprintf("%-50s  lme4=[%s]  mm=MISSING\n",
                label,
                paste(sprintf(fmt, lme4_val), collapse = ", ")))
    return(invisible(NA_real_))
  }
  if (length(lme4_val) != length(mm_val)) {
    cat(sprintf("%-50s  LENGTH MISMATCH: lme4=%d mm=%d\n",
                label, length(lme4_val), length(mm_val)))
    return(invisible(NA_real_))
  }
  diff   <- max(abs(lme4_val - mm_val), na.rm = TRUE)
  status <- if (diff <= tol) "WITHIN-TOL" else "EXCEEDS-TOL"
  cat(sprintf("%-50s  maxAbsDiff=%.4e  tol=%.0e  [%s]\n",
              label, diff, tol, status))
  invisible(diff)
}

## Helper: extract summary columns from an emmGrid summary
extract_emm <- function(emm_obj, col) {
  if (is.null(emm_obj)) return(NULL)
  s <- as.data.frame(summary(emm_obj))
  if (col %in% names(s)) s[[col]] else NULL
}

## ── recipe marginal means ──
cat("\n-- recipe marginal means --\n")
compare_vec("recipe emmean: estimate",
            extract_emm(emm_lme4_recipe, "emmean"),
            extract_emm(emm_mm_recipe,   "emmean"),
            tols$estimate)
compare_vec("recipe emmean: SE",
            extract_emm(emm_lme4_recipe, "SE"),
            extract_emm(emm_mm_recipe,   "SE"),
            tols$se)
compare_vec("recipe emmean: df",
            extract_emm(emm_lme4_recipe, "df"),
            extract_emm(emm_mm_recipe,   "df"),
            tols$df)

## ── temperature marginal means ──
cat("\n-- temperature marginal means --\n")
compare_vec("temperature emmean: estimate",
            extract_emm(emm_lme4_temp, "emmean"),
            extract_emm(emm_mm_temp,   "emmean"),
            tols$estimate)
compare_vec("temperature emmean: SE",
            extract_emm(emm_lme4_temp, "SE"),
            extract_emm(emm_mm_temp,   "SE"),
            tols$se)
compare_vec("temperature emmean: df",
            extract_emm(emm_lme4_temp, "df"),
            extract_emm(emm_mm_temp,   "df"),
            tols$df)

## ── recipe * temperature cell means ──
cat("\n-- recipe * temperature cell means --\n")
compare_vec("recipe*temp emmean: estimate",
            extract_emm(emm_lme4_int, "emmean"),
            extract_emm(emm_mm_int,   "emmean"),
            tols$estimate)
compare_vec("recipe*temp emmean: SE",
            extract_emm(emm_lme4_int, "SE"),
            extract_emm(emm_mm_int,   "SE"),
            tols$se)
compare_vec("recipe*temp emmean: df",
            extract_emm(emm_lme4_int, "df"),
            extract_emm(emm_mm_int,   "df"),
            tols$df)

## ── recipe pairwise contrasts ──
cat("\n-- recipe pairwise contrasts --\n")
compare_vec("pairs(recipe): estimate",
            extract_emm(pairs_lme4_recipe, "estimate"),
            extract_emm(pairs_mm_recipe,   "estimate"),
            tols$estimate)
compare_vec("pairs(recipe): SE",
            extract_emm(pairs_lme4_recipe, "SE"),
            extract_emm(pairs_mm_recipe,   "SE"),
            tols$se)
compare_vec("pairs(recipe): df",
            extract_emm(pairs_lme4_recipe, "df"),
            extract_emm(pairs_mm_recipe,   "df"),
            tols$df)
compare_vec("pairs(recipe): p.value",
            extract_emm(pairs_lme4_recipe, "p.value"),
            extract_emm(pairs_mm_recipe,   "p.value"),
            tols$p)

## ── temperature pairwise contrasts ──
cat("\n-- temperature pairwise contrasts --\n")
compare_vec("pairs(temp): estimate",
            extract_emm(pairs_lme4_temp, "estimate"),
            extract_emm(pairs_mm_temp,   "estimate"),
            tols$estimate)
compare_vec("pairs(temp): SE",
            extract_emm(pairs_lme4_temp, "SE"),
            extract_emm(pairs_mm_temp,   "SE"),
            tols$se)
compare_vec("pairs(temp): df",
            extract_emm(pairs_lme4_temp, "df"),
            extract_emm(pairs_mm_temp,   "df"),
            tols$df)
compare_vec("pairs(temp): p.value",
            extract_emm(pairs_lme4_temp, "p.value"),
            extract_emm(pairs_mm_temp,   "p.value"),
            tols$p)

## ── 7. mm_means / mm_comparisons native API ─────────────────────────────────
cat("\n=== MIXEFF NATIVE API: mm_means / mm_comparisons ===\n")

cat("-- mm_means(fit_mm, ~ recipe) --\n")
mm_means_recipe <- tryCatch(
  mm_means(fit_mm, ~ recipe),
  error = function(e) {
    cat("mm_means error:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(mm_means_recipe)) print(mm_means_recipe)

cat("\n-- mm_means(fit_mm, ~ temperature) --\n")
mm_means_temp <- tryCatch(
  mm_means(fit_mm, ~ temperature),
  error = function(e) {
    cat("mm_means error:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(mm_means_temp)) print(mm_means_temp)

cat("\n-- mm_comparisons(fit_mm, ~ recipe) --\n")
mm_comp_recipe <- tryCatch(
  mm_comparisons(fit_mm, ~ recipe),
  error = function(e) {
    cat("mm_comparisons error:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(mm_comp_recipe)) print(mm_comp_recipe)

cat("\n-- mm_comparisons(fit_mm, ~ temperature) --\n")
mm_comp_temp <- tryCatch(
  mm_comparisons(fit_mm, ~ temperature),
  error = function(e) {
    cat("mm_comparisons error:", conditionMessage(e), "\n")
    NULL
  }
)
if (!is.null(mm_comp_temp)) print(mm_comp_temp)

## ── 8. Speed comparison ──────────────────────────────────────────────────────
cat("\n=== SPEED COMPARISON ===\n")
NREPS <- 5
t_lme4_emm <- system.time(
  for (i in seq_len(NREPS))
    emmeans(fit_lme4, ~ recipe)
)
t_mm_emm <- system.time(
  for (i in seq_len(NREPS)) {
    tryCatch(emmeans(fit_mm, ~ recipe), error = function(e) NULL)
  }
)
cat(sprintf("emmeans(lme4, ~recipe) mean/call: %.4f s\n",
            t_lme4_emm["elapsed"] / NREPS))
cat(sprintf("emmeans(mixeff, ~recipe) mean/call: %.4f s\n",
            t_mm_emm["elapsed"] / NREPS))
if (t_lme4_emm["elapsed"] > 0) {
  cat(sprintf("ratio (mm/lme4): %.2fx\n",
              t_mm_emm["elapsed"] / t_lme4_emm["elapsed"]))
}

cat("\n=== PROBE COMPLETE ===\n")
