#!/usr/bin/env Rscript

# Run from the package root:
#   Rscript examples/sdamr-lme4-companion.R
#
# The examples mirror Chapter 9 of:
# https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("This example requires the lme4 package.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(mixeff)
})

run_slow <- identical(tolower(Sys.getenv("MIXEFF_EXAMPLE_RUN_SLOW")), "true")

load_sdamr_dataset <- function(name, fixture) {
  if (requireNamespace("sdamr", quietly = TRUE)) {
    env <- new.env(parent = emptyenv())
    utils::data(list = name, package = "sdamr", envir = env)
    if (exists(name, envir = env, inherits = FALSE)) {
      return(get(name, envir = env, inherits = FALSE))
    }
  }

  path <- file.path("tests", "fixtures", fixture)
  if (!file.exists(path)) {
    stop(
      sprintf("Install sdamr or run from the package root with `%s` present.", path),
      call. = FALSE
    )
  }
  utils::read.csv(path, stringsAsFactors = TRUE)
}

center <- function(x) {
  x - mean(x, na.rm = TRUE)
}

prepare_anchoring <- function() {
  dat <- load_sdamr_dataset("anchoring", "sdamr_anchoring.csv")
  dat$anchor <- factor(dat$anchor)
  contrasts(dat$anchor) <- c(1 / 2, -1 / 2)
  dat$referrer <- factor(dat$referrer)
  if (!"anchor_contrast" %in% names(dat)) {
    dat$anchor_contrast <- ifelse(dat$anchor == "high", 1 / 2, -1 / 2)
  }
  dat
}

prepare_speeddate <- function() {
  dat <- load_sdamr_dataset("speeddate", "sdamr_speeddate_lmm.csv")
  dat$iid <- factor(dat$iid)
  dat$pid <- factor(dat$pid)
  if (!"other_attr_c" %in% names(dat)) {
    dat$other_attr_c <- center(dat$other_attr)
  }
  if (!"other_intel_c" %in% names(dat)) {
    dat$other_intel_c <- center(dat$other_intel)
  }
  keep <- c("other_like", "other_attr_c", "other_intel_c", "iid", "pid")
  stats::na.omit(dat[keep])
}

comparison_key <- function(x) {
  gsub(": ", "", as.character(x), fixed = TRUE)
}

expected_fixed_effects <- function(lme4_fit, basis = NULL) {
  beta <- lme4::fixef(lme4_fit)
  if (identical(basis, "anchor_treatment") &&
      all(c("(Intercept)", "anchor1") %in% names(beta))) {
    beta <- c(
      "(Intercept)" = unname(beta[["(Intercept)"]] + 0.5 * beta[["anchor1"]]),
      "anchor: low" = unname(-beta[["anchor1"]])
    )
  }
  beta
}

fixed_effects <- function(mm_fit, lme4_fit, basis = NULL) {
  mm_beta <- fixef(mm_fit)
  lme4_beta <- expected_fixed_effects(lme4_fit, basis)
  names(mm_beta) <- comparison_key(names(mm_beta))
  names(lme4_beta) <- comparison_key(names(lme4_beta))
  terms <- union(names(mm_beta), names(lme4_beta))
  data.frame(
    term = terms,
    mixeff = unname(mm_beta[terms]),
    lme4 = unname(lme4_beta[terms]),
    difference = unname(mm_beta[terms] - lme4_beta[terms]),
    row.names = NULL,
    check.names = FALSE
  )
}

fit_stats <- function(mm_fit, lme4_fit) {
  data.frame(
    statistic = c("sigma", "logLik", "AIC", "BIC"),
    mixeff = c(
      sigma(mm_fit),
      as.numeric(logLik(mm_fit)),
      AIC(mm_fit),
      BIC(mm_fit)
    ),
    lme4 = c(
      sigma(lme4_fit),
      as.numeric(stats::logLik(lme4_fit)),
      AIC(lme4_fit),
      BIC(lme4_fit)
    ),
    row.names = NULL,
    check.names = FALSE
  )
}

round_numeric_columns <- function(x, digits = 7) {
  numeric_cols <- vapply(x, is.numeric, logical(1))
  x[numeric_cols] <- lapply(x[numeric_cols], signif, digits = digits)
  x
}

print_case <- function(id, formula, data, reml = TRUE, fixed_basis = NULL) {
  mm_fit <- mixeff::lmm(
    formula, data, REML = reml,
    control = mixeff::mm_control(verbose = -1)
  )
  lme4_fit <- suppressMessages(suppressWarnings(
    lme4::lmer(formula, data = data, REML = reml)
  ))

  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(id, "\n", sep = "")
  cat("Formula: ", deparse1(formula), "\n", sep = "")
  cat("Rows: ", nobs(mm_fit), "; REML: ", reml, "\n", sep = "")

  cat("\nFit statistics\n")
  print(round_numeric_columns(fit_stats(mm_fit, lme4_fit)), row.names = FALSE)

  cat("\nFixed effects\n")
  print(round_numeric_columns(fixed_effects(mm_fit, lme4_fit, fixed_basis)),
        row.names = FALSE)

  cat("\nVarCorr: mixeff\n")
  print(VarCorr(mm_fit))

  invisible(list(mixeff = mm_fit, lme4 = lme4_fit))
}

anchoring <- prepare_anchoring()

print_case(
  "anchoring random intercept",
  everest_feet ~ anchor + (1 | referrer),
  anchoring,
  fixed_basis = "anchor_treatment"
)

print_case(
  "anchoring correlated random slope",
  everest_feet ~ anchor + (1 + anchor | referrer),
  anchoring,
  fixed_basis = "anchor_treatment"
)

print_case(
  "anchoring independent numeric slope",
  everest_feet ~ anchor_contrast + (1 | referrer) +
    (0 + anchor_contrast | referrer),
  anchoring
)

correlated_ml <- print_case(
  "anchoring correlated random slope, ML",
  everest_feet ~ anchor + (1 + anchor | referrer),
  anchoring,
  reml = FALSE,
  fixed_basis = "anchor_treatment"
)

uncorrelated_ml <- print_case(
  "anchoring independent numeric slope, ML",
  everest_feet ~ anchor_contrast + (1 | referrer) +
    (0 + anchor_contrast | referrer),
  anchoring,
  reml = FALSE
)

mm_lrt <- 2 * (
  as.numeric(logLik(correlated_ml$mixeff)) -
    as.numeric(logLik(uncorrelated_ml$mixeff))
)
lme4_lrt <- stats::anova(uncorrelated_ml$lme4, correlated_ml$lme4)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("anchoring likelihood-ratio comparison, ML\n")
print(data.frame(
  engine = c("mixeff", "lme4"),
  chisq = c(mm_lrt, lme4_lrt$Chisq[[2L]]),
  p_value = c(
    stats::pchisq(mm_lrt, df = 1L, lower.tail = FALSE),
    lme4_lrt$`Pr(>Chisq)`[[2L]]
  ),
  check.names = FALSE
), row.names = FALSE)

if (run_slow) {
  speeddate <- prepare_speeddate()
  print_case(
    "speeddate crossed independent random effects",
    other_like ~ other_attr_c * other_intel_c +
      (1 + other_attr_c + other_intel_c || iid) +
      (1 + other_attr_c + other_intel_c || pid),
    speeddate
  )
} else {
  cat(
    "\nSkipping speeddate crossed model. Set MIXEFF_EXAMPLE_RUN_SLOW=true ",
    "to run it.\n",
    sep = ""
  )
}
