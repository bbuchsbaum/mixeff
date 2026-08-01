# Model-support registry consistency (Audit1 WI-4.1/4.2).
#
# inst/support/model-support.json is the single source of truth for the
# supported envelope. These tests hold three surfaces to it: the fitting
# code's own family/link table, the refusal behavior of unsupported
# requests, and the roxygen prose that enumerates the surface. Drift in any
# of them fails here rather than shipping as contradictory documentation.

test_that("supported_models() reflects the shipped registry", {
  sm <- supported_models()
  expect_s3_class(sm, "mm_model_support")
  expect_true(all(c("glmm_families", "features", "non_goals", "registry")
                  %in% names(sm)))
  expect_identical(sm$registry$schema_name, "mixeff.model_support")
  expect_true(nrow(sm$glmm_families) >= 4L)
  expect_true(nrow(sm$features) >= 5L)
  expect_output(print(sm), "supported model envelope")
})

test_that("the registry and the fitting code agree on family/link cells", {
  sm <- supported_models()
  code_table <- mixeff:::mm_glmm_supported_family_links()
  reg <- setNames(
    lapply(sm$registry$glmm_families, function(f) sort(unlist(f$links))),
    vapply(sm$registry$glmm_families, function(f) f$family, character(1))
  )
  expect_setequal(names(reg), names(code_table))
  for (fam in names(reg)) {
    expect_identical(reg[[fam]], sort(unname(unlist(code_table[[fam]]))),
                     info = sprintf("family `%s` links diverge", fam))
  }
})

test_that("unsupported family/link requests refuse and name the supported set", {
  set.seed(2)
  d <- data.frame(y = rbinom(60, 1, 0.5), x = rnorm(60),
                  g = factor(rep(1:6, each = 10)))
  err <- tryCatch(
    glmm(y ~ x + (1 | g), d, family = inverse.gaussian(),
         control = mm_control(verbose = -1)),
    error = function(e) e
  )
  expect_s3_class(err, "mm_inference_unavailable")
  expect_match(conditionMessage(err), "Supported families", fixed = TRUE)

  err2 <- tryCatch(
    glmm(y ~ x + (1 | g), d, family = binomial("cauchit"),
         control = mm_control(verbose = -1)),
    error = function(e) e
  )
  expect_s3_class(err2, "mm_inference_unavailable")
})

test_that("registry-declared refusals carry their declared classes and codes", {
  set.seed(3)
  d <- data.frame(y = rbinom(60, 1, 0.5), x = rnorm(60),
                  g = factor(rep(1:6, each = 10)))
  # NB + joint_laplace: declared reason code.
  err <- tryCatch(
    glmm(y ~ x + (1 | g), d, family = mm_negative_binomial(theta = 1),
         method = "joint_laplace", control = mm_control(verbose = -1)),
    error = function(e) e
  )
  expect_s3_class(err, "mm_inference_unavailable")

  fit <- glmm(y ~ x + (1 | g), d, family = binomial(),
              control = mm_control(verbose = -1))
  # Feature refusals: classes declared in the registry.
  expect_error(simulate(fit), class = "mm_inference_unavailable")
  expect_error(refit(fit, fitted(fit)), class = "mm_inference_unavailable")
  expect_error(
    glmm(y ~ x + (1 | g), d, family = binomial(), subset = 1:50,
         control = mm_control(verbose = -1)),
    class = "mm_fit_error"
  )
  expect_error(mm_means(fit, ~x), class = "mm_inference_unavailable")
})

test_that("the glmm() roxygen family prose matches the registry", {
  # Source-tree-only guard against prose drift: the enumerated family/link
  # text in R/glmm.R must mention every registry family and link.
  src <- test_path("../../R/glmm.R")
  skip_if_not(file.exists(src), "package source is not available")
  text <- paste(readLines(src, warn = FALSE)[1:60], collapse = " ")
  sm <- supported_models()
  for (f in sm$registry$glmm_families) {
    for (link in unlist(f$links)) {
      expect_match(text, link, fixed = TRUE,
                   info = sprintf("link `%s` missing from glmm() roxygen", link))
    }
  }
  expect_match(text, "negative binomial", ignore.case = TRUE)
})

test_that("the migration vignette renders its support table from the registry", {
  src <- test_path("../../vignettes/lme4-migration.Rmd")
  skip_if_not(file.exists(src), "vignette source is not available")
  text <- readLines(src, warn = FALSE)
  expect_true(any(grepl("supported_models()", text, fixed = TRUE)))
})
