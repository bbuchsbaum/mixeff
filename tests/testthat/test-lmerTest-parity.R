test_that("scalar fixed-effect contrasts match lmerTest references", {
  mm_skip_if_no_lmerTest()

  cases <- mm_lme4_parity_cases(c(
    "sleepstudy_random_intercept",
    "sleepstudy_random_intercept_slope",
    "dyestuff_random_intercept",
    "penicillin_crossed_intercepts",
    "cake_recipe_temp"
  ))
  methods <- c("satterthwaite", "kenward_roger", "asymptotic")

  for (case in cases) {
    for (method in methods) {
      mm_expect_scalar_lmerTest_parity(case, method)
    }
  }
})

test_that("single-df term rows match lmerTest ANOVA F equivalents", {
  mm_skip_if_no_lmerTest()

  cases <- mm_lme4_parity_cases(c(
    "sleepstudy_random_intercept",
    "sleepstudy_random_intercept_slope"
  ))
  methods <- c("satterthwaite", "kenward_roger")
  sources <- c("test_effect", "anova")

  for (case in cases) {
    for (method in methods) {
      for (source in sources) {
        mm_expect_term_lmerTest_parity(case, "Days", method, source)
      }
    }
  }
})

test_that("multi-df Kenward-Roger term rows match lmerTest ANOVA", {
  mm_skip_if_no_lmerTest()
  mm_skip_if_no_pbkrtest()

  case <- mm_lme4_parity_cases("cake_recipe_temp")[[1L]]
  terms <- c("recipe", "recipe:temp")
  sources <- c("test_effect", "anova")

  for (term in terms) {
    for (source in sources) {
      mm_expect_term_lmerTest_parity(case, term, "kenward_roger", source)
    }
  }
})


test_that("anova(type=...) routes to Rust-owned typed term hypotheses", {
  df <- expand.grid(
    subject = factor(seq_len(18)),
    x = factor(c("x0", "x1")),
    z = factor(c("z0", "z1")),
    KEEP.OUT.ATTRS = FALSE
  )
  df <- df[!(as.integer(df$subject) %% 3L == 0L & df$x == "x1" & df$z == "z1"), ]
  df$subject_u <- seq(-0.7, 0.7, length.out = 18)[as.integer(df$subject)]
  df$y <- 1 +
    0.5 * (df$x == "x1") +
    0.7 * (df$z == "z1") +
    1.1 * (df$x == "x1" & df$z == "z1") +
    df$subject_u +
    seq(-0.2, 0.2, length.out = nrow(df))

  fit <- lmm(y ~ x * z + (1 | subject), df, REML = FALSE,
             control = mm_control(verbose = -1))
  by_type <- list(
    type_i = stats::anova(fit, type = "I", method = "asymptotic"),
    type_ii = stats::anova(fit, type = "II", method = "asymptotic"),
    type_iii = stats::anova(fit, type = "III", method = "asymptotic")
  )

  term_notes <- function(x, term) {
    row <- x$table[x$table$term == term, , drop = FALSE]
    expect_equal(nrow(row), 1L)
    unlist(row$notes[[1L]], use.names = FALSE)
  }

  expect_true(any(grepl("fixed-effect term test type: type_i",
                        term_notes(by_type$type_i, "x"), fixed = TRUE)))
  expect_true(any(grepl("fixed-effect term test type: type_ii",
                        term_notes(by_type$type_ii, "x"), fixed = TRUE)))
  expect_true(any(grepl("fixed-effect term test type: type_iii",
                        term_notes(by_type$type_iii, "x"), fixed = TRUE)))
})
