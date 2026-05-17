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
