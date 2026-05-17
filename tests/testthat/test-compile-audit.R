# Smoke fixture covering Phase 1.A: compile_model() builds an mm_spec and
# audit_design() returns the upstream-rendered audit text for a few worked
# shapes drawn from compiler_contract_v0_prd.md lines 907-922 plus the
# §9.5.7 syntax coverage list. The full snapshot suite (with locked text)
# lands in 1.C alongside explain_model().

mk_sleepstudy_like <- function(n_subjects = 8L, n_per = 4L, seed = 1L) {
  set.seed(seed)
  data.frame(
    Reaction = rnorm(n_subjects * n_per, mean = 250, sd = 30),
    Days     = rep(seq_len(n_per) - 1, n_subjects),
    Subject  = factor(rep(seq_len(n_subjects), each = n_per))
  )
}

mk_crossed <- function(n_subjects = 4L, n_items = 3L, seed = 2L) {
  set.seed(seed)
  grid <- expand.grid(
    Subject = factor(seq_len(n_subjects)),
    Item    = factor(seq_len(n_items))
  )
  grid$y <- rnorm(nrow(grid))
  grid
}

diagnostic_codes <- function(x) {
  vapply(x$diagnostics %||% list(), function(d) d$code, character(1))
}

test_that("compile_model() returns an mm_spec for a sleepstudy-shape design", {
  df <- mk_sleepstudy_like()
  spec <- compile_model(Reaction ~ Days + (1 + Days | Subject), df)

  expect_s3_class(spec, "mm_spec")
  expect_true(is.list(spec$artifact))
  expect_identical(spec$vars, c("Reaction", "Days", "Subject"))
  expect_s3_class(spec$model_frame, "data.frame")
  expect_identical(names(spec$model_frame), spec$vars)

  schema <- spec$artifact$schema
  expect_identical(as.character(schema$schema_name),
                   "mixedmodels.compiled_model_artifact")
  expect_identical(as.character(schema$schema_version), "1")

  # The semantic model recorded one random term grouped by Subject.
  rt <- spec$artifact$semantic_model$random_terms
  expect_length(rt, 1L)
})

test_that("compile_model() handles crossed grouping factors", {
  df <- mk_crossed()
  spec <- compile_model(y ~ 1 + (1 | Subject) + (1 | Item), df)

  expect_s3_class(spec, "mm_spec")
  expect_length(spec$artifact$semantic_model$random_terms, 2L)
})

test_that("audit_design() round-trips upstream text and exposes design_audit", {
  df <- mk_sleepstudy_like()
  spec <- compile_model(Reaction ~ Days + (1 + Days | Subject), df)
  audit <- audit_design(spec)

  expect_s3_class(audit, "mm_audit")
  expect_type(audit$text, "character")
  expect_length(audit$text, 1L)
  expect_match(audit$text, "Audit Summary", fixed = TRUE)
  expect_match(audit$text, "Random Effects",   fixed = TRUE)
  expect_match(audit$text, "Random-Effect Information Budget", fixed = TRUE)

  expect_true(is.list(audit$design_audit))
  expect_true(is.list(audit$report))
  expect_identical(audit$report$schema_name, "mixedmodels.model_audit_report")
  expect_identical(as.character(audit$report$schema_version), "2")
  expect_length(audit$random_term_cards, 1L)
  expect_identical(audit$random_term_cards[[1L]]$schema_name,
                   "mixedmodels.random_term_card")
  expect_identical(audit$random_term_cards[[1L]]$term_id, "r0")
})

test_that("compile/audit surface the upstream pedagogical DiagnosticCode variants", {
  sleep <- mk_sleepstudy_like()
  double_bar <- compile_model(Reaction ~ Days + (1 + Days || Subject), sleep)
  expect_true("covariance_assumption" %in% diagnostic_codes(double_bar$artifact))

  nested <- compile_model(Reaction ~ Days + (1 | Subject/Days), sleep)
  expect_true("syntax_expansion" %in% diagnostic_codes(nested$artifact))

  sparse <- data.frame(
    y       = rnorm(8),
    x       = rep(0:1, 4),
    between = rep(seq_len(4), each = 2),
    g       = factor(rep(seq_len(4), each = 2))
  )

  punt <- compile_model(y ~ x + (1 | g), sparse)
  punt_codes <- diagnostic_codes(audit_design(punt))
  expect_true("support_note" %in% punt_codes)
  expect_true("scope_note" %in% punt_codes)

  refused <- compile_model(y ~ between + (1 + between | g), sparse)
  expect_true("structural_refusal" %in% diagnostic_codes(refused$artifact))
})

# Negative-case tests for scope_note. The upstream generator (audit.rs::
# scope_note_diagnostics) has dedicated Rust tests for the positive trigger
# and split-block suppression; these R tests pin the three remaining filters
# (numeric-only, no-interaction, varies-within-group) so a future Rust change
# cannot regress them silently from the R surface.

test_that("scope_note does not fire for categorical (factor) fixed effects", {
  set.seed(101)
  n_subj <- 8L
  n_per  <- 4L
  # `treat` varies within Subject but is a factor; the upstream generator
  # restricts scope_note candidates to numeric fixed effects.
  df <- data.frame(
    y       = rnorm(n_subj * n_per),
    treat   = factor(rep(c("a", "b"), length.out = n_subj * n_per)),
    Subject = factor(rep(seq_len(n_subj), each = n_per))
  )
  spec  <- compile_model(y ~ treat + (1 | Subject), df)
  codes <- diagnostic_codes(audit_design(spec))
  expect_false("scope_note" %in% codes)
})

test_that("scope_note fires for numeric main effects but not for their interaction", {
  set.seed(102)
  n_subj <- 8L
  n_per  <- 4L
  df <- data.frame(
    y       = rnorm(n_subj * n_per),
    x       = rnorm(n_subj * n_per),
    z       = rnorm(n_subj * n_per),
    Subject = factor(rep(seq_len(n_subj), each = n_per))
  )
  spec  <- compile_model(y ~ x * z + (1 | Subject), df)
  diags <- audit_design(spec)$diagnostics %||% list()
  scope_terms <- vapply(diags, function(d) {
    if (identical(d$code, "scope_note")) d$payload$fixed_effect else NA_character_
  }, character(1))
  scope_terms <- scope_terms[!is.na(scope_terms)]
  expect_true("x" %in% scope_terms)
  expect_true("z" %in% scope_terms)
  expect_false("x:z" %in% scope_terms)
})

test_that("scope_note does not fire for fixed effects constant within group", {
  set.seed(103)
  n_subj <- 8L
  n_per  <- 4L
  # `dose` is a between-subject covariate: one value per Subject, constant
  # across the within-subject rows. The basis.supported check upstream
  # rejects within-group-constant fixed effects as scope_note candidates.
  dose_per_subj <- rnorm(n_subj)
  df <- data.frame(
    y       = rnorm(n_subj * n_per),
    dose    = rep(dose_per_subj, each = n_per),
    Subject = factor(rep(seq_len(n_subj), each = n_per))
  )
  spec  <- compile_model(y ~ dose + (1 | Subject), df)
  codes <- diagnostic_codes(audit_design(spec))
  expect_false("scope_note" %in% codes)
})

test_that("compile_model() refuses unknown formula variables", {
  df <- mk_sleepstudy_like()
  expect_error(
    compile_model(Reaction ~ Missing + (1 | Subject), df),
    class = "mm_data_error",
    regexp = "not found in `data`.*Missing"
  )
})

test_that("compile_model() refuses NA in design variables", {
  df <- mk_sleepstudy_like()
  df$Days[3L] <- NA_real_
  expect_error(
    compile_model(Reaction ~ Days + (1 | Subject), df),
    class = "mm_data_error",
    regexp = "Missing values"
  )
})

test_that("compile_model() refuses one-sided formulas", {
  df <- mk_sleepstudy_like()
  expect_error(
    compile_model(~ Days + (1 | Subject), df),
    class = "mm_formula_error"
  )
})

test_that("compile_model() refuses non-data.frame data", {
  df <- mk_sleepstudy_like()
  expect_error(
    compile_model(Reaction ~ Days + (1 | Subject),
                  list(Reaction = df$Reaction, Days = df$Days,
                       Subject = df$Subject)),
    class = "mm_data_error"
  )
})

test_that("audit_design() refuses inputs that are not mm_spec / mm_fit", {
  expect_error(audit_design(list(foo = 1)), class = "mm_schema_error")
})

test_that("compile_model() preserves input column order through the FFI", {
  # Re-order the data.frame so Subject comes first; the artifact's
  # requested_formula and semantic model should still describe the
  # original formula. The translator preserves data.frame column order.
  df <- mk_sleepstudy_like()
  df <- df[, c("Subject", "Days", "Reaction")]
  spec <- compile_model(Reaction ~ Days + (1 | Subject), df)
  expect_identical(spec$vars, c("Reaction", "Days", "Subject"))
  # `requested_formula` is the upstream canonical Display form (which makes
  # the implicit intercept explicit), not the original input string.
  expect_identical(spec$artifact$requested_formula,
                   mm_parse_formula(Reaction ~ Days + (1 | Subject)))
})
