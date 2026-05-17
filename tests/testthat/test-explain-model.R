mk_explain_design <- function(seed = 3L) {
  set.seed(seed)
  df <- expand.grid(
    s = factor(seq_len(6)),
    i = factor(seq_len(4)),
    b = factor(seq_len(2))
  )
  df$a <- factor(rep(seq_len(3), length.out = nrow(df)))
  df$t <- rep(seq_len(4), length.out = nrow(df))
  df$y <- rnorm(nrow(df))
  df
}

explain_text <- function(formula, data = mk_explain_design()) {
  explain_model(compile_model(formula, data))$text
}

test_that("explain_model() returns a printable mm_explanation", {
  x <- explain_model(compile_model(y ~ t + (1 | s), mk_explain_design()))

  expect_s3_class(x, "mm_explanation")
  expect_type(x$text, "character")
  expect_length(x$text, 1L)
  expect_true(is.list(x$cards))
  expect_true(is.list(x$report))

  printed <- paste(capture.output(print(x)), collapse = "\n")
  expect_match(printed, "Random effects explanation", fixed = TRUE)
})

test_that("explain_model() covers the eight Phase 1.C syntax patterns", {
  texts <- list(
    punt        = explain_text(y ~ t + (1 | s)),
    slope       = explain_text(y ~ t + (0 + t | s)),
    full        = explain_text(y ~ t + (1 + t | s)),
    double_bar  = explain_text(y ~ t + (1 + t || s)),
    split       = explain_text(y ~ t + (1 | s) + (0 + t | s)),
    nested      = explain_text(y ~ t + (1 | a / b)),
    interaction = explain_text(y ~ t + (1 | s:i)),
    crossed     = explain_text(y ~ t + (1 | s) + (1 | i))
  )

  expect_match(texts$punt,
               're(group = s, intercept = TRUE, slopes = NULL, cov = "scalar")',
               fixed = TRUE)
  expect_match(texts$punt, "scope_note:", fixed = TRUE)

  expect_match(texts$slope,
               're(group = s, intercept = FALSE, slopes = t, cov = "scalar")',
               fixed = TRUE)
  expect_match(texts$slope, "`s` units may differ in their `t` slope.",
               fixed = TRUE)

  expect_match(texts$full,
               're(group = s, intercept = TRUE, slopes = t, cov = "full")',
               fixed = TRUE)
  expect_match(texts$full, "theta parameters: 3", fixed = TRUE)

  expect_match(texts$double_bar, "s has 2 separate random-effect blocks.",
               fixed = TRUE)
  expect_match(texts$double_bar,
               "double-bar syntax fixes the covariance between `Intercept` and `t` to zero.",
               fixed = TRUE)
  expect_match(texts$double_bar, "covariance_assumption:", fixed = TRUE)

  expect_match(texts$split, "separate random-effect blocks", fixed = TRUE)
  expect_match(texts$split,
               "r0 <-> r1 (Intercept <-> t): separate random-effect blocks fix the covariance",
               fixed = TRUE)
  expect_match(texts$split, "separate random-effect blocks fix the covariance",
               fixed = TRUE)

  expect_match(texts$nested, "canonical:  (1 | a)", fixed = TRUE)
  expect_match(texts$nested, "canonical:  (1 | a:b)", fixed = TRUE)
  expect_match(texts$nested, "syntax_expansion:", fixed = TRUE)

  expect_match(texts$interaction, "re(group = s:i", fixed = TRUE)

  expect_match(texts$crossed, "re(group = s", fixed = TRUE)
  expect_match(texts$crossed, "re(group = i", fixed = TRUE)
  expect_false(grepl("No random slopes were added", texts$crossed, fixed = TRUE))
})

test_that("explain_model() renders structural refusals as possible repairs", {
  df <- data.frame(
    y       = rnorm(8),
    x       = rep(0:1, 4),
    between = rep(seq_len(4), each = 2),
    g       = factor(rep(seq_len(4), each = 2))
  )

  text <- explain_text(y ~ between + (1 + between | g), df)

  expect_match(text, "structural_refusal", fixed = TRUE)
  expect_match(text, "Possible repairs, not applied automatically:", fixed = TRUE)
  expect_match(text, "`between` does not vary within `g`", fixed = TRUE)
})

test_that("explain_model() keeps the R9 forbidden advice phrases out of output", {
  texts <- c(
    explain_text(y ~ t + (1 | s)),
    explain_text(y ~ t + (1 + t || s)),
    explain_text(y ~ t + (1 | s) + (0 + t | s)),
    explain_text(y ~ t + (1 | a / b))
  )
  output <- tolower(paste(texts, collapse = "\n"))

  forbidden <- c(
    "suggested starting model",
    "we recommend",
    "you should",
    "try .* instead",
    "drop the random slope",
    "same model, different font"
  )
  for (pattern in forbidden) {
    expect_false(grepl(pattern, output, perl = TRUE),
                 info = sprintf("forbidden phrase matched: %s", pattern))
  }
})

test_that("explain_model() refuses non-compiled inputs", {
  expect_error(explain_model(list()), class = "mm_schema_error")
})
