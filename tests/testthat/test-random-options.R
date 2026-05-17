mk_random_options_design <- function(seed = 4L) {
  set.seed(seed)
  data.frame(
    y = rnorm(36),
    t = rep(seq_len(6), 6),
    s = factor(rep(seq_len(6), each = 6))
  )
}

random_output <- function(x) paste(capture.output(print(x)), collapse = "\n")

test_that("random_options() renders a nearby map with current marker", {
  spec <- compile_model(y ~ t + (1 | s), mk_random_options_design())
  opts <- random_options(spec, group = s)

  expect_s3_class(opts, "mm_random_options")
  expect_s3_class(opts$options, "data.frame")
  expect_identical(
    names(opts$options),
    c("formula", "varying_coefficients", "covariance_family",
      "theta_parameters", "design_status", "plain_meaning", "note", "current")
  )
  expect_false("recommended" %in% tolower(names(opts$options)))
  expect_false(any(tolower(names(opts$options)) %in% c("rank", "score", "preference")))
  expect_identical(opts$options$formula[[1L]], "(1 | s)")
  expect_true(opts$options$current[[1L]])
  expect_true("(0 + t | s)" %in% opts$options$formula)
  expect_true("(1 | s) + (0 + t | s)" %in% opts$options$formula)
  expect_true("(1 + t || s)" %in% opts$options$formula)
  expect_true("(1 + t | s)" %in% opts$options$formula)
  expect_false(any(nzchar(opts$options$note)))
  expect_true(length(opts$constraints$split_uncorrelated) >= 1L)
  expect_true(length(opts$constraints$double_bar_synonym) >= 1L)
  expect_match(opts$options$plain_meaning[opts$options$formula == "(1 + t || s)"],
               "fixes the covariance", fixed = TRUE)

  printed <- random_output(opts)
  expect_match(printed, "Random-effect options for group: s", fixed = TRUE)
  expect_match(printed, "(1 | s) <- this is what you wrote", fixed = TRUE)
  expect_match(printed, "Nearby options:", fixed = TRUE)
})

test_that("random_options() marks split and double-bar spellings by exact current syntax", {
  df <- mk_random_options_design()

  split <- random_options(compile_model(y ~ t + (1 | s) + (0 + t | s), df), "s", "t")
  expect_true(split$options$current[split$options$formula == "(1 | s) + (0 + t | s)"])
  expect_false(split$options$current[split$options$formula == "(1 + t || s)"])

  dbl <- random_options(compile_model(y ~ t + (1 + t || s), df), "s", "t")
  expect_true(dbl$options$current[dbl$options$formula == "(1 + t || s)"])
  expect_false(dbl$options$current[dbl$options$formula == "(1 | s) + (0 + t | s)"])
})

test_that("random_options() stores upstream cards that round-trip as JSON", {
  opts <- random_options(
    compile_model(y ~ t + (1 + t | s), mk_random_options_design()),
    group = "s",
    slope = "t"
  )
  encoded <- jsonlite::toJSON(opts$cards$full, auto_unbox = TRUE)
  decoded <- jsonlite::fromJSON(encoded, simplifyVector = FALSE)

  expect_identical(decoded[[1L]]$schema_name, "mixedmodels.random_term_card")
  expect_identical(decoded[[1L]]$blocks[[1L]]$english,
                   opts$cards$full[[1L]]$blocks[[1L]]$english)
})

test_that("compare_covariance() lays out full diagonal scalar without recommendations", {
  cmp <- compare_covariance(
    compile_model(y ~ t + (1 + t | s), mk_random_options_design())
  )

  expect_s3_class(cmp, "mm_compare_covariance")
  expect_s3_class(cmp$table, "data.frame")
  expect_identical(unique(cmp$table$covariance_family),
                   c("full", "diagonal", "scalar"))
  expect_false("recommended" %in% tolower(names(cmp$table)))
  expect_true(cmp$table$current[cmp$table$covariance_family == "full"])
  expect_equal(cmp$table$theta_parameters[cmp$table$covariance_family == "full"], 3)
  expect_equal(cmp$table$theta_parameters[cmp$table$covariance_family == "diagonal"], 2)
  expect_equal(cmp$table$theta_parameters[cmp$table$covariance_family == "scalar"], 1)
  expect_true("assumes_zero" %in% names(cmp$table))
  expect_true(is.list(cmp$cross_card_constraints))

  printed <- random_output(cmp)
  expect_match(printed, "Covariance comparison:", fixed = TRUE)
  expect_match(printed, "assumes zero:", fixed = TRUE)
})

test_that("random_options() and compare_covariance() keep R9 forbidden advice phrases out", {
  spec <- compile_model(y ~ t + (1 | s), mk_random_options_design())
  output <- tolower(paste(
    random_output(random_options(spec, "s", "t")),
    random_output(compare_covariance(spec)),
    collapse = "\n"
  ))

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

test_that("random_options() refuses unusable inputs clearly", {
  expect_error(random_options(list(), "s"), class = "mm_schema_error")
})
