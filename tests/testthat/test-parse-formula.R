test_that("mm_parse_formula round-trips canonical strings unchanged", {
  expect_identical(
    mm_parse_formula("y ~ 1 + x + (1 | g)"),
    "y ~ 1 + x + (1 | g)"
  )
  expect_identical(
    mm_parse_formula("y ~ 1 + x + (1 | a) + (1 | a:b)"),
    "y ~ 1 + x + (1 | a) + (1 | a:b)"
  )
})

test_that("mm_parse_formula inserts the implicit intercept", {
  # Bare RHS with no explicit `1` should canonicalize with `1 +` prefixed.
  expect_identical(
    mm_parse_formula("y ~ x + (1|g)"),
    "y ~ 1 + x + (1 | g)"
  )
})

test_that("mm_parse_formula expands nested grouping (1|a/b) -> (1|a)+(1|a:b)", {
  expect_identical(
    mm_parse_formula("y ~ x + (1|a/b)"),
    "y ~ 1 + x + (1 | a) + (1 | a:b)"
  )
})

test_that("mm_parse_formula preserves the || zero-correlation marker", {
  expect_identical(
    mm_parse_formula("y ~ x + (1+x||g)"),
    "y ~ 1 + x + (1 + x || g)"
  )
})

test_that("mm_parse_formula accepts an R formula object", {
  expect_identical(
    mm_parse_formula(y ~ x + (1 | g)),
    "y ~ 1 + x + (1 | g)"
  )
})

test_that("mm_parse_formula multi-line formula objects collapse to one line", {
  f <- y ~ x +
    (1 | g) +
    (0 + x | h)
  expect_identical(
    mm_parse_formula(f),
    "y ~ 1 + x + (1 | g) + (0 + x | h)"
  )
})

test_that("mm_parse_formula raises mm_formula_error for invalid input", {
  expect_error(
    mm_parse_formula("y ~ x ((((( z"),
    class = "mm_formula_error"
  )
  expect_error(
    mm_parse_formula("y ~ x @@@@"),
    class = "mm_formula_error"
  )
})

test_that("mm_formula_error inherits from mm_condition for generic catch", {
  # Callers should be able to catch all package conditions via mm_condition
  # without enumerating subclasses.
  caught <- tryCatch(
    mm_parse_formula("y ~ x @@@@"),
    mm_condition = function(cnd) cnd
  )
  expect_s3_class(caught, "mm_condition")
  expect_s3_class(caught, "mm_formula_error")
})

test_that("mm_formula_error carries the offending input on the condition", {
  caught <- tryCatch(
    mm_parse_formula("y ~ x @@@@"),
    mm_formula_error = function(cnd) cnd
  )
  expect_identical(caught$formula, "y ~ x @@@@")
})

test_that("mm_parse_formula validates input shape with a typed condition", {
  expect_error(
    mm_parse_formula(NULL),
    class = "mm_formula_error"
  )
  expect_error(
    mm_parse_formula(NA_character_),
    class = "mm_formula_error"
  )
  expect_error(
    mm_parse_formula(c("y ~ x", "z ~ x")),
    class = "mm_formula_error"
  )
  expect_error(
    mm_parse_formula(""),
    class = "mm_formula_error"
  )
  expect_error(
    mm_parse_formula("   "),
    class = "mm_formula_error"
  )
})
