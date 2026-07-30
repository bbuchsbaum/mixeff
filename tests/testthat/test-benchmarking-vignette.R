test_that("benchmarking vignette stays light and points at benchmark scripts", {
  skip_on_cran() # duplicates R CMD check's own vignette build
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available(), "pandoc is unavailable")

  path <- test_path("../../vignettes/benchmarking.Rmd")
  # Source-tree-only test: the vignette source is not on the installed
  # tests path during R CMD check, only during devtools::test().
  skip_if_not(file.exists(path), "vignette source path is not available")
  text <- readLines(path, warn = FALSE)
  expect_true(any(grepl("inst/benchmarks/lme4-scaling.R", text, fixed = TRUE)))
  expect_true(any(grepl("inst/benchmarks/bootstrap-inference.R", text, fixed = TRUE)))
  expect_true(any(grepl("eval = FALSE", text, fixed = TRUE)))

  out <- rmarkdown::render(path, output_file = tempfile(fileext = ".html"),
                           quiet = TRUE)
  expect_true(file.exists(out))
})
