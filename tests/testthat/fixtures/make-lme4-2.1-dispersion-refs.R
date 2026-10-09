# Regenerate lme4-2.1-dispersion-refs.json: lme4 >= 2.1-0 glmer() references
# for the free-dispersion GLMM cases in ../helper-lme4-dispersion.R, used by
# the tests when the installed lme4 predates 2.1-0 (CI runs lme4 2.0-6).
#
# Run from the package root with lme4 >= 2.1-0 first on the library path:
#   Rscript tests/testthat/fixtures/make-lme4-2.1-dispersion-refs.R
if (utils::packageVersion("lme4") < "2.1-0") {
  stop("lme4 >= 2.1-0 is required to regenerate the dispersion references")
}
test_path <- function(...) file.path("tests", "testthat", ...)
testthat_env <- new.env()
sys.source(test_path("helper-lme4-dispersion.R"), envir = testthat_env)
cases <- names(testthat_env$mm_disp_cases())
refs <- lapply(cases, testthat_env$mm_disp_ref_live)
names(refs) <- cases
out <- list(
  provenance = list(
    generator = "tests/testthat/fixtures/make-lme4-2.1-dispersion-refs.R",
    R = R.version.string,
    lme4 = format(utils::packageVersion("lme4")),
    Matrix = format(utils::packageVersion("Matrix")),
    platform = R.version$platform,
    generated_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ")
  ),
  cases = refs
)
path <- test_path("fixtures", "lme4-2.1-dispersion-refs.json")
jsonlite::write_json(out, path, digits = NA, auto_unbox = TRUE, pretty = TRUE)
cat("wrote", path, "\n")
