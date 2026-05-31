#!/usr/bin/env Rscript

pkg_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
r_bin <- file.path(R.home("bin"), "R")
work <- tempfile("mixeff-offline-install-")
lib <- file.path(work, "lib")
dir.create(work)
dir.create(lib)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

run <- function(args, env = character()) {
  out <- system2(r_bin, args, stdout = TRUE, stderr = TRUE, env = env)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    cat(out, sep = "\n")
    stop(sprintf("command failed: R %s", paste(args, collapse = " ")),
         call. = FALSE)
  }
  out
}

old <- setwd(work)
on.exit(setwd(old), add = TRUE)
run(c("CMD", "build", "--no-build-vignettes", "--no-manual", pkg_root))
tarballs <- list.files(work, pattern = "^mixeff_.*[.]tar[.]gz$", full.names = TRUE)
if (length(tarballs) != 1L) {
  stop("expected one built mixeff tarball, found: ",
       paste(basename(tarballs), collapse = ", "), call. = FALSE)
}

offline_env <- c("CARGO_NET_OFFLINE=true")
run(c("CMD", "INSTALL", "-l", lib, tarballs[[1L]]), env = offline_env)
load_check <- file.path(work, "load-check.R")
writeLines(c(
  sprintf("library(mixeff, lib.loc = %s)", dQuote(lib)),
  "stopifnot(packageVersion('mixeff') >= '0.1.0')"
), load_check)
run(c("--vanilla", "-s", "-f", load_check), env = offline_env)

cat(sprintf("offline source install passed: %s\n", basename(tarballs[[1L]])))
