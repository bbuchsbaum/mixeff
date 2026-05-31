#!/usr/bin/env Rscript

pkg_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
vendor_script <- file.path(pkg_root, "tools", "vendor-rust.R")
lock_path <- file.path(pkg_root, "src", "rust", "upstream", "mixeff-rs.lock")
cargo_toml <- file.path(pkg_root, "src", "rust", "Cargo.toml")
vendor_config <- file.path(pkg_root, "src", "rust", "vendor-config.toml")
vendor_tar <- file.path(pkg_root, "src", "rust", "vendor.tar.xz")
upstream <- file.path(pkg_root, "src", "rust", "upstream", "mixeff-rs")

fail <- function(...) stop(sprintf(...), call. = FALSE)

read_value <- function(path, key) {
  if (!file.exists(path)) fail("missing file: %s", path)
  lines <- readLines(path, warn = FALSE)
  pat <- sprintf("^%s\\s*=\\s*\"?([^\"#]+)\"?", key)
  hit <- grep(pat, lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  trimws(sub(pat, "\\1", hit[[1L]]))
}

vendor_lines <- readLines(vendor_script, warn = FALSE)
pin_line <- grep("^PINNED_REV\\s*<-", vendor_lines, value = TRUE)
if (!length(pin_line)) fail("tools/vendor-rust.R does not declare PINNED_REV")
pinned <- sub('^PINNED_REV\\s*<-\\s*"([^"]+)".*$', "\\1", pin_line[[1L]])
resolved <- read_value(lock_path, "resolved_sha")
if (!identical(pinned, resolved)) {
  fail("PINNED_REV (%s) does not match mixeff-rs.lock resolved_sha (%s)",
       pinned, resolved)
}

for (path in c(cargo_toml, vendor_config, vendor_tar,
               file.path(upstream, "Cargo.toml"),
               file.path(upstream, "Cargo.lock"),
               file.path(upstream, "LICENSE"),
               file.path(upstream, "src"))) {
  if (!file.exists(path) && !dir.exists(path)) fail("missing vendored input: %s", path)
}

cargo_lines <- readLines(cargo_toml, warn = FALSE)
dep_line <- grep("^\\s*mixeff-rs\\s*=", cargo_lines, value = TRUE)
if (!length(dep_line) || !grepl("upstream/mixeff-rs", dep_line[[1L]], fixed = TRUE)) {
  fail("src/rust/Cargo.toml must depend on upstream/mixeff-rs")
}
if (!grepl("unstable-internals", paste(dep_line, collapse = " "), fixed = TRUE)) {
  fail("src/rust/Cargo.toml mixeff-rs dependency must enable unstable-internals")
}

config_lines <- readLines(vendor_config, warn = FALSE)
if (!any(grepl('replace-with = "vendored-sources"', config_lines, fixed = TRUE)) ||
    !any(grepl('directory = "vendor"', config_lines, fixed = TRUE))) {
  fail("src/rust/vendor-config.toml does not point Cargo at the vendored sources")
}

entries <- utils::untar(vendor_tar, list = TRUE)
if (!any(entries %in% c("vendor", "vendor/"))) fail("vendor.tar.xz lacks vendor/ root")
if (!any(grepl("^vendor/.+/Cargo\\.toml$", entries))) {
  fail("vendor.tar.xz contains no vendored Cargo.toml files")
}
if (!any(grepl("^vendor/.+/(LICENSE|COPYING|UNLICENSE)", entries))) {
  fail("vendor.tar.xz contains no recognizable license files")
}
if (dir.exists(file.path(upstream, ".git"))) {
  fail("vendored upstream snapshot must not contain a .git directory")
}

cat(sprintf("vendor snapshot current at %s with %d archive entries\n",
            substr(resolved, 1, 12), length(entries)))
