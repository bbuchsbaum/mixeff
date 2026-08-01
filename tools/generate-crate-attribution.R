#!/usr/bin/env Rscript

# Generate inst/AUTHORS and inst/COPYRIGHTS from the vendored Rust sources
# (Audit1 WI-6.1; CRAN policy on bundled Rust code requires authorship and
# copyright information for the bundled crates, including dependencies).
#
#   inst/AUTHORS     one entry per bundled crate: name, version, declared
#                    authors (from Cargo.toml [package].authors; crates that
#                    omit the field are attributed to "the <name> crate
#                    developers" with the repository URL), license.
#   inst/COPYRIGHTS  the distinct "Copyright ..." notice lines preserved in
#                    each crate's shipped LICENSE/COPYING/NOTICE files. The
#                    complete license texts ship inside src/rust/vendor.tar.xz
#                    next to each crate (cargo vendor preserves them);
#                    inst/LICENSE.note records the per-crate license set.
#
# Cargo author fields are not guaranteed complete (CRAN's own caveat), so
# both outputs point back at the preserved in-archive texts as the
# authoritative record.
#
# Usage:
#   Rscript tools/generate-crate-attribution.R           # (re)write both files
#   Rscript tools/generate-crate-attribution.R --check   # fail if stale
#
# Re-run (via tools/vendor-rust.R, which calls this) on every PINNED_REV
# bump and commit the outputs.

`%||%` <- function(x, y) if (is.null(x)) y else x

args <- commandArgs(TRUE)
check_only <- "--check" %in% args

root <- normalizePath(".")
archive <- file.path(root, "src", "rust", "vendor.tar.xz")
upstream_dir <- file.path(root, "src", "rust", "upstream", "mixeff-rs")
if (!file.exists(archive)) {
  stop("src/rust/vendor.tar.xz not found; run tools/vendor-rust.R first.")
}

work <- file.path(tempdir(), "mixeff-attribution-vendor")
unlink(work, recursive = TRUE, force = TRUE)
dir.create(work, recursive = TRUE)
untar(archive, exdir = work)
vendor_root <- file.path(work, "vendor")
if (!dir.exists(vendor_root)) {
  # Archive may unpack the crate dirs at top level.
  vendor_root <- work
}
crate_dirs <- list.dirs(vendor_root, recursive = FALSE)
crate_dirs <- crate_dirs[file.exists(file.path(crate_dirs, "Cargo.toml"))]
if (!length(crate_dirs)) {
  stop("No vendored crates found inside vendor.tar.xz.")
}

toml_field <- function(lines, field) {
  # Minimal TOML scraping for [package] scalar/array fields; adequate for
  # cargo-generated manifests. Arrays: authors = ["a", "b"].
  pkg_start <- grep("^\\[package\\]", lines)
  if (!length(pkg_start)) return(character())
  section_ends <- grep("^\\[", lines)
  pkg_end <- min(c(section_ends[section_ends > pkg_start[1L]] - 1L,
                   length(lines)))
  block <- lines[pkg_start[1L]:pkg_end]
  hit <- grep(sprintf("^%s\\s*=", field), block, value = TRUE)
  if (!length(hit)) return(character())
  value <- sub(sprintf("^%s\\s*=\\s*", field), "", hit[1L])
  if (startsWith(value, "[")) {
    # Single-line array (cargo normalizes vendored manifests to this form).
    inner <- sub("^\\[", "", sub("\\]\\s*$", "", value))
    items <- regmatches(inner, gregexpr('"[^"]*"', inner))[[1L]]
    gsub('"', "", items, fixed = TRUE)
  } else {
    gsub('"', "", value, fixed = TRUE)
  }
}

first_or <- function(x, default) {
  if (length(x) && !is.na(x[1L]) && nzchar(x[1L])) x[1L] else default
}

crate_info <- function(dir) {
  lines <- readLines(file.path(dir, "Cargo.toml"), warn = FALSE)
  license <- first_or(toml_field(lines, "license"), "")
  if (!nzchar(license)) {
    license_file <- first_or(toml_field(lines, "license-file"), "")
    license <- if (nzchar(license_file)) {
      paste("license-file:", license_file)
    } else {
      "no declared license field; see crate sources"
    }
  }
  list(
    name = first_or(toml_field(lines, "name"), basename(dir)),
    version = first_or(toml_field(lines, "version"), ""),
    authors = toml_field(lines, "authors"),
    license = license,
    repository = first_or(toml_field(lines, "repository"), ""),
    dir = dir
  )
}

copyright_lines <- function(dir) {
  files <- list.files(dir, pattern = "^(LICEN[SC]E|COPYING|NOTICE)",
                      ignore.case = TRUE, full.names = TRUE, recursive = FALSE)
  files <- files[!dir.exists(files)]
  hits <- character()
  for (f in files) {
    text <- tryCatch(readLines(f, warn = FALSE, encoding = "UTF-8"),
                     error = function(e) character())
    text <- iconv(text, to = "UTF-8", sub = "")
    m <- trimws(grep("copyright", text, ignore.case = TRUE, value = TRUE))
    # Keep actual notices ("Copyright (c) 2019 Jane Doe", "Copyright 2014
    # The Rust Project Developers"), not license-body prose or the Apache
    # template placeholder line.
    m <- m[grepl("^Copyright\\s+(\\(c\\)|©|[0-9]{4})", m,
                 ignore.case = TRUE)]
    m <- m[!grepl("\\[yyyy\\]", m, fixed = TRUE)]
    hits <- c(hits, m)
  }
  unique(hits)
}

infos <- lapply(crate_dirs, crate_info)
# radix = C-locale collation, so the committed output is identical across
# machines regardless of LC_COLLATE (en_* vs C sort "-" and "_" differently,
# which would make the vendor-drift comparison machine-dependent).
infos <- infos[order(vapply(infos, function(i) i$name, character(1)),
                     method = "radix")]

upstream_lines <- readLines(file.path(upstream_dir, "Cargo.toml"), warn = FALSE)
upstream_authors <- toml_field(upstream_lines, "authors")
upstream_version <- toml_field(upstream_lines, "version")[1L] %||% ""
upstream_license <- toml_field(upstream_lines, "license")[1L] %||% "MIT"

fmt_authors <- function(info) {
  if (length(info$authors)) {
    # Strip email addresses (CRAN prefers names; the manifests keep them).
    paste(sub("\\s*<[^>]*>", "", info$authors), collapse = ", ")
  } else if (nzchar(info$repository)) {
    sprintf("The %s crate developers (%s)", info$name, info$repository)
  } else {
    sprintf("The %s crate developers", info$name)
  }
}

authors_out <- c(
  "Authors of mixeff and the bundled Rust sources",
  "==============================================",
  "",
  "mixeff (R package): Bradley Buchsbaum",
  sprintf("mixeff-rs %s (bundled engine, src/rust/upstream/): Bradley Buchsbaum (%s)",
          upstream_version, upstream_license),
  "",
  sprintf("Vendored Cargo registry dependencies (src/rust/vendor.tar.xz, %d crates).",
          length(infos)),
  "Declared authors come from each crate's Cargo.toml; the manifests are",
  "not guaranteed to list every contributor, so the repository URLs and the",
  "license texts preserved inside the archive are the authoritative record.",
  "",
  vapply(infos, function(i) {
    sprintf("- %s %s: %s [%s]", i$name, i$version, fmt_authors(i), i$license)
  }, character(1)),
  ""
)

copyrights_out <- c(
  "Copyright notices for mixeff and the bundled Rust sources",
  "==========================================================",
  "",
  "mixeff (R package) and the bundled mixeff-rs engine:",
  sprintf("  Copyright (c) Bradley Buchsbaum (%s + file LICENSE)", upstream_license),
  "",
  "Vendored Cargo registry dependencies: the notices below are the distinct",
  "'Copyright ...' lines preserved in each crate's shipped LICENSE/COPYING/",
  "NOTICE files (cargo vendor keeps the full texts next to each crate inside",
  "src/rust/vendor.tar.xz; per-crate license identifiers are recorded in",
  "inst/LICENSE.note). Crates listing no explicit notice line publish under",
  "their declared license with copyright retained by their authors (see",
  "inst/AUTHORS).",
  "",
  unlist(lapply(infos, function(i) {
    notices <- copyright_lines(i$dir)
    if (!length(notices)) return(character())
    c(sprintf("%s %s:", i$name, i$version),
      paste0("  ", notices),
      "")
  })),
  ""
)

write_or_check <- function(lines, path) {
  if (check_only) {
    current <- if (file.exists(path)) readLines(path, warn = FALSE) else NULL
    if (!identical(current, as.character(lines))) {
      stop(sprintf("%s is stale; re-run tools/generate-crate-attribution.R", path))
    }
    message(sprintf("%s is current", path))
  } else {
    writeLines(lines, path)
    message(sprintf("wrote %s (%d lines)", path, length(lines)))
  }
}

write_or_check(authors_out, file.path(root, "inst", "AUTHORS"))
write_or_check(copyrights_out, file.path(root, "inst", "COPYRIGHTS"))
unlink(work, recursive = TRUE, force = TRUE)
