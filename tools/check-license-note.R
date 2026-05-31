#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
write_mode <- "--write" %in% args

pkg_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
note_path <- file.path(pkg_root, "inst", "LICENSE.note")
vendor_tar <- file.path(pkg_root, "src", "rust", "vendor.tar.xz")
lock_path <- file.path(pkg_root, "src", "rust", "upstream", "mixeff-rs.lock")
upstream_license <- file.path(pkg_root, "src", "rust", "upstream", "mixeff-rs", "LICENSE")

read_lock_value <- function(key) {
  if (!file.exists(lock_path)) return(NA_character_)
  lines <- readLines(lock_path, warn = FALSE)
  pat <- sprintf("^%s\\s*=\\s*\"([^\"]+)\"", key)
  hit <- grep(pat, lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(pat, "\\1", hit[[1L]])
}

field_value <- function(lines, key) {
  in_package <- FALSE
  pat <- sprintf("^%s\\s*=\\s*\"([^\"]+)\"", key)
  for (line in lines) {
    trimmed <- trimws(line)
    if (identical(trimmed, "[package]")) {
      in_package <- TRUE
      next
    }
    if (startsWith(trimmed, "[") && !identical(trimmed, "[package]")) {
      in_package <- FALSE
    }
    if (in_package && grepl(pat, trimmed)) {
      return(sub(pat, "\\1", trimmed))
    }
  }
  NA_character_
}

vendor_crates <- function() {
  if (!file.exists(vendor_tar)) {
    stop("Missing vendor tarball: ", vendor_tar, call. = FALSE)
  }
  entries <- utils::untar(vendor_tar, list = TRUE)
  cargo_entries <- entries[grepl("/Cargo\\.toml$", entries)]
  if (!length(cargo_entries)) {
    stop("vendor.tar.xz contains no vendored Cargo.toml files", call. = FALSE)
  }

  exdir <- tempfile("mixeff-vendor-license-")
  dir.create(exdir)
  on.exit(unlink(exdir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::untar(vendor_tar, files = cargo_entries, exdir = exdir)

  crates <- lapply(cargo_entries, function(entry) {
    lines <- readLines(file.path(exdir, entry), warn = FALSE)
    name <- field_value(lines, "name")
    version <- field_value(lines, "version")
    license <- field_value(lines, "license")
    license_file <- field_value(lines, "license-file")
    crate_dir <- dirname(entry)
    list(
      name = name,
      version = version,
      license = license,
      license_file = license_file,
      dir = crate_dir
    )
  })
  crates <- Filter(function(x) !is.na(x$name) && nzchar(x$name), crates)
  crates <- crates[order(vapply(crates, function(x) x$name, character(1)),
                         vapply(crates, function(x) x$version, character(1)))]
  crates
}

render_note <- function() {
  source <- read_lock_value("source")
  sha <- read_lock_value("resolved_sha")
  crates <- vendor_crates()
  upstream_license_status <- if (file.exists(upstream_license)) {
    "included at src/rust/upstream/mixeff-rs/LICENSE"
  } else {
    "MISSING: src/rust/upstream/mixeff-rs/LICENSE"
  }

  crate_lines <- vapply(crates, function(crate) {
    license <- crate$license
    if (is.na(license) || !nzchar(license)) {
      license <- if (!is.na(crate$license_file) && nzchar(crate$license_file)) {
        sprintf("license-file: %s", crate$license_file)
      } else {
        "license field not declared; see preserved license files"
      }
    }
    sprintf("- %s %s: %s (%s)",
            crate$name, crate$version, license, crate$dir)
  }, character(1))

  c(
    "mixeff license notice",
    "",
    "mixeff is licensed as MIT + file LICENSE.",
    "",
    "Bundled Rust source",
    sprintf("- crate: mixeff-rs"),
    sprintf("- source: %s", source),
    sprintf("- pinned commit: %s", sha),
    sprintf("- upstream license: %s", upstream_license_status),
    "",
    "Vendored Cargo registry crates",
    "",
    "The package vendors Cargo registry dependencies under src/rust/vendor.tar.xz",
    "for offline R CMD INSTALL and CRAN/R-Universe builds. License files from",
    "cargo vendor are preserved inside that archive next to each crate. The",
    "declared crate licenses at generation time are:",
    "",
    crate_lines,
    "",
    "Regeneration",
    "",
    "Run `Rscript tools/vendor-rust.R` after changing the pinned Rust revision,",
    "then run `Rscript tools/check-license-note.R --write` and commit the",
    "updated inst/LICENSE.note."
  )
}

expected <- render_note()

if (write_mode) {
  dir.create(dirname(note_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(expected, note_path)
  cat(sprintf("wrote %s\n", note_path))
  quit(status = 0L)
}

if (!file.exists(note_path)) {
  stop("Missing inst/LICENSE.note; run tools/check-license-note.R --write",
       call. = FALSE)
}
actual <- readLines(note_path, warn = FALSE)
if (!identical(actual, expected)) {
  stop("inst/LICENSE.note is stale; run tools/check-license-note.R --write",
       call. = FALSE)
}
cat("inst/LICENSE.note is current\n")
