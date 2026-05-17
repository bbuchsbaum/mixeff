# tools/vendor-rust.R
#
# Bundle the upstream `mixedmodels` Rust crate into the package source tree
# and vendor every transitive registry dependency, so the resulting tarball
# is self-contained and `R CMD check` (and CRAN / R-Universe) can build the
# package without reaching outside the unpacked tree.
#
# Layout produced (all paths relative to the package root):
#
#   src/rust/upstream/mixedmodels/   <-- bundled snapshot of the upstream crate
#   src/vendor/                      <-- cargo-vendored registry deps
#   src/rust/vendor-config.toml      <-- the [source.crates-io] stanza, copied
#                                        to src/.cargo/config.toml at build
#                                        time by src/Makevars.in
#
# Run before building or `R CMD check`. Re-run whenever the upstream crate or
# its `Cargo.lock` changes.
#
# Environment variables:
#   MIXEFF_UPSTREAM_PATH  override path to the upstream `mixedmodels` repo.
#                         Defaults to /Users/bbuchsbaum/code/rust/mixedmodels.
#
# Exit on first error so partial states are visible.

options(error = function() {
  message("[vendor-rust.R] aborting on error")
  quit(status = 1L, save = "no")
})

# ---- locate upstream ------------------------------------------------------

upstream <- Sys.getenv(
  "MIXEFF_UPSTREAM_PATH",
  unset = "/Users/bbuchsbaum/code/rust/mixedmodels"
)
if (!dir.exists(upstream)) {
  stop(sprintf(
    "Upstream `mixedmodels` not found at: %s\nSet MIXEFF_UPSTREAM_PATH to override.",
    upstream
  ))
}

pkg_root <- normalizePath(".", winslash = "/")
src_rust <- file.path(pkg_root, "src", "rust")
dest_upstream <- file.path(src_rust, "upstream", "mixedmodels")
dest_vendor   <- file.path(pkg_root, "src", "vendor")
dest_cargo    <- file.path(pkg_root, "src", ".cargo")
vendor_config <- file.path(src_rust, "vendor-config.toml")

cat(sprintf("[vendor-rust.R] upstream:    %s\n", upstream))
cat(sprintf("[vendor-rust.R] dest:        %s\n", dest_upstream))
cat(sprintf("[vendor-rust.R] vendor:      %s\n", dest_vendor))

# ---- 1. clean previous outputs --------------------------------------------

for (p in c(dest_upstream, dest_vendor, dest_cargo, vendor_config)) {
  if (file.exists(p) || dir.exists(p)) {
    unlink(p, recursive = TRUE, force = TRUE)
  }
}

# ---- 2. copy upstream snapshot --------------------------------------------

# Whitelist of files/dirs we bundle. Everything not listed (target/, tests/,
# benches/, examples/, comparison/, datasets/, MixedModels.jl/, docs/,
# scripts/, .git/, .mote/, .omc/) is excluded. Compile-time inputs only.
keep <- c("Cargo.toml", "Cargo.lock", "src")

dir.create(dest_upstream, recursive = TRUE, showWarnings = FALSE)
for (entry in keep) {
  src_path <- file.path(upstream, entry)
  if (!file.exists(src_path) && !dir.exists(src_path)) {
    stop(sprintf("[vendor-rust.R] missing required upstream entry: %s", entry))
  }
  ok <- file.copy(
    from      = src_path,
    to        = dest_upstream,
    recursive = TRUE,
    copy.date = TRUE
  )
  if (!isTRUE(ok)) {
    stop(sprintf("[vendor-rust.R] failed to copy %s", entry))
  }
}

# Synthesize a LICENSE file from the Cargo.toml `license` field. Upstream
# does not ship a LICENSE file at the repo root, but its Cargo.toml declares
# `license = "MIT"`. Bundle a minimal MIT notice so transitive license audits
# can find one.
upstream_cargo <- readLines(file.path(dest_upstream, "Cargo.toml"))
license_line <- grep('^\\s*license\\s*=', upstream_cargo, value = TRUE)[1]
license_id <- if (length(license_line) && !is.na(license_line)) {
  sub('^.*"([^"]+)".*$', "\\1", license_line)
} else {
  NA_character_
}
if (identical(license_id, "MIT")) {
  writeLines(c(
    "MIT License",
    "",
    "Copyright (c) the mixedmodels authors",
    "",
    "Permission is hereby granted, free of charge, to any person obtaining a copy",
    "of this software and associated documentation files (the \"Software\"), to deal",
    "in the Software without restriction, including without limitation the rights",
    "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell",
    "copies of the Software, and to permit persons to whom the Software is",
    "furnished to do so, subject to the following conditions:",
    "",
    "The above copyright notice and this permission notice shall be included in",
    "all copies or substantial portions of the Software.",
    "",
    "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR",
    "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,",
    "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE",
    "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER",
    "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,",
    "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN",
    "THE SOFTWARE."
  ), file.path(dest_upstream, "LICENSE-MIT"))
}

cat(sprintf(
  "[vendor-rust.R] bundled upstream snapshot (%s files)\n",
  length(list.files(dest_upstream, recursive = TRUE))
))

# ---- 3. point Cargo.toml at the bundled path ------------------------------

cargo_toml_path <- file.path(src_rust, "Cargo.toml")
cargo_toml <- readLines(cargo_toml_path)
old_path <- "../../../rust/mixedmodels"
new_path <- "upstream/mixedmodels"
hit <- grepl(old_path, cargo_toml, fixed = TRUE)
if (any(hit)) {
  cargo_toml[hit] <- sub(old_path, new_path, cargo_toml[hit], fixed = TRUE)
  writeLines(cargo_toml, cargo_toml_path)
  cat(sprintf("[vendor-rust.R] rewrote path dep to %s\n", new_path))
} else if (any(grepl(new_path, cargo_toml, fixed = TRUE))) {
  cat("[vendor-rust.R] Cargo.toml already references bundled path\n")
} else {
  stop("[vendor-rust.R] could not find mixedmodels path dep in Cargo.toml")
}

# ---- 4. run cargo vendor --------------------------------------------------

# `cargo vendor` writes its config snippet to stdout and the vendored crates
# to the target dir. Run from src/rust so the relative target ../vendor maps
# to src/vendor (which is what src/Makevars.in expects).
cargo <- Sys.which("cargo")
if (!nzchar(cargo)) {
  cargo <- file.path(Sys.getenv("HOME"), ".cargo", "bin", "cargo")
}
if (!file.exists(cargo)) stop("[vendor-rust.R] cargo not found on PATH")

cat("[vendor-rust.R] running cargo vendor (this can take a minute)...\n")
res <- system2(
  command = cargo,
  args    = c("vendor", "--manifest-path", file.path(src_rust, "Cargo.toml"),
              file.path(pkg_root, "src", "vendor")),
  stdout  = TRUE,
  stderr  = TRUE
)
if (!is.null(attr(res, "status")) && attr(res, "status") != 0L) {
  cat(res, sep = "\n")
  stop("[vendor-rust.R] cargo vendor failed")
}

# `cargo vendor` writes a config snippet to stdout, but it embeds whatever
# vendor path we passed (absolute, in our case). The Makevars copies our
# vendor-config.toml to src/.cargo/config.toml at build time, and cargo
# then resolves the directory relative to the config file's location. From
# src/.cargo/, the vendor dir at src/vendor/ is `../vendor`. Always write
# the canonical relative form, regardless of what cargo printed — the
# captured stdout is consulted only as a sanity check.
canonical_config <- c(
  "[source.crates-io]",
  'replace-with = "vendored-sources"',
  "",
  "[source.vendored-sources]",
  "# Resolved relative to CARGO_HOME's parent (= the package src/ dir at",
  "# build time per src/Makevars.in: CARGO_HOME=src/.cargo, vendored sources",
  "# at src/vendor). Survives tarball relocation because both are inside",
  "# src/.",
  'directory = "vendor"'
)
writeLines(canonical_config, vendor_config)
cat(sprintf("[vendor-rust.R] wrote %s\n", vendor_config))

# ---- 5. write src/rust/vendor.tar.xz for tarball distribution -------------

# `R CMD build` runs the package's Makevars `clean` target, which (per the
# rextendr template) removes src/vendor/ as a build artifact. The directory
# form is correct for dev workflows (devtools::test), but the *tarball*
# distributed via R CMD INSTALL needs an immutable form that survives the
# build cleanup. The Makevars `if [ -d ./vendor ] ... elif [ -f
# ./rust/vendor.tar.xz ]` block honors either; we ship both:
#   - src/vendor/                  -> in-place dev (faster: no untar step)
#   - src/rust/vendor.tar.xz       -> tarball distribution (survives cleanup)
#
# The tarball must contain entries with a leading `vendor/` segment because
# `tar xf` is invoked from src/ at install time.

vendor_tarball <- file.path(src_rust, "vendor.tar.xz")
if (file.exists(vendor_tarball)) unlink(vendor_tarball, force = TRUE)
src_dir <- file.path(pkg_root, "src")
old_wd <- setwd(src_dir)
on.exit(setwd(old_wd), add = TRUE)
tar_status <- system2(
  command = "tar",
  args    = c("-cJf", vendor_tarball, "vendor"),
  stdout  = TRUE,
  stderr  = TRUE
)
if (!is.null(attr(tar_status, "status")) && attr(tar_status, "status") != 0L) {
  cat(tar_status, sep = "\n")
  stop("[vendor-rust.R] failed to create vendor.tar.xz")
}
setwd(old_wd)
on.exit()
cat(sprintf(
  "[vendor-rust.R] wrote %s (%.1f MB)\n",
  vendor_tarball,
  file.info(vendor_tarball)$size / 1024 / 1024
))

# ---- 6. summary -----------------------------------------------------------

vendored_crates <- list.dirs(dest_vendor, recursive = FALSE, full.names = FALSE)
cat(sprintf(
  "[vendor-rust.R] done. %d vendored crates in src/vendor/\n",
  length(vendored_crates)
))
cat("[vendor-rust.R] next: devtools::check(document = FALSE)\n")
