# tools/vendor-rust.R
#
# Bundle a *pinned* snapshot of the upstream `mixeff-rs` Rust crate into the
# package source tree and vendor every transitive registry dependency, so the
# resulting tarball is self-contained and `R CMD check` (and CRAN /
# R-Universe) can build the package without reaching outside the unpacked
# tree.
#
# Provenance model (see CLAUDE.md: "JSON artifacts are source of truth,
# external pointers are caches"):
#
#   * The GitHub repo `bbuchsbaum/mixeff-rs` is the canonical provenance.
#   * The vendored snapshot under src/rust/upstream/mixeff-rs/ is a cache.
#   * That cache is always extracted from a single pinned commit SHA
#     (PINNED_REV below), never from a working-tree copy, so the bundled
#     crate is reproducible and auditable.
#
# The local peer checkout (MIXEFF_RS_PATH) is used only as a git object
# cache to avoid the network when it already contains PINNED_REV; the tree
# is still materialized via `git archive PINNED_REV`, not by copying the
# working directory.
#
# Layout produced (all paths relative to the package root):
#
#   src/rust/upstream/mixeff-rs/        <-- snapshot of the crate @ PINNED_REV
#   src/rust/upstream/mixeff-rs.lock    <-- provenance manifest (url + sha)
#   src/vendor/                         <-- cargo-vendored registry deps
#   src/rust/vendor-config.toml         <-- [source.crates-io] stanza, copied
#                                           to src/.cargo/config.toml at build
#                                           time by src/Makevars.in
#   src/rust/vendor.tar.xz              <-- immutable vendor form for tarballs
#
# Run before building or `R CMD check`. Re-run whenever PINNED_REV changes.
#
# To bump the pin: update PINNED_REV to a commit SHA (or tag) that exists on
# `origin/main` of bbuchsbaum/mixeff-rs, then re-run this script and commit
# the resulting Cargo.lock / parity changes.
#
# Environment variables (all optional; defaults are the committed pin):
#   MIXEFF_RS_REV   commit SHA or tag to vendor. Default: PINNED_REV.
#   MIXEFF_RS_URL   git URL to clone if the rev is not already local.
#   MIXEFF_RS_PATH  local peer checkout used as a git object cache.
#
# Exit on first error so partial states are visible.

options(error = function() {
  message("[vendor-rust.R] aborting on error")
  quit(status = 1L, save = "no")
})

# ---- the pin --------------------------------------------------------------
#
# This constant is the committed source of truth for which `mixeff-rs` ships.
# It must be a full 40-char commit SHA reachable from origin/main of
# bbuchsbaum/mixeff-rs (or a tag, once the crate starts tagging releases).
PINNED_REV <- "4a2abb39bb21901541285ae90d9e4159a02930c1"

rev <- Sys.getenv("MIXEFF_RS_REV", unset = PINNED_REV)
url <- Sys.getenv(
  "MIXEFF_RS_URL",
  unset = "https://github.com/bbuchsbaum/mixeff-rs.git"
)
peer <- Sys.getenv(
  "MIXEFF_RS_PATH",
  unset = "/Users/bbuchsbaum/code/rust/mixeff-rs"
)

pkg_root      <- normalizePath(".", winslash = "/")
src_rust      <- file.path(pkg_root, "src", "rust")
upstream_root <- file.path(src_rust, "upstream")
dest_upstream <- file.path(upstream_root, "mixeff-rs")
manifest_path <- file.path(upstream_root, "mixeff-rs.lock")
dest_vendor   <- file.path(pkg_root, "src", "vendor")
dest_cargo    <- file.path(pkg_root, "src", ".cargo")
vendor_config <- file.path(src_rust, "vendor-config.toml")

git <- Sys.which("git")
if (!nzchar(git)) stop("[vendor-rust.R] git not found on PATH")

run_git <- function(args, dir = NULL, ...) {
  full <- if (is.null(dir)) args else c("-C", dir, args)
  system2(git, full, stdout = TRUE, stderr = TRUE, ...)
}
git_ok <- function(res) {
  st <- attr(res, "status")
  is.null(st) || identical(st, 0L)
}

cat(sprintf("[vendor-rust.R] requested rev: %s\n", rev))
cat(sprintf("[vendor-rust.R] provenance:    %s\n", url))

# ---- 1. resolve a git repo that contains `rev` ----------------------------
#
# Preference order, all of which yield the *same* tree because we extract by
# SHA, not by working copy:
#   (a) the local peer checkout, if it already has the rev (no network);
#   (b) the local peer checkout after `git fetch` (uses its remote);
#   (c) a fresh clone of MIXEFF_RS_URL into a tempdir.

has_rev <- function(dir) {
  git_ok(run_git(c("cat-file", "-e", paste0(rev, "^{commit}")), dir = dir))
}

src_repo <- NULL
clone_dir <- NULL

if (dir.exists(file.path(peer, ".git"))) {
  if (has_rev(peer)) {
    cat(sprintf("[vendor-rust.R] using local peer (has rev): %s\n", peer))
    src_repo <- peer
  } else {
    cat("[vendor-rust.R] peer missing rev; git fetch...\n")
    run_git(c("fetch", "--quiet", "--all", "--tags"), dir = peer)
    if (has_rev(peer)) {
      cat(sprintf("[vendor-rust.R] using local peer (post-fetch): %s\n", peer))
      src_repo <- peer
    }
  }
}

if (is.null(src_repo)) {
  clone_dir <- file.path(tempdir(), "mixeff-rs-pin")
  unlink(clone_dir, recursive = TRUE, force = TRUE)
  cat(sprintf("[vendor-rust.R] cloning %s ...\n", url))
  res <- run_git(c("clone", "--quiet", url, clone_dir))
  if (!git_ok(res)) {
    cat(res, sep = "\n")
    stop("[vendor-rust.R] git clone failed")
  }
  if (!has_rev(clone_dir)) {
    # Rev may be a non-default-branch commit; fetch it explicitly.
    run_git(c("fetch", "--quiet", "origin", rev), dir = clone_dir)
  }
  if (!has_rev(clone_dir)) {
    stop(sprintf(
      "[vendor-rust.R] rev %s not found in %s after clone+fetch", rev, url
    ))
  }
  src_repo <- clone_dir
}

resolved_sha <- run_git(
  c("rev-parse", paste0(rev, "^{commit}")), dir = src_repo
)[1]
if (!grepl("^[0-9a-f]{40}$", resolved_sha)) {
  stop(sprintf("[vendor-rust.R] could not resolve rev to a SHA: %s", rev))
}
cat(sprintf("[vendor-rust.R] resolved SHA:  %s\n", resolved_sha))

# ---- 2. clean previous outputs --------------------------------------------
#
# `upstream_root` is wiped wholesale (not just `dest_upstream`/`manifest_path`)
# so any stray hand-placed or pre-rename directory next to the snapshot —
# e.g. an old `upstream/mixedmodels/` — cannot survive a re-vendor and
# masquerade as bundled provenance. The snapshot is regenerated below.
for (p in c(upstream_root, dest_vendor, dest_cargo, vendor_config)) {
  if (file.exists(p) || dir.exists(p)) {
    unlink(p, recursive = TRUE, force = TRUE)
  }
}
dir.create(dest_upstream, recursive = TRUE, showWarnings = FALSE)

# ---- 3. materialize the snapshot via `git archive` ------------------------
#
# Whitelist: compile-time inputs only. Everything else upstream ships
# (tests/, benches/, examples/, comparison/, datasets/, scripts/,
# .github/, CHANGELOG, etc.) is excluded from the bundle. `build.rs` is
# REQUIRED: Cargo auto-detects and runs it even though it only emits a link
# directive under the (unused) `prima` feature; omitting it breaks the
# offline build. `LICENSE` is bundled so transitive license audits resolve.
# `docs/guide/` is REQUIRED: src/guide/mod.rs uses `include_str!` to embed
# its tutorial pages at compile time, so the .md files must be on disk.
tree_entries <- run_git(
  c("ls-tree", "--name-only", resolved_sha), dir = src_repo
)
required <- c("Cargo.toml", "Cargo.lock", "src", "build.rs")
missing <- setdiff(required, tree_entries)
if (length(missing)) {
  stop(sprintf(
    "[vendor-rust.R] pinned tree missing required entries: %s",
    paste(missing, collapse = ", ")
  ))
}
keep <- c(required, intersect(c("LICENSE", "README.md"), tree_entries))
if ("docs" %in% tree_entries) {
  keep <- c(keep, "docs/guide")
}

archive_tar <- file.path(tempdir(), "mixeff-rs-archive.tar")
unlink(archive_tar, force = TRUE)
res <- run_git(c(
  "archive", "--format=tar", "-o", archive_tar, resolved_sha, keep
), dir = src_repo)
if (!git_ok(res) || !file.exists(archive_tar)) {
  cat(res, sep = "\n")
  stop("[vendor-rust.R] git archive failed")
}
res <- system2(
  "tar", c("-xf", archive_tar, "-C", dest_upstream),
  stdout = TRUE, stderr = TRUE
)
if (!git_ok(res)) {
  cat(res, sep = "\n")
  stop("[vendor-rust.R] tar extract failed")
}
unlink(archive_tar, force = TRUE)

cat(sprintf(
  "[vendor-rust.R] bundled mixeff-rs @ %s (%d files)\n",
  substr(resolved_sha, 1, 12),
  length(list.files(dest_upstream, recursive = TRUE))
))

# ---- 4. write the provenance manifest -------------------------------------

manifest <- c(
  "# Generated by tools/vendor-rust.R -- do not edit.",
  "# Provenance for the bundled src/rust/upstream/mixeff-rs/ snapshot.",
  sprintf("crate        = \"mixeff-rs\""),
  sprintf("source       = \"%s\"", url),
  sprintf("requested    = \"%s\"", rev),
  sprintf("resolved_sha = \"%s\"", resolved_sha),
  sprintf("vendored_utc = \"%s\"",
          format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"))
)
writeLines(manifest, manifest_path)
cat(sprintf("[vendor-rust.R] wrote %s\n", manifest_path))

# ---- 5. verify Cargo.toml points at the bundled path ----------------------
#
# No silent surgery: this script does not rewrite the path dependency. It
# only asserts the committed Cargo.toml already targets the bundled snapshot.
cargo_toml_path <- file.path(src_rust, "Cargo.toml")
cargo_toml <- readLines(cargo_toml_path)
expected_path <- "upstream/mixeff-rs"
dep_line <- grep("^\\s*mixeff-rs\\s*=", cargo_toml, value = TRUE)
if (!length(dep_line)) {
  stop("[vendor-rust.R] no `mixeff-rs = ...` dependency line in Cargo.toml")
}
if (!grepl(expected_path, dep_line[1], fixed = TRUE)) {
  stop(sprintf(
    paste0("[vendor-rust.R] Cargo.toml mixeff-rs dep does not point at %s\n",
           "  found: %s\n",
           "  fix Cargo.toml manually (no silent rewrite)."),
    expected_path, trimws(dep_line[1])
  ))
}
cat(sprintf("[vendor-rust.R] Cargo.toml dep verified -> %s\n", expected_path))

# ---- 6. run cargo vendor --------------------------------------------------

cargo <- Sys.which("cargo")
if (!nzchar(cargo)) {
  cargo <- file.path(Sys.getenv("HOME"), ".cargo", "bin", "cargo")
}
if (!file.exists(cargo)) stop("[vendor-rust.R] cargo not found on PATH")

cat("[vendor-rust.R] running cargo vendor (this can take a minute)...\n")
res <- system2(
  command = cargo,
  args    = c("vendor", "--manifest-path", cargo_toml_path, dest_vendor),
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
# the canonical relative form, regardless of what cargo printed -- the
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

# ---- 7. write src/rust/vendor.tar.xz for tarball distribution -------------
#
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

# ---- 8. cleanup + summary -------------------------------------------------

if (!is.null(clone_dir) && dir.exists(clone_dir)) {
  unlink(clone_dir, recursive = TRUE, force = TRUE)
}

vendored_crates <- list.dirs(dest_vendor, recursive = FALSE, full.names = FALSE)
cat(sprintf(
  "[vendor-rust.R] done. mixeff-rs @ %s, %d vendored crates in src/vendor/\n",
  substr(resolved_sha, 1, 12),
  length(vendored_crates)
))
cat("[vendor-rust.R] next: devtools::check(document = FALSE)\n")
cat(paste0(
  "[vendor-rust.R] REMINDER: you just pulled new engine code. Consult\n",
  "                planning/upstream-blocked.md and unblock any downstream\n",
  "                beads whose upstream mixeff-rs feature has now shipped\n",
  "                (mote ls --tag upstream-blocked).\n"
))
