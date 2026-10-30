args <- commandArgs(TRUE)
builder <- args[[1L]]
tarball <- normalizePath("candidate/.release-gate/mixeff_0.2.0.tar.gz")
expected <- "722f88db7b71561870c7979a7e34e469f7fcbfbc3e51d14166fe3c89dad575c7"
stopifnot(unname(tools::sha256sum(tarball)) == expected)
receipt <- list(builder = builder, submitted_utc = format(Sys.time(), tz = "UTC"),
                tarball_sha256 = expected, bytes = file.size(tarball))
if (startsWith(builder, "win-")) {
  flavor <- sub("^win-", "", builder)
  destination <- paste0("ftp://win-builder.r-project.org/", flavor,
                        "/", basename(tarball))
  result <- devtools:::upload_ftp(tarball, destination)
  receipt$destination <- destination
  receipt$status_code <- result$status_code
} else {
  request <- httr2::request("https://mac.r-project.org/macbuilder/v1/submit") |>
    httr2::req_body_multipart(pkgfile = curl::form_file(tarball), rflavor = "r-devel") |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_timeout(120)
  result <- httr2::req_perform(request)
  receipt$status_code <- httr2::resp_status(result)
  receipt$response <- httr2::resp_body_json(result)
}
jsonlite::write_json(receipt, paste0(builder, "-submission.json"), auto_unbox = TRUE,
                     pretty = TRUE)
print(receipt)
