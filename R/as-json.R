#' Serialize a mixeff spec or fit to JSON
#'
#' `as_json()` returns a JSON string that carries the parsed object's public
#' R-side state and the raw compiler artifact JSON. `saveRDS()` / `readRDS()`
#' remains the primary persistence path for fitted objects; [revive()] restores
#' process-local caches after deserialization.
#'
#' @param x A compiled `mm_spec` or fitted `mm_fit`.
#' @param pretty Logical; pretty-print JSON when `TRUE`.
#' @param ... Reserved for future methods.
#'
#' @return A length-one character string containing JSON.
#'
#' @export
as_json <- function(x, pretty = FALSE, ...) {
  UseMethod("as_json")
}

#' @rdname as_json
#' @export
as_json.mm_compiled <- function(x, pretty = FALSE, ...) {
  artifact <- mm_compiled_artifact(x)
  payload <- list(
    schema = list(schema_name = "mixeff.r_object", schema_version = 1L),
    object_class = class(x),
    object_type = if (inherits(x, "mm_fit")) "fit" else "spec",
    formula = deparse1(x$formula),
    vars = x$vars %||% character(),
    artifact_json = attr(artifact, "raw_json") %||% jsonlite::toJSON(
      artifact,
      auto_unbox = TRUE,
      null = "null"
    )
  )
  if (inherits(x, "mm_fit")) {
    payload$fit <- list(
      REML = x$REML,
      beta = as.list(unname(x$beta)),
      beta_names = names(x$beta),
      theta = as.list(unname(x$theta)),
      sigma = x$sigma,
      logLik = x$logLik,
      deviance = x$deviance,
      df_residual = x$df_residual,
      fit_status = fit_status(x),
      nobs = x$nobs,
      fitted = as.list(unname(x$fitted)),
      residuals = as.list(unname(x$residuals))
    )
    payload$fit$schema <- x$schema %||% list()
  }
  jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", pretty = pretty)
}
