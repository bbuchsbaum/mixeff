#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

find_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        dir.exists(file.path(path, "R"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not locate mixeff repository root.", call. = FALSE)
    }
    path <- parent
  }
}

fixture_hash <- function(x, salt) {
  vapply(x, function(value) {
    tmp <- tempfile("mixeff-aphantasia-hash-")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(paste0(salt, ":", value), tmp, useBytes = TRUE)
    paste0("p_", substr(unname(tools::md5sum(tmp)), 1L, 16L))
  }, character(1), USE.NAMES = FALSE)
}

hash_participants <- function(ids, salt) {
  ids <- sort(unique(as.character(ids)))
  hashed <- fixture_hash(ids, salt)
  names(hashed) <- ids
  hashed
}

build_folder_map <- function(root) {
  files <- list.files(
    file.path(root, "data"),
    pattern = "Aphantasia_v.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  rows <- lapply(files, function(path) {
    pid <- tryCatch({
      dat <- read.csv(path, stringsAsFactors = FALSE)
      as.character(na.omit(dat$participant)[1L])
    }, error = function(cnd) NA_character_)
    rel <- sub(".*/data/", "", path)
    data.frame(
      participant = pid,
      source_folder = strsplit(rel, "/", fixed = TRUE)[[1L]][1L],
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows) |>
    filter(!is.na(participant)) |>
    distinct(participant, .keep_all = TRUE)
}

model_reference <- function(model, formula, model_type, family = NULL,
                            subset = NULL, source = "raw_trial_analysis.rds") {
  beta <- lme4::fixef(model)
  theta <- tryCatch(lme4::getME(model, "theta"), error = function(cnd) numeric())
  list(
    model_type = model_type,
    formula = formula,
    family = family,
    subset = subset,
    source = source,
    nobs = unname(stats::nobs(model)),
    logLik = unname(as.numeric(stats::logLik(model))),
    AIC = unname(stats::AIC(model)),
    BIC = unname(stats::BIC(model)),
    fixef = as.list(beta),
    theta = as.list(theta)
  )
}

fit_reference <- function(formula, data, model_type, family = NULL,
                          subset = NULL) {
  if (identical(model_type, "glmm")) {
    fit <- lme4::glmer(
      stats::as.formula(formula),
      data = data,
      family = stats::binomial(),
      control = glmm_control
    )
  } else {
    fit <- lme4::lmer(
      stats::as.formula(formula),
      data = data,
      REML = FALSE,
      control = lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 1e5)
      )
    )
  }
  model_reference(fit, formula, model_type, family, subset,
                  source = "tools/build-aphantasia-fixture.R")
}

read_cache <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("Required manuscript cache `%s` is missing: %s", label, path),
         call. = FALSE)
  }
  readRDS(path)
}

compact_data_frame <- function(x) {
  if (is.data.frame(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  x
}

main <- function() {
  repo_root <- find_repo_root()
  manuscript_root <- Sys.getenv(
    "MIXEFF_APHANTASIA_ROOT",
    "/Users/bbuchsbaum/Dropbox/manuscripts/Loo_aphantasia/revision3"
  )
  manuscript_root <- normalizePath(manuscript_root, winslash = "/", mustWork = TRUE)
  manuscript_dir <- file.path(manuscript_root, "manuscript")
  raw_script <- file.path(manuscript_dir, "raw_trial_analysis.R")
  if (!file.exists(raw_script)) {
    stop("Could not find manuscript raw_trial_analysis.R.", call. = FALSE)
  }

  old_wd <- setwd(manuscript_dir)
  on.exit(setwd(old_wd), add = TRUE)
  source(raw_script, chdir = TRUE)

  cache_dir <- file.path(manuscript_dir, ".cache")
  analysis <- read_cache(file.path(cache_dir, "raw_trial_analysis.rds"),
                         "raw_trial_analysis.rds")
  supp <- read_cache(file.path(cache_dir, "supplemental.rds"),
                     "supplemental.rds")
  supp_s3 <- read_cache(file.path(cache_dir, "supp_S3.rds"),
                        "supp_S3.rds")
  supp_s9 <- read_cache(file.path(cache_dir, "supp_S9.rds"),
                        "supp_S9.rds")

  loaded <- load_trial_data(manuscript_root)
  folder_map <- build_folder_map(manuscript_root)
  salt <- "mixeff-aphantasia-revision3-v1"
  id_map <- hash_participants(loaded$trials$participant, salt)

  trials <- loaded$trials |>
    left_join(folder_map, by = "participant") |>
    mutate(participant = unname(id_map[participant])) |>
    transmute(
      participant,
      bubbled,
      back_masked,
      SOA,
      block_num,
      trial_image,
      category,
      correct,
      rt,
      aphantasia,
      age,
      vviq_standard,
      source,
      source_folder
    )

  metadata <- loaded$meta |>
    left_join(folder_map, by = "participant") |>
    mutate(
      participant = unname(id_map[participant]),
      group = ifelse(aphantasia == "yes", "aphant", "control")
    ) |>
    transmute(
      participant,
      group,
      age,
      vviq_standard,
      source,
      source_folder
    ) |>
    distinct(participant, .keep_all = TRUE)

  exclude_hashed <- unname(id_map[analysis$exclude_ids])

  primary_formula <- paste(
    "correct ~ group * mask * soa_s + block +",
    "(1 + mask + soa_s || participant) + (1 | item)"
  )
  combined_formula <- paste(
    "correct ~ group * mask * soa_s * stimtype + block +",
    "(1 + mask + soa_s || participant) + (1 | item)"
  )
  rt_formula <- paste(
    "log_rt ~ group * mask * soa_s + block +",
    "(1 | participant) + (1 | item)"
  )

  dat_primary <- analysis$primary$data
  dat_age <- dat_primary |> filter(!is.na(age))
  dat_folder <- dat_primary |>
    mutate(participant_raw = as.character(participant)) |>
    left_join(folder_map, by = c("participant_raw" = "participant"))
  age_controls <- dat_folder |>
    filter(aphantasia != "yes",
           source_folder == "prolific_control_age_match") |>
    distinct(participant) |>
    pull(participant)
  aphants <- dat_folder |>
    filter(aphantasia == "yes") |>
    distinct(participant) |>
    pull(participant)
  dat_matched <- dat_folder |>
    filter(participant %in% c(aphants, age_controls)) |>
    select(-participant_raw, -source_folder)
  dat_age$age_z <- as.numeric(scale(dat_age$age))
  dat_matched_age <- dat_matched
  dat_matched_age$age_z <- as.numeric(scale(dat_matched_age$age))

  s1_specs <- list(
    "S1_intercept_only" = paste(
      "correct ~ group * mask * soa_s + block +",
      "(1 | participant) + (1 | item)"
    ),
    "S1_current_uncorrelated_slopes" = primary_formula,
    "S1_correlated_slopes" = paste(
      "correct ~ group * mask * soa_s + block +",
      "(1 + mask + soa_s | participant) + (1 | item)"
    ),
    "S1_item_mask_slope" = paste(
      "correct ~ group * mask * soa_s + block +",
      "(1 + mask + soa_s || participant) + (1 + mask | item)"
    ),
    "S1_maximal" = paste(
      "correct ~ group * mask * soa_s + block +",
      "(1 + mask * soa_s | participant) + (1 + group | item)"
    )
  )

  models <- list(
    primary = model_reference(
      analysis$primary$model,
      primary_formula,
      "glmm",
      family = "binomial(logit)",
      subset = "bubbled == yes, non-missing correct, intermediate VVIQ controls excluded"
    ),
    sensitivity = model_reference(
      analysis$sensitivity$model,
      primary_formula,
      "glmm",
      family = "binomial(logit)",
      subset = "bubbled == yes, non-missing correct, intermediate VVIQ controls assigned to control"
    ),
    intact = model_reference(
      analysis$intact$model,
      primary_formula,
      "glmm",
      family = "binomial(logit)",
      subset = "bubbled == no, non-missing correct, intermediate VVIQ controls excluded"
    ),
    combined = model_reference(
      analysis$combined$model,
      combined_formula,
      "glmm",
      family = "binomial(logit)",
      subset = "all stimuli, non-missing correct, intermediate VVIQ controls excluded"
    ),
    rt = model_reference(
      analysis$rt$model,
      rt_formula,
      "lmm",
      family = "gaussian",
      subset = "primary occluded correct trials with finite positive RT"
    )
  )

  for (id in names(s1_specs)) {
    models[[id]] <- fit_reference(
      s1_specs[[id]],
      dat_primary,
      "glmm",
      family = "binomial(logit)",
      subset = id
    )
  }
  models$S7_age_covariate <- fit_reference(
    paste(primary_formula, "+ age_z"),
    dat_age,
    "glmm",
    family = "binomial(logit)",
    subset = "primary occluded trials with non-missing age plus age_z"
  )
  models$S9_age_matched_subset <- fit_reference(
    primary_formula,
    dat_matched,
    "glmm",
    family = "binomial(logit)",
    subset = "aphantasia plus prolific_control_age_match controls"
  )
  models$S9_age_matched_subset_age_covariate <- fit_reference(
    paste(primary_formula, "+ age_z"),
    dat_matched_age,
    "glmm",
    family = "binomial(logit)",
    subset = "age-matched subset plus age_z"
  )

  reference <- list(
    schema = list(
      name = "mixeff.aphantasia_fixture_reference",
      version = "1"
    ),
    provenance = list(
      manuscript_root = manuscript_root,
      raw_cache = ".cache/raw_trial_analysis.rds",
      supplemental_cache = ".cache/supplemental.rds",
      supplemental_s3_cache = ".cache/supp_S3.rds",
      supplemental_s9_cache = ".cache/supp_S9.rds",
      anonymization = list(
        participant_hash = "md5",
        salt_label = salt,
        hash_prefix = "p_",
        hash_chars = 16L
      )
    ),
    excluded_participants = as.list(exclude_hashed),
    counts = list(
      trials = nrow(trials),
      participants = length(unique(trials$participant)),
      metadata_rows = nrow(metadata),
      primary_trials = nrow(analysis$primary$data),
      sensitivity_trials = nrow(analysis$sensitivity$data),
      intact_trials = nrow(analysis$intact$data),
      combined_trials = nrow(analysis$combined$data),
      rt_trials = nrow(analysis$rt$data)
    ),
    tolerances = list(
      glmm = list(fixef_abs = 2.5e-2, logLik_rel = 1e-3, AIC_rel = 1e-3),
      lmm = list(fixef_abs = 1e-4, logLik_rel = 1e-5, AIC_rel = 1e-5)
    ),
    sample = compact_data_frame(analysis$sample),
    quality = compact_data_frame(analysis$quality),
    models = models,
    inference = list(
      primary_dd = compact_data_frame(analysis$primary$dd),
      sensitivity_dd = compact_data_frame(analysis$sensitivity$dd),
      intact_dd = compact_data_frame(analysis$intact$dd),
      combined_did_diff = compact_data_frame(analysis$combined$did_diff),
      S1 = compact_data_frame(supp$S1$result),
      S3 = compact_data_frame(supp_s3$S3$result$summary),
      S6 = compact_data_frame(supp$S6$result),
      S7 = compact_data_frame(supp$S7$result),
      S9 = compact_data_frame(supp_s9$S9$result)
    )
  )

  out_dir <- file.path(repo_root, "inst", "extdata", "aphantasia")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(trials, file.path(out_dir, "trials.rds"), version = 2)
  saveRDS(metadata, file.path(out_dir, "metadata.rds"), version = 2)
  writeLines(
    jsonlite::toJSON(reference, pretty = TRUE, auto_unbox = TRUE,
                     null = "null", digits = 15),
    file.path(out_dir, "reference.json")
  )

  message("Wrote aphantasia fixture to ", out_dir)
  invisible(out_dir)
}

main()
