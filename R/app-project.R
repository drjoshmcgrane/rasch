# rasch :: Shiny project files

# The graphical interface stores complete analysis sessions as ordinary RDS
# files with a small identifying header. Validation lives here so it can be
# exercised without a running Shiny session. Schema 2 adds an integrity seal:
# schema 1 stored data and fits side by side but had no way to establish that
# they came from the same saved analysis.

.app_fit_family <- function(fit) {
  if (inherits(fit, "rasch_btl")) return("btl")
  if (inherits(fit, "rasch_efrm")) return("efrm")
  if (inherits(fit, "rasch_mfrm")) return("mfrm")
  if (inherits(fit, "rasch")) return("rasch")
  NA_character_
}

.app_scalar_text <- function(x) {
  is.character(x) && is.null(dim(x)) && is.null(oldClass(x)) &&
    length(x) == 1L && !is.na(x) && nzchar(trimws(x))
}

.validate_app_frame_calibration <- function(fit, label) {
  if (!inherits(fit, c("rasch_efrm", "rasch_btl_efrm")))
    return(invisible(NULL))
  if (!identical(fit$calibration_algorithm, "frame-likelihood-1"))
    stop(label, " ", paste("has an unverified frame calibration; refit this analysis",
      "with the current likelihood checks before reopening it.",
      "The saved file is unchanged. Recover its source data and settings",
      "with readRDS(file)$data and readRDS(file)$settings."), call. = FALSE)
  # An EFRM also stores score curves, whose designs and labels now come from
  # the shared design enumeration. The guard is on the EFRM class alone:
  # btl_efrm() shares the calibration tag but has no score curves, and must
  # keep reading.
  if (inherits(fit, "rasch_efrm") &&
      !identical(fit$efrm_results_algorithm, "efrm-results-2"))
    stop(label, " ", paste("holds score curves from a superseded design",
      "enumeration; refit this analysis before reopening it.",
      "The saved file is unchanged. Recover its source data and settings",
      "with readRDS(file)$data and readRDS(file)$settings."), call. = FALSE)
  invisible(NULL)
}

# Derived results can contain their own frame fits, outside base/history/kept
# slots. Walk the authenticated bundle so those fits cannot bypass the guard.
.validate_app_frame_fits <- function(x, label = "the saved analysis") {
  if (inherits(x, c("rasch", "rasch_btl"))) {
    .validate_app_frame_calibration(x, label)
    return(invisible(NULL))
  }
  if (is.list(x) && !is.data.frame(x)) for (i in seq_along(x)) {
    nm <- names(x)[i]
    part <- if (length(nm) && !is.na(nm) && nzchar(nm)) paste0("$", nm)
            else paste0("[[", i, "]]")
    .validate_app_frame_fits(x[[i]], paste0(label, part))
  }
  invisible(NULL)
}

.validate_app_fit <- function(fit, what = "fitted model", check_algorithm = TRUE) {
  fail <- function(message) stop(what, " ", message, call. = FALSE)
  family <- .app_fit_family(fit)
  if (is.na(family)) fail("is not a fitted rasch model")
  if (check_algorithm) .validate_app_frame_calibration(fit, what)

  if (identical(family, "btl")) {
    if (!is.data.frame(fit$objects) || nrow(fit$objects) < 2L ||
        !all(c("object", "location") %in% names(fit$objects)) ||
        !is.character(fit$objects$object) || anyNA(fit$objects$object) ||
        any(!nzchar(trimws(fit$objects$object))) ||
        anyDuplicated(fit$objects$object))
      fail("has an invalid object calibration")
    if (!is.data.frame(fit$comparisons) || !nrow(fit$comparisons) ||
        !all(c("object_a", "object_b", "response", "weight", "judge") %in%
             names(fit$comparisons)))
      fail("has an invalid comparison design")
    if (!is.data.frame(fit$pairs) ||
        !all(c("object_a", "object_b") %in% names(fit$pairs)))
      fail("has an invalid pair-fit table")
    if (length(fit$m) != 1L || !is.numeric(fit$m) || is.complex(fit$m) ||
        !is.null(dim(fit$m)) || !is.null(oldClass(fit$m)) ||
        !is.finite(fit$m) ||
        fit$m < 1 || fit$m != floor(fit$m))
      fail("has an invalid response scale")
    if (length(fit$n_comparisons) != 1L ||
        !is.numeric(fit$n_comparisons) || is.complex(fit$n_comparisons) ||
        !is.null(dim(fit$n_comparisons)) ||
        !is.null(oldClass(fit$n_comparisons)) ||
        !is.finite(fit$n_comparisons) ||
        fit$n_comparisons < 1)
      fail("has an invalid comparison count")
    if (length(fit$converged) != 1L || !is.logical(fit$converged) ||
        !is.null(dim(fit$converged)) || !is.null(oldClass(fit$converged)) ||
        is.na(fit$converged))
      fail("has no valid convergence record")

    # The extended-frame result has its own fitted-probability machinery.
    # Ordinary and explanatory BTL fits must retain the row-wise probabilities
    # and public refit specification used by the fitted-model bootstrap.
    if (inherits(fit, "rasch_btl_efrm")) {
      if (!is.data.frame(fit$phi_table) || !is.data.frame(fit$alpha_table) ||
          !is.data.frame(fit$kappa_table) || !is.data.frame(fit$frames))
        fail("has an incomplete frame calibration")
    } else {
      if (!is.matrix(fit$fitted_prob) ||
          !identical(dim(fit$fitted_prob),
                     c(nrow(fit$comparisons), as.integer(fit$m) + 1L)) ||
          any(!is.finite(fit$fitted_prob)))
        fail("has invalid fitted comparison probabilities")
      if (!is.list(fit$refit_spec))
        fail("has no valid refit specification")
    }
    return(invisible(family))
  }

  X <- fit$X
  if (!is.matrix(X) || length(dim(X)) != 2L || any(dim(X) < 1L))
    fail("has no valid response matrix")
  N <- nrow(X); L <- ncol(X)
  if (!.app_scalar_text(fit$model)) fail("has no valid model name")
  if (!is.numeric(fit$m) || length(fit$m) != L ||
      any(!is.finite(fit$m)) || any(fit$m < 1) || any(fit$m != floor(fit$m)))
    fail("has an invalid item response scale")
  if (!is.data.frame(fit$items) || nrow(fit$items) != L ||
      !all(c("item", "location") %in% names(fit$items)) ||
      !is.character(fit$items$item) || anyNA(fit$items$item) ||
      any(!nzchar(trimws(fit$items$item))) || anyDuplicated(fit$items$item))
    fail("has an invalid item calibration")
  if (!is.data.frame(fit$thresholds) ||
      !all(c("item", "k", "tau") %in% names(fit$thresholds)) ||
      nrow(fit$thresholds) != sum(fit$m))
    fail("has an invalid threshold calibration")
  if (!is.list(fit$tau_list) || length(fit$tau_list) != L ||
      !identical(as.integer(lengths(fit$tau_list)), as.integer(fit$m)) ||
      any(!vapply(fit$tau_list, function(z)
        is.numeric(z) && all(is.finite(z)), logical(1))))
    fail("has invalid item threshold vectors")
  if (!is.data.frame(fit$person) || nrow(fit$person) != N ||
      !all(c("id", "raw", "theta") %in% names(fit$person)))
    fail("has an invalid person calibration")
  if (!is.matrix(fit$residuals) || !identical(dim(fit$residuals), dim(X)))
    fail("has an invalid residual matrix")
  if (!is.list(fit$moments) || !is.matrix(fit$moments$E) ||
      !is.matrix(fit$moments$V) ||
      !identical(dim(fit$moments$E), dim(X)) ||
      !identical(dim(fit$moments$V), dim(X)))
    fail("has invalid fitted response moments")
  if (!is.list(fit$est) || length(fit$est$converged) != 1L ||
      !is.logical(fit$est$converged) || is.na(fit$est$converged))
    fail("has no valid estimation record")
  if (!is.null(fit$factors) &&
      (!is.data.frame(fit$factors) || nrow(fit$factors) != N))
    fail("has an invalid person-factor design")
  invisible(family)
}

.app_project_binding <- function(project) {
  x <- project
  attr(x, "rasch_project_legacy") <- NULL
  attr(x, "rasch_project_legacy_dropped") <- NULL
  x$binding <- NULL
  .fit_boot_md5(x)
}

.seal_app_project <- function(project) {
  if (!is.list(project)) stop("`project` must be a list", call. = FALSE)
  project$schema <- 2L
  project$binding <- .app_project_binding(project)
  project
}

# Older files did not identify their person-scoring algorithm. Check them
# numerically without altering the saved fit or breaking signed downstream
# results. A changed score requires a refit of all dependent diagnostics, not
# just replacement of the person table. Current fits carry the solver stamp.
.validate_app_person_scoring <- function(fit, label) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    return(invisible(NULL))
  # Scoring an old PC calibration again cannot repair a different threshold
  # model. Fits that requested kurtosis need the corrected calibration too.
  if (isTRUE(fit$est$n_components == 4L) && any(fit$m >= 4L) &&
      !identical(fit$est$pc_algorithm, "guttman-four-1"))
    stop(label, paste("predates the corrected principal-component kurtosis;",
      "refit this analysis before reopening it. The saved file is unchanged;",
      "recover its source data with readRDS(file)$data."), call. = FALSE)
  stamp <- fit$person_scoring_algorithm
  if (identical(stamp, "max-wle-1")) return(invisible(NULL))
  fail <- function() stop(label, paste(
    "uses superseded or unverified person scoring; refit this analysis",
    "before reopening it. The saved file has not been changed.",
    "Its source data can be recovered in R with readRDS(file)$data."),
    call. = FALSE)
  if (!is.null(stamp)) fail()
  disc <- fit$disc %||% rep(1, ncol(fit$X))
  if (length(disc) == 1L) disc <- rep(disc, ncol(fit$X))
  if (!is.numeric(disc) || length(disc) != ncol(fit$X) ||
      any(!is.finite(disc) | disc <= 0)) fail()
  unit <- max(disc)
  same <- function(a, b) {
    if (!is.numeric(a) || !is.numeric(b) || length(a) != length(b) ||
        !identical(is.na(a), is.na(b))) return(FALSE)
    ok <- !is.na(b)
    a <- a[ok] * unit; b <- b[ok] * unit
    finite <- is.finite(b)
    all(is.finite(a) == finite) &&
      all(a[!finite] == b[!finite]) &&
      all(abs(a[finite] - b[finite]) <= 1e-7 * pmax(1, abs(b[finite])))
  }
  equal <- length(unique(disc)) == 1L
  expected <- tryCatch(if (equal)
    .person_estimates(fit$X, fit$tau_list, disc = disc[1L]) else
    .efrm_person_estimates(fit$X, fit$tau_list, disc),
    error = function(e) NULL)
  if (is.null(expected)) fail()
  if (!isTRUE(fit$est$converged)) expected$se[] <- NA_real_
  if (!same(fit$person$theta, expected$theta) ||
      !same(fit$person$se, expected$se)) fail()
  if (!is.null(fit$score_table)) {
    if (!equal) fail()
    sc <- person_wle(fit$tau_list, disc = disc[1L])
    if (!isTRUE(fit$est$converged)) sc$se[] <- NA_real_
    if (!identical(as.numeric(fit$score_table$score),
                   as.numeric(0:sum(fit$m))) ||
        !same(fit$score_table$theta, unname(sc$theta)) ||
        !same(fit$score_table$se, unname(sc$se))) fail()
  }
  invisible(NULL)
}

# One statement of what an embedded display resource must satisfy, shared by
# the loader and by the app that writes the file. The app already refuses an
# unusable equating reference or panel map on its own page; without the same
# test at the writing end an upload the app tolerates would make the whole
# analysis unsaveable, and the failure would surface only as a broken
# download. `m` is the fitted response-scale size a polytomous comparison
# bank must carry, and NULL wherever that check does not apply.
.app_project_resource_problem <- function(name, value, m = NULL) {
  if (is.null(value)) return(NULL)
  switch(name,
    eq_reference = tryCatch({
      ref <- .equate_ref(value)
      .equate_bank_cov(value, ref$item)
      NULL
    }, error = function(e) conditionMessage(e)),
    bt_eq_bank = tryCatch({
      ref <- .btl_equate_ref(value)
      .btl_equate_bank_cov(value, ref$object)
      .btl_equate_cov_df(value)
      if (!is.null(m)) {
        bm <- attr(value, "m", exact = TRUE)
        if (!is.numeric(bm) || is.complex(bm) || length(bm) != 1L ||
            !is.null(dim(bm)) || !is.null(oldClass(bm)) || !is.finite(bm) ||
            bm < 1L || bm > .Machine$integer.max || bm != floor(bm) ||
            !identical(as.integer(bm), as.integer(m)))
          stop("its response-scale metadata do not match the saved fit",
               call. = FALSE)
      }
      NULL
    }, error = function(e) conditionMessage(e)),
    wright_item_map = {
      z <- value
      ok <- is.data.frame(z) && !anyDuplicated(names(z)) &&
        all(c("item", "panel") %in% names(z)) && nrow(z) > 0L &&
        !anyNA(z$item) && !anyNA(z$panel) &&
        all(nzchar(trimws(as.character(z$item)))) &&
        all(nzchar(trimws(as.character(z$panel)))) &&
        !anyDuplicated(trimws(as.character(z$item)))
      if (ok) NULL else "it is not a usable item panel map"
    },
    NULL)
}

.validate_app_project <- function(project) {
  fail <- function(message) stop(message, call. = FALSE)
  if (!is.list(project) || !identical(project$format, "rasch-shiny-project"))
    fail("not a rasch analysis file")
  # Read the schema as stored: coercion would accept "2", TRUE, 2.5 or a
  # factor's level code as schema 2.
  schema <- project$schema
  if (length(schema) != 1L || !is.numeric(schema) || is.complex(schema) ||
      !is.null(dim(schema)) || !is.null(oldClass(schema)) || is.na(schema) ||
      schema != 2L)
    fail(paste("unsupported rasch analysis-file schema; this version needs",
               "a schema-2 file with data-to-fit integrity information"))
  if (!(is.data.frame(project$data) || is.matrix(project$data)) ||
      nrow(project$data) < 1L || ncol(project$data) < 1L)
    fail("the analysis file does not contain a valid source dataset")
  data_names <- colnames(project$data)
  if (is.null(data_names) || anyNA(data_names) ||
      any(!nzchar(trimws(data_names))) || anyDuplicated(data_names))
    fail("the analysis file has invalid source-data column names")

  base_problem <- tryCatch({
    .validate_app_fit(project$base_fit, "the saved base fit", check_algorithm = FALSE)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(base_problem)) fail(base_problem)
  base_family <- .app_fit_family(project$base_fit)
  if (!is.null(project$model_type)) {
    # A factor would pass a membership test and then reach switch() as its
    # integer code, selecting a branch by level order.
    if (!.app_scalar_text(project$model_type) ||
        !(project$model_type %in% c("rasch", "mfrm", "efrm", "btl")))
      fail("the analysis file names an unsupported model type")
    if (!identical(base_family, project$model_type))
      fail(sprintf(paste("the analysis file declares model type '%s' but",
                         "contains a fit of class '%s'"),
                   project$model_type, class(project$base_fit)[1]))
  }
  for (field in c("rasch_steps", "btl_steps", "kept_fits", "kept_fit_code",
                  "simulation", "results", "settings", "resources"))
    if (!is.null(project[[field]]) && !is.list(project[[field]]))
      fail(sprintf("the analysis file has an invalid %s field", field))

  # File-backed display analyses are embedded in the project because their
  # upload paths expire with the Shiny session. Validate the parsed resources
  # before any reactive equating or Wright-map code receives them.
  resources <- project$resources %||% list()
  if (length(resources)) {
    resource_names <- names(resources)
    if (is.null(resource_names) || anyNA(resource_names) ||
        any(!nzchar(trimws(resource_names))) || anyDuplicated(resource_names))
      fail("the analysis file has invalid resource names")
  }
  problem <- .app_project_resource_problem("eq_reference",
                                           resources[["eq_reference"]])
  if (!is.null(problem))
    fail(paste("the analysis file has an invalid item-equating reference:",
               problem))
  problem <- .app_project_resource_problem(
    "bt_eq_bank", resources[["bt_eq_bank"]],
    m = if (identical(base_family, "btl") && isTRUE(project$base_fit$m > 1L))
      project$base_fit$m else NULL)
  if (!is.null(problem))
    fail(paste("the analysis file has an invalid object-equating bank:",
               problem))
  if (!is.null(.app_project_resource_problem("wright_item_map",
                                             resources[["wright_item_map"]])))
    fail("the analysis file has an invalid Wright-map item panel map")

  # Current app fits retain the exact data, controls and uploaded metadata
  # used for their base calibration. The enclosing project must reproduce
  # that source. Older projects have no such attribute and remain readable.
  fit_source <- attr(project$base_fit, "rasch_app_source", exact = TRUE)
  if (!is.null(fit_source)) {
    source_ok <- is.list(fit_source) && is.data.frame(fit_source$data) &&
      nrow(fit_source$data) > 0L && ncol(fit_source$data) > 0L &&
      is.list(fit_source$settings) && is.list(fit_source$resources) &&
      is.list(fit_source$simulation)
    if (!source_ok)
      fail("the saved base fit has invalid app run metadata")
    source_names <- colnames(fit_source$data)
    if (is.null(source_names) || anyNA(source_names) ||
        any(!nzchar(trimws(source_names))) || anyDuplicated(source_names))
      fail("the saved base fit has invalid source-data column names")
    if (!identical(fit_source$data,
                   as.data.frame(project$data, check.names = FALSE)))
      fail("the saved base fit does not belong to the source dataset")

    validate_source_names <- function(x, what) {
      if (!length(x)) return(invisible(NULL))
      nm <- names(x)
      if (is.null(nm) || anyNA(nm) || any(!nzchar(trimws(nm))) ||
          anyDuplicated(nm))
        fail(sprintf("the saved base fit has invalid %s names", what))
      invisible(NULL)
    }
    validate_source_names(fit_source$settings, "run-setting")
    validate_source_names(fit_source$resources, "run-resource")
    validate_source_names(fit_source$simulation, "simulation-metadata")

    for (nm in names(fit_source$settings))
      if (!nm %in% names(project$settings) ||
          !identical(project$settings[[nm]], fit_source$settings[[nm]]))
        fail(sprintf(paste("the saved base fit's `%s` setting does not",
                           "match the analysis file"), nm))
    for (nm in names(fit_source$resources))
      if (!nm %in% names(project$resources) ||
          !identical(project$resources[[nm]], fit_source$resources[[nm]]))
        fail(sprintf(paste("the saved base fit's `%s` resource does not",
                           "match the analysis file"), nm))
    if (!identical(project$simulation %||% list(), fit_source$simulation))
      fail(paste("the saved base fit's simulation metadata does not match",
                 "the analysis file"))
  }

  # Kept fits are used directly by comparison and equating after a project is
  # reopened. Validate each one here rather than allowing a malformed entry to
  # fail later inside a table or plot. Names are part of the selector state and
  # therefore must be stable and unambiguous.
  kept <- project$kept_fits %||% list()
  if (length(kept)) {
    kept_names <- names(kept)
    if (is.null(kept_names) || anyNA(kept_names) ||
        any(!nzchar(trimws(kept_names))) || anyDuplicated(kept_names))
      fail("the analysis file has invalid kept-fit names")
    for (i in seq_along(kept)) {
      problem <- tryCatch({
        .validate_app_fit(kept[[i]],
                          sprintf("the kept fit '%s'", kept_names[i]),
                          check_algorithm = FALSE)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem)) fail(problem)
    }
  }

  validate_history <- function(history, family, field) {
    if (!length(history)) return(invisible(NULL))
    for (i in seq_along(history)) {
      entry <- history[[i]]
      metadata_ok <- is.list(entry) && .app_scalar_text(entry$type) &&
        .app_scalar_text(entry$label) && is.list(entry$details) &&
        (is.null(entry$code) || .app_scalar_text(entry$code)) &&
        .app_scalar_text(entry$created)
      if (!metadata_ok)
        fail(sprintf("the analysis file has invalid %s history metadata at entry %d",
                     field, i))
      problem <- tryCatch({
        .validate_app_fit(entry$fit,
                          sprintf("the %s history fit at entry %d", field, i),
                          check_algorithm = FALSE)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem) ||
          !identical(.app_fit_family(entry$fit), family))
        fail(sprintf("the analysis file has an invalid fitted-model history in %s at entry %d%s",
                     field, i,
                     if (is.null(problem)) "" else paste0(": ", problem)))
    }
    invisible(NULL)
  }

  is_btl <- identical(base_family, "btl")
  rasch_history <- project$rasch_steps %||% list()
  btl_history <- project$btl_steps %||% list()
  if ((is_btl && length(rasch_history)) ||
      (!is_btl && length(btl_history)))
    fail("the analysis file has fitted-model history for the inactive model family")
  validate_history(if (is_btl) btl_history else rasch_history,
                   base_family, if (is_btl) "btl_steps" else "rasch_steps")
  history <- if (is_btl) btl_history else rasch_history
  active_fit <- if (length(history)) history[[length(history)]]$fit
                else project$base_fit

  # Externally weighted person estimates are derived from the active
  # calibration. Authenticate both the stored table and that relationship;
  # otherwise a table calculated before a split, superitem or other refit can
  # be reopened beside a different calibration while still looking plausible.
  person_weights <- project$results[["person_weights"]]
  if (!is.null(person_weights)) {
    problem <- tryCatch({
      .validate_weighted_person_result(person_weights, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the analysis file has invalid saved weighted person estimates:",
                 problem))
  }

  guessing <- project$results[["guessing"]]
  if (!is.null(guessing)) {
    problem <- tryCatch({
      .validate_tailored_result(guessing, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved tailored analysis does not belong to the active fit:",
                 problem))
  }

  # Results belong to the active fit, which is the final structural change.
  # In particular, a bootstrap null from an earlier fit must not be restored
  # beside later DIF splits, superitems or paired-comparison frame changes.
  bootstrap <- project$results$bootstrap
  if (!is.null(bootstrap)) {
    if (!is.list(bootstrap) || is.null(bootstrap$bs))
      fail("the analysis file has an invalid saved bootstrap result")
    problem <- tryCatch({
      .validate_fit_bootstrap(bootstrap$bs, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved bootstrap result does not belong to the active fit:",
                 problem))
  }

  # Exact indexing matters here: `$dif` partially matches `dif_bootstrap`
  # when the primary result is absent. A primary result is tied to the active
  # fit even when no bootstrap was requested; the project seal proves only
  # that the bundle has not changed since it was written, not that two
  # separately supplied fitted objects belong together.
  primary_dif <- if (is_btl) project$results[["btl_dif"]] else
    project$results[["dif"]]
  if (!is.null(primary_dif)) {
    problem <- tryCatch({
      if (is_btl) .validate_btl_dif_result(primary_dif, active_fit)
      else .validate_dif_result(primary_dif, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved primary DIF analysis does not belong to the active fit:",
                 problem))
    # Ordinary Rasch DIF is recomputed reactively after the controls are
    # restored. Pin the stored analysis to those controls even when it has no
    # bootstrap, so the displayed run is the one that was validated.
    if (!is_btl) {
      settings <- project$settings %||% list()
      effects <- settings$dif_effects %||% "main"
      alpha <- settings$dif_alpha %||% 0.05
      if (!.app_scalar_text(effects) ||
          !identical(primary_dif$effects, effects) ||
          length(alpha) != 1L || !is.numeric(alpha) || !is.finite(alpha) ||
          !isTRUE(all.equal(primary_dif$alpha, alpha, tolerance = 0)))
        fail(paste("the saved primary DIF analysis does not match the",
                   "restored DIF settings"))
    }
  }

  btl_dif_meta <- project$results[["btl_dif_meta"]]
  if (!is.null(btl_dif_meta)) {
    if (!is_btl || is.null(primary_dif) || !is.list(btl_dif_meta) ||
        !.app_scalar_text(btl_dif_meta$judge_col) ||
        !btl_dif_meta$judge_col %in% data_names)
      fail("the analysis file has invalid saved Comparative Judgement DIF display metadata")
  }
  if (is_btl && !is.null(primary_dif)) {
    # A current CJ DIF result is keyed to the judge role used by the fit, not
    # merely to the column selected for display. Without the source binding,
    # a project can be signed while its maps silently refer to another role.
    source <- attr(project$base_fit, "rasch_app_source", exact = TRUE)
    role_ok <- is.list(source) && is.data.frame(source$data) &&
      is.list(source$settings) && .app_scalar_text(source$settings$bt_judge) &&
      is.list(btl_dif_meta) &&
      identical(btl_dif_meta$fitted_judge_col, source$settings$bt_judge) &&
      .app_scalar_text(btl_dif_meta$judge_col) &&
      btl_dif_meta$fitted_judge_col %in% names(source$data) &&
      btl_dif_meta$judge_col %in% data_names &&
      identical(as.character(source$data[[btl_dif_meta$fitted_judge_col]]),
                as.character(project$data[[btl_dif_meta$judge_col]]))
    if (!isTRUE(role_ok))
      fail(paste("the saved Comparative Judgement DIF lacks authenticated",
                 "fitted judge-role provenance"))
  }

  dif_bootstrap <- project$results$dif_bootstrap
  if (!is.null(dif_bootstrap)) {
    if (!is.list(dif_bootstrap) || is.null(dif_bootstrap$db))
      fail("the analysis file has an invalid saved DIF bootstrap result")
    if (is.null(primary_dif))
      fail(paste("the saved DIF bootstrap has no accompanying primary DIF",
                 "analysis"))
    problem <- tryCatch({
      .validate_dif_bootstrap(dif_bootstrap$db, active_fit, primary_dif)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved DIF bootstrap result does not belong to the active fit and DIF analysis:",
                 problem))
  }

  contrasts <- project$results[["contrasts"]]
  if (!is.null(contrasts) &&
      (is_btl || !inherits(contrasts, "rasch_dif_contrasts") ||
       !identical(contrasts$algorithm, "complete-contrast-cells-2")))
    fail("the saved planned DIF contrasts use a superseded calculation; recompute them")

  resolution <- project$results[["resolve"]]
  if (!is.null(resolution) &&
      (is_btl || !inherits(resolution, "rasch_resolve_dif") ||
       !identical(resolution$algorithm, "factor-design-resolution-2") ||
       !.app_scalar_text(resolution$effects) ||
       !resolution$effects %in% c("main", "factorial")))
    fail(paste("the saved automatic DIF resolution uses a superseded",
               "calculation; recompute it"))

  dimensionality <- project$results$dimensionality
  if (!is.null(dimensionality)) {
    problem <- tryCatch({
      if (is_btl) .validate_btl_dimensionality(dimensionality, active_fit)
      else .validate_scree_result(dimensionality, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved dimensionality analysis does not belong to the active fit:",
                 problem))
  }
  subtest <- project$results$subtest
  if (!is.null(subtest)) {
    if (is_btl)
      fail("the saved person-subset dimensionality test accompanies a paired-comparison fit")
    if (!is.list(subtest) ||
        !identical(subtest$algorithm, "person-subset-comparison-1"))
      fail(paste("the saved person-subset dimensionality test uses a superseded",
                 "person-subset comparison; recompute it"))
    problem <- tryCatch({
      .validate_dimensionality_test(subtest, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved person-subset dimensionality test does not belong to the active fit:",
                 problem))
  }
  invariance <- project$results$frame_invariance
  if (!is.null(invariance)) {
    problem <- tryCatch({
      .validate_frame_invariance(invariance, active_fit)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(problem))
      fail(paste("the saved frame-invariance analysis does not belong to the active fit:",
                 problem))
  }

  if (!.app_scalar_text(project$binding))
    fail("the analysis file has no valid data-to-fit integrity information")
  unsigned_project <- project
  attr(unsigned_project, "rasch_project_legacy") <- NULL
  attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
  unsigned_project$binding <- NULL
  if (!.fit_boot_hash_matches(project$binding, unsigned_project))
    fail(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"))
  # Check only after structural validation and authentication. Inspect every
  # retained fit, including inactive history and comparison/equating fits.
  .validate_app_frame_fits(project)
  .validate_app_person_scoring(project$base_fit, "the saved base fit")
  for (field in c("rasch_steps", "btl_steps"))
    for (i in seq_along(project[[field]]))
      .validate_app_person_scoring(project[[field]][[i]]$fit,
        sprintf("the %s history fit at entry %d", field, i))
  for (nm in names(project$kept_fits))
    .validate_app_person_scoring(project$kept_fits[[nm]],
                               sprintf("the kept fit '%s'", nm))
  invisible(project)
}

.save_app_project <- function(project, file) {
  .validate_app_project(project)
  saveRDS(project, file, version = 3, compress = "xz")
  invisible(file)
}

.read_app_project <- function(file) {
  project <- readRDS(file)
  legacy <- is.list(project) &&
    identical(project$format, "rasch-shiny-project") &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 1L
  dropped <- character(0)
  recomputed_weights <- FALSE
  old_contrasts <- is.list(project) && is.list(project$results) &&
    !is.null(project$results[["contrasts"]]) &&
    (!is.list(project$results[["contrasts"]]) ||
     !identical(project$results[["contrasts"]]$algorithm,
                "complete-contrast-cells-2"))
  old_resolution <- is.list(project) && is.list(project$results) &&
    !is.null(project$results[["resolve"]]) &&
    (!is.list(project$results[["resolve"]]) ||
     !identical(project$results[["resolve"]]$algorithm,
                "factor-design-resolution-2") ||
     !is.character(project$results[["resolve"]]$effects) ||
     length(project$results[["resolve"]]$effects) != 1L ||
     anyNA(project$results[["resolve"]]$effects) ||
     !nzchar(trimws(project$results[["resolve"]]$effects)) ||
     !project$results[["resolve"]]$effects %in% c("main", "factorial"))
  # The person-subset comparison withdrew its paired t-test of the two subset
  # means: that test reads the targeting of the split rather than its
  # dimensionality. A saved result without the current stamp carries that
  # reading, so it is dropped rather than restored beside the current verdict.
  old_subtest <- is.list(project) && is.list(project$results) &&
    !is.null(project$results[["subtest"]]) &&
    (!is.list(project$results[["subtest"]]) ||
     !identical(project$results[["subtest"]]$algorithm,
                "person-subset-comparison-1"))
  old_btl_dimensionality <- FALSE
  # The residual decomposition changed from row/count residuals to the
  # pooled expected-score definition. A saved result without the current
  # stamp is not comparable, even when it has no inferential reference.
  old_btl_residual_method <- is.list(project) && is.list(project$results) &&
    inherits(project$results$dimensionality, "rasch_btl_dim") &&
    inherits(project$base_fit, "rasch_btl") &&
    !identical(project$results$dimensionality$residual_method,
               "pooled-expected-score-1")
  legacy_btl_dimension <- is.list(project) && is.list(project$results) &&
    inherits(project$results$dimensionality, "rasch_btl_dim") &&
    is.list(project$results$dimensionality$reference) &&
    inherits(project$base_fit, "rasch_btl") && {
      ref <- project$results$dimensionality$reference
      finite_reference <- any(is.finite(c(
        ref$mean, ref$p95, ref$p, ref$p_adj,
        project$results$dimensionality$bimensions$ref_mean,
        project$results$dimensionality$bimensions$ref_p95)))
      finite_reference && !isTRUE(old_btl_residual_method) &&
        (is.null(ref$inference_available) ||
        ((isTRUE(project$base_fit$clustered) ||
          inherits(project$base_fit, "rasch_btl_efrm")) &&
         !isTRUE(ref$independent_comparisons)))
    }
  old_dim_magnitude <- !legacy && is.list(project) &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 2L && is.list(project$results) &&
    !is.null(project$results$dimension_magnitude) &&
    (!is.list(project$results$dimension_magnitude) ||
     !identical(project$results$dimension_magnitude$algorithm, "complete-panel-1"))
  # Early schema-2 projects contain valid signed bootstrap arrays but predate
  # external (leave-one-out) maxT standardisation. Their adjusted
  # probabilities must not be displayed under the corrected algorithm. The
  # original project seal is checked before the derived result is removed.
  # Rasch loo-maxt-1 results also predate requested-interval preservation.
  old_interval_bootstrap <- is.list(project) && is.list(project$results) &&
    is.list(project$results$bootstrap) &&
    is.list(project$results$bootstrap$bs) &&
    identical(project$results$bootstrap$bs$algorithm, "loo-maxt-1") &&
    !identical(project$results$bootstrap$bs$model_kind, "btl")
  old_maxt <- !legacy && is.list(project) &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 2L &&
    is.list(project$results) &&
    is.list(project$results$bootstrap) &&
    is.list(project$results$bootstrap$bs) &&
    (is.null(project$results$bootstrap$bs$algorithm) || old_interval_bootstrap)
  old_dif_bootstrap <- !legacy && is.list(project) &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 2L &&
    is.list(project$results) &&
    is.list(project$results$dif_bootstrap) &&
    is.list(project$results$dif_bootstrap$db) &&
    is.null(project$results$dif_bootstrap$db$algorithm)
  # Mixed-panel DIF results from before the joint between-person adjustment
  # carry a finite-looking table but are rejected by the current validator.
  # They must be authenticated against the active fit before being omitted;
  # the saved fit/history themselves remain valid and are retained.
  old_mixed_dif <- is.list(project) && is.list(project$results) &&
    inherits(project$results[["dif"]], "rasch_dif") &&
    length(project$results[["dif"]]$within) > 0L &&
    !identical(project$results[["dif"]]$algorithm, "joint-between-1")
  old_dif_followups <- is.list(project) && is.list(project$results) &&
    inherits(project$results[["dif"]], "rasch_dif") &&
    !.dif_followups_current(project$results[["dif"]])
  # Judge-group DIF became tied to the fitted judge role after older app
  # sessions could build maps from a changed sidebar column. A signed result
  # remains structurally valid in that case, so require the source metadata
  # and exact fitted/display judge-ID correspondence before restoring it.
  old_btl_dif_role <- is.list(project) && is.list(project$results) &&
    inherits(project$results[["btl_dif"]], "rasch_btl_dif") &&
    inherits(project$base_fit, "rasch_btl") && {
      meta <- project$results[["btl_dif_meta"]]
      source <- attr(project$base_fit, "rasch_app_source", exact = TRUE)
      role_ok <- is.list(meta) && .app_scalar_text(meta$fitted_judge_col) &&
        .app_scalar_text(meta$judge_col) && is.list(source) &&
        is.data.frame(source$data) && is.list(source$settings) &&
        .app_scalar_text(source$settings$bt_judge) &&
        identical(meta$fitted_judge_col, source$settings$bt_judge) &&
        meta$fitted_judge_col %in% names(source$data) &&
        meta$judge_col %in% names(project$data)
      if (isTRUE(role_ok)) {
        fitted_ids <- as.character(source$data[[meta$fitted_judge_col]])
        display_ids <- as.character(project$data[[meta$judge_col]])
        role_ok <- length(fitted_ids) == length(display_ids) &&
          identical(fitted_ids, display_ids)
      }
      !isTRUE(role_ok)
    }
  old_tailored <- !legacy && is.list(project) &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 2L &&
    is.list(project$results) &&
    !is.null(project$results$guessing) &&
    (!is.list(project$results$guessing) ||
       is.null(project$results$guessing$result_signature) ||
       is.null(project$results$guessing$algorithm) ||
       (identical(project$results$guessing$algorithm, "tailored-four-stage-1") &&
        identical(project$results$guessing$se_method, "bootstrap")) ||
       !"anchor_items_requested" %in% names(project$results$guessing))
  old_frame_invariance <- !legacy && is.list(project) &&
    length(project$schema) == 1L && is.numeric(project$schema) &&
    !is.na(project$schema) && project$schema == 2L &&
    is.list(project$results) &&
    is.list(project$results$frame_invariance) &&
    !all(c("algorithm", "family_n", "boot_reps", "boot_reps_used",
           "boot_reps_nonconverged", "boot_reps_errors",
           "boot_minimum_usable", "bootstrap_stratified") %in%
         names(project$results$frame_invariance))
  # Earlier CJ dimensionality results can retain finite probabilities and
  # reference bands for unsupported comparison designs. Authenticate the
  # original bundle and the result's active-fit binding before omitting only
  # that analysis. Complete, supported legacy references remain usable.
  if (!legacy && isTRUE(old_btl_residual_method) &&
      is.numeric(project$schema) && length(project$schema) == 1L &&
      isTRUE(project$schema == 2L)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    history <- project$btl_steps
    active_fit <- if (length(history)) history[[length(history)]]$fit else
      project$base_fit
    .authenticate_btl_dimensionality(project$results$dimensionality,
                                     active_fit)
    project$results$dimensionality <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped,
                 "Comparative Judgement dimensionality (superseded residual definition; rerun dimensionality)")
  }
  if (!legacy && isTRUE(legacy_btl_dimension) &&
      is.numeric(project$schema) && length(project$schema) == 1L &&
      isTRUE(project$schema == 2L)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    history <- project$btl_steps
    active_fit <- if (length(history)) history[[length(history)]]$fit else
      project$base_fit
    .authenticate_btl_dimensionality(project$results$dimensionality, active_fit)
    old_btl_dimensionality <- .btl_dimensionality_unsupported_reference(
      project$results$dimensionality, active_fit)
    if (old_btl_dimensionality) {
      project$results$dimensionality <- NULL
      project <- .seal_app_project(project)
      dropped <- c(dropped,
                   "Comparative Judgement dimensionality (unsupported reference)")
    }
  }
  if (!legacy && isTRUE(old_mixed_dif) && is.numeric(project$schema) &&
      length(project$schema) == 1L && isTRUE(project$schema == 2L)) {
    # Authenticate the enclosing bundle and the saved primary tables before
    # omitting the obsolete result. The saved fit/history remain intact.
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    history <- project$rasch_steps
    active_fit <- if (length(history)) history[[length(history)]]$fit else
      project$base_fit
    old_dif <- project$results[["dif"]]
    .validate_primary_dif_tables(old_dif, "item", c("Residuals", "ci"), ":ci")
    if (is.null(old_dif$fit_signature) ||
        !.fit_boot_signature_matches(old_dif$fit_signature, active_fit))
      stop("the saved mixed-panel DIF was computed from a different fitted model",
           call. = FALSE)
    project$results$dif <- NULL
    project$results$dif_bootstrap <- NULL
    project$results$resolve <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped,
                 "DIF (superseded mixed-panel adjustment)",
                 "DIF bootstrap (dependent on superseded DIF)",
                 "automatic DIF resolution (dependent on superseded DIF)")
  }
  if (isTRUE(old_dif_followups) && !isTRUE(old_mixed_dif) &&
      is.numeric(project$schema) && length(project$schema) == 1L &&
      isTRUE(project$schema == 2L)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    history <- project$rasch_steps
    active_fit <- if (length(history)) history[[length(history)]]$fit else
      project$base_fit
    saved_dif <- project$results[["dif"]]
    if (!is.null(saved_dif$fit_signature)) {
      problem <- tryCatch({
        .validate_primary_dif_tables(saved_dif, "item", c("Residuals", "ci"), ":ci")
        if (!.fit_boot_signature_matches(saved_dif$fit_signature, active_fit))
          stop("`dif` was computed from a different fitted model")
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem))
        stop(paste("the saved DIF cannot be authenticated before migration:",
                   problem), call. = FALSE)
    }
    has_bootstrap <- !is.null(project$results$dif_bootstrap)
    has_resolution <- !is.null(project$results$resolve)
    project$results$dif <- NULL
    project$results$dif_bootstrap <- NULL
    project$results$resolve <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped,
                 "DIF (superseded normalized-factor follow-ups)",
                 if (has_bootstrap)
                   "DIF bootstrap (dependent on superseded DIF)",
                 if (has_resolution)
                   "automatic DIF resolution (dependent on superseded DIF)")
  }
  if (isTRUE(old_btl_dif_role) && is.numeric(project$schema) &&
      length(project$schema) == 1L && isTRUE(project$schema == 2L)) {
    # Authenticate the enclosing project before omitting the old derived
    # result. A fit-mismatched signed result is not made acceptable merely by
    # classifying its judge metadata as obsolete; an unsigned result is simply
    # unverifiable and is omitted with the other legacy derived results.
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    history <- project$btl_steps
    active_fit <- if (length(history)) history[[length(history)]]$fit else
      project$base_fit
    saved_btl_dif <- project$results[["btl_dif"]]
    if (!is.null(saved_btl_dif$fit_signature)) {
      problem <- tryCatch({
        .validate_btl_dif_result(saved_btl_dif, active_fit)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem))
        stop(paste("the saved Comparative Judgement DIF cannot be",
                   "authenticated before migration:", problem), call. = FALSE)
    }
    has_bootstrap <- !is.null(project$results$dif_bootstrap)
    project$results$btl_dif <- NULL
    project$results$btl_dif_meta <- NULL
    project$results$dif_bootstrap <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped,
                 "Comparative Judgement DIF (missing fitted judge-role provenance)",
                 if (has_bootstrap)
                   "DIF bootstrap (dependent on superseded Comparative Judgement DIF)")
  }
  # Earlier planned contrasts could renormalise away unresolved weighted
  # cells. Keep the source and fitted models, but do not restore estimates
  # for a different comparison. Check the original bundle before resealing.
  if (!legacy && (isTRUE(old_contrasts) || isTRUE(old_resolution)) &&
      is.numeric(project$schema) &&
      length(project$schema) == 1L && isTRUE(project$schema == 2L)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    if (isTRUE(old_contrasts)) {
      project$results$contrasts <- NULL
      dropped <- c(dropped, "planned DIF contrasts (superseded cell support)")
    }
    if (isTRUE(old_resolution)) {
      project$results$resolve <- NULL
      dropped <- c(dropped,
                   "automatic DIF resolution (superseded factor model or reporting)")
    }
    project <- .seal_app_project(project)
  }
  if (isTRUE(old_maxt)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    project$results$bootstrap <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped, if (old_interval_bootstrap)
      "fit bootstrap (superseded class-interval allocation)" else
      "fit bootstrap (superseded maxT adjustment)")
  }
  # The earlier conditional DIF bootstrap used raw F values for its marginal
  # empirical probability. Degrees of freedom can change across sparse refits,
  # so those values are not on a common reference scale. Authenticate the
  # complete project before omitting only that derived result.
  if (isTRUE(old_dif_bootstrap) && !isTRUE(old_mixed_dif) &&
      !isTRUE(old_dif_followups) && !isTRUE(old_btl_dif_role)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    project$results$dif_bootstrap <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped,
                 "DIF bootstrap (superseded raw-F marginal reference)")
  }
  # Omit unverifiable or superseded tailored results, retaining their source
  # data and fits after checking the complete project seal.
  if (isTRUE(old_tailored)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    tailored_reason <- if (
      is.list(project$results$guessing) &&
      identical(project$results$guessing$algorithm, "tailored-four-stage-1") &&
      identical(project$results$guessing$se_method, "bootstrap"))
      "superseded bootstrap" else "unverifiable fitted model"
    project$results$guessing <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped, paste0("tailored analysis (", tailored_reason, ")"))
  }
  # Earlier schema-2 projects can carry frame-invariance results produced
  # before complete multiplicity/bootstrap accounting and the strict
  # all-frame comparison rule. Authenticate the bundle before omitting that
  # derived result.
  if (isTRUE(old_frame_invariance)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    project$results$frame_invariance <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped, "frame-invariance analysis (superseded inference)")
  }
  # Older magnitudes can compare PSI from different response samples. Verify
  # the original bundle before removing that result and retaining its fits.
  if (isTRUE(old_dim_magnitude)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    project$results$dimension_magnitude <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped, "dimensionality magnitude (unmatched reliability samples)")
  }
  # Keep the source data, fits and history of a project whose person-subset
  # test predates the current comparison; omit only that result. Verify the
  # original bundle before resealing.
  if (!legacy && isTRUE(old_subtest) && is.numeric(project$schema) &&
      length(project$schema) == 1L && isTRUE(project$schema == 2L)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    project$results$subtest <- NULL
    project <- .seal_app_project(project)
    dropped <- c(dropped,
                 "person-subset dimensionality test (superseded subset-mean comparison)")
  }
  # Earlier weighted solvers have no algorithm stamp, use pattern-wle-1, or
  # use the intermediate pattern-unit-wle-2 implementation.
  # They can differ in their last bits, or fail for tiny observed weights
  # or large changes of measurement unit.
  # Authenticate the saved bundle and result before replacing that derived
  # table. Current results still require exact reproduction.
  old_weights <- is.list(project) && is.list(project$results) &&
    is.list(project$results$person_weights) &&
    (is.null(attr(project$results$person_weights$table, "algorithm", exact = TRUE)) ||
     identical(attr(project$results$person_weights$table, "algorithm", exact = TRUE),
               "pattern-wle-1") ||
     identical(attr(project$results$person_weights$table, "algorithm", exact = TRUE),
               "pattern-unit-wle-2"))
  if (!legacy && isTRUE(old_weights) && is.numeric(project$schema) &&
      length(project$schema) == 1L && isTRUE(project$schema == 2L)) {
    unsigned_project <- project
    attr(unsigned_project, "rasch_project_legacy") <- NULL
    attr(unsigned_project, "rasch_project_legacy_dropped") <- NULL
    unsigned_project$binding <- NULL
    if (!.app_scalar_text(project$binding) ||
        !.fit_boot_hash_matches(project$binding, unsigned_project))
      stop(paste("the analysis file's source data, fitted models or results have",
                 "changed since they were saved"), call. = FALSE)
    history <- if (identical(.app_fit_family(project$base_fit), "btl"))
      project$btl_steps else project$rasch_steps
    active_fit <- if (length(history)) history[[length(history)]]$fit else
      project$base_fit
    weighted <- project$results$person_weights
    .authenticate_weighted_person_result(weighted, active_fit)
    weighted$table <- weighted_person_estimates(
      active_fit, weighted$weights, by = weighted$by, sets = weighted$sets)
    weighted$result_signature <- NULL
    weighted$result_signature <- .fit_boot_md5(weighted)
    project$results$person_weights <- weighted
    project <- .seal_app_project(project)
    recomputed_weights <- TRUE
  }
  if (legacy) {
    # Schema 1 did not record an integrity binding. It can be checked
    # structurally and upgraded, but its original data-to-fit relationship
    # cannot be authenticated retrospectively. Results written before the
    # current result fingerprints cannot be validated against the retained
    # fit, so omit those results while preserving the data, fit and history.
    results <- project$results
    has_signature <- function(x)
      is.list(x) && is.character(x$result_signature) &&
        length(x$result_signature) == 1L && !is.na(x$result_signature)
    history <- if (identical(.app_fit_family(project$base_fit), "btl"))
      project$btl_steps else
      project$rasch_steps
    active_fit <- if (length(history)) history[[length(history)]]$fit else
      project$base_fit
    # Schema 1 has no enclosing binding, but a result written by a current
    # run can still carry an authenticated fit signature. Do not silently
    # hide a signed result whose fit or primary tables were tampered with;
    # only genuinely unverified legacy results are omitted below.
    if (!isTRUE(old_mixed_dif) && isTRUE(old_dif_followups) &&
        is.list(results$dif) && !is.null(results$dif$fit_signature)) {
      problem <- tryCatch({
        .validate_primary_dif_tables(results$dif, "item",
                                     c("Residuals", "ci"), ":ci")
        if (!.fit_boot_signature_matches(results$dif$fit_signature,
                                         active_fit))
          stop("`dif` was computed from a different fitted model")
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem))
        stop(paste("the saved DIF cannot be authenticated before migration:",
                   problem), call. = FALSE)
    }
    if (isTRUE(old_btl_dif_role) && is.list(results$btl_dif) &&
        !is.null(results$btl_dif$fit_signature)) {
      problem <- tryCatch({
        .validate_btl_dif_result(results$btl_dif, active_fit)
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(problem))
        stop(paste("the saved Comparative Judgement DIF cannot be",
                   "authenticated before migration:", problem), call. = FALSE)
    }
    if (is.list(results)) {
      if (isTRUE(old_contrasts)) {
        results$contrasts <- NULL
        dropped <- c(dropped, "planned DIF contrasts")
      }
      if (isTRUE(old_resolution)) {
        results$resolve <- NULL
        dropped <- c(dropped, "automatic DIF resolution")
      }
      if (isTRUE(old_btl_residual_method) &&
          !is.null(results$dimensionality)) {
        signature <- attr(results$dimensionality, "result_signature", exact = TRUE)
        if (.app_scalar_text(signature)) {
          problem <- tryCatch({
            .authenticate_btl_dimensionality(results$dimensionality,
                                             active_fit)
            NULL
          }, error = function(e) conditionMessage(e))
          if (!is.null(problem))
            stop(paste("the saved Comparative Judgement dimensionality cannot be",
                       "authenticated before migration:", problem),
                 call. = FALSE)
        }
        results$dimensionality <- NULL
        dropped <- c(dropped,
                     "Comparative Judgement dimensionality (superseded residual definition; rerun dimensionality)")
      }
      if (isTRUE(legacy_btl_dimension)) {
        signature <- attr(results$dimensionality, "result_signature", exact = TRUE)
        omit <- !.app_scalar_text(signature)
        if (!omit) {
          history <- project$btl_steps
          active_fit <- if (length(history)) history[[length(history)]]$fit else
            project$base_fit
          .authenticate_btl_dimensionality(results$dimensionality, active_fit)
          omit <- .btl_dimensionality_unsupported_reference(
            results$dimensionality, active_fit)
        }
        if (omit) {
          results$dimensionality <- NULL
          dropped <- c(dropped, "Comparative Judgement dimensionality")
        }
      }
      if (isTRUE(old_weights)) {
        results$person_weights <- NULL
        dropped <- c(dropped, "externally weighted person estimates")
      }
      if (isTRUE(old_mixed_dif) && !is.null(results$dif)) {
        results$dif <- NULL
        results$dif_bootstrap <- NULL
        results$resolve <- NULL
        dropped <- c(dropped,
                     "DIF (superseded mixed-panel adjustment)",
                     "DIF bootstrap (dependent on superseded DIF)",
                     "automatic DIF resolution (dependent on superseded DIF)")
      }
      if (isTRUE(old_dif_followups) && !isTRUE(old_mixed_dif) &&
          !is.null(results$dif)) {
        has_bootstrap <- !is.null(results$dif_bootstrap)
        has_resolution <- !is.null(results$resolve)
        results$dif <- NULL
        results$dif_bootstrap <- NULL
        results$resolve <- NULL
        dropped <- c(dropped,
                     "DIF (superseded normalized-factor follow-ups)",
                     if (has_bootstrap)
                       "DIF bootstrap (dependent on superseded DIF)",
                     if (has_resolution)
                       "automatic DIF resolution (dependent on superseded DIF)")
      }
      if (!is.null(results$dimension_magnitude)) {
        results$dimension_magnitude <- NULL
        dropped <- c(dropped, "dimensionality magnitude")
      }
      if (isTRUE(old_subtest)) {
        results$subtest <- NULL
        dropped <- c(dropped, "person-subset dimensionality test")
      }
      if (!is.null(results$bootstrap) &&
          (old_interval_bootstrap || !is.list(results$bootstrap) ||
           !has_signature(results$bootstrap$bs))) {
        results$bootstrap <- NULL
        dropped <- c(dropped, "fit bootstrap")
      }
      primary_dropped <- FALSE
      if (isTRUE(old_btl_dif_role) && !is.null(results$btl_dif)) {
        has_bootstrap <- !is.null(results$dif_bootstrap)
        results$btl_dif <- NULL
        results$btl_dif_meta <- NULL
        results$dif_bootstrap <- NULL
        dropped <- c(dropped,
                     "Comparative Judgement DIF (missing fitted judge-role provenance)",
                     if (has_bootstrap)
                       "DIF bootstrap (dependent on superseded Comparative Judgement DIF)")
        primary_dropped <- TRUE
      }
      for (nm in c("dif", "btl_dif")) {
        if (!is.null(results[[nm]]) && !has_signature(results[[nm]])) {
          results[[nm]] <- NULL
          dropped <- c(dropped, if (nm == "btl_dif")
            "Comparative Judgement DIF" else "DIF")
          primary_dropped <- TRUE
        }
      }
      if (isTRUE(primary_dropped)) results$btl_dif_meta <- NULL
      if (!is.null(results$dif_bootstrap) &&
          (isTRUE(primary_dropped) || !is.list(results$dif_bootstrap) ||
           !has_signature(results$dif_bootstrap$db))) {
        results$dif_bootstrap <- NULL
        dropped <- c(dropped, "DIF bootstrap")
      }
      if (!is.null(results$guessing) && !has_signature(results$guessing)) {
        results$guessing <- NULL
        dropped <- c(dropped, "tailored analysis")
      }
      project$results <- results
    }
    project <- .seal_app_project(project)
  }
  .validate_app_project(project)
  if (recomputed_weights)
    warning(paste("saved weighted person estimates were recomputed using",
                  "the current scoring algorithm; save the analysis again",
                  "to retain the updated table"), call. = FALSE)
  if (legacy) {
    attr(project, "rasch_project_legacy") <- TRUE
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("this schema-1 analysis predates data-to-fit integrity",
                  "checks; it passed structural validation and has been",
                  "upgraded in memory; save it again to retain schema 2",
                  if (length(dropped)) paste0("; unverifiable saved results ",
                    "were omitted: ", paste(unique(dropped), collapse = ", "))
                  else ""),
            call. = FALSE)
  }
  if (isTRUE(old_maxt)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved fit bootstrap used the earlier",
                  if (old_interval_bootstrap) "class-interval allocation" else
                    "maxT standardisation",
                  "and was omitted; recompute it before",
                  "reporting adjusted bootstrap probabilities"),
            call. = FALSE)
  }
  if (isTRUE(old_dif_bootstrap) && !isTRUE(old_mixed_dif) &&
      !isTRUE(old_dif_followups) && !isTRUE(old_btl_dif_role)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved DIF bootstrap used the earlier raw-F marginal",
                  "reference and was omitted; recompute it before reporting",
                  "bootstrap DIF probabilities"), call. = FALSE)
  }
  if (isTRUE(old_tailored)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved tailored analysis predates current result",
                  "provenance and was omitted; recompute it before",
                  "reporting tailored item shifts"), call. = FALSE)
  }
  if (isTRUE(old_frame_invariance)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved frame-invariance analysis used earlier",
                  "comparison or bootstrap-accounting rules and was omitted;",
                  "recompute it before reporting frame-invariance inference"),
            call. = FALSE)
  }
  if (isTRUE(old_dim_magnitude)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved dimensionality magnitude used earlier reliability",
                  "samples and was omitted; recompute it on matched response rows"),
            call. = FALSE)
  }
  if (!legacy && isTRUE(old_subtest)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved person-subset dimensionality test reported the",
                  "superseded paired t-test of the subset means and was",
                  "omitted; recompute it with dimensionality_test() before",
                  "reporting its verdict"), call. = FALSE)
  }
  if (!legacy && isTRUE(old_btl_residual_method)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved Comparative Judgement dimensionality used an",
                  "earlier residual definition and was omitted; rerun",
                  "dimensionality before reporting its decomposition or inference"),
            call. = FALSE)
  }
  if (isTRUE(old_btl_dimensionality)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved Comparative Judgement dimensionality reference",
                  "used an unsupported comparison design and was omitted;",
                  "recompute it to retain the observed decomposition",
                  "without unsupported inference"), call. = FALSE)
  }
  if (!legacy && isTRUE(old_mixed_dif)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved mixed-panel DIF analysis predates the joint",
                  "between-person adjustment and was omitted, along with",
                  "its dependent bootstrap and resolution; recompute DIF",
                  "before reporting inference"), call. = FALSE)
  }
  if (!legacy && isTRUE(old_dif_followups) && !isTRUE(old_mixed_dif)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved DIF follow-ups may use earlier factor values;",
                  "recompute the DIF analysis before reporting inference"),
            call. = FALSE)
  }
  if (!legacy && isTRUE(old_btl_dif_role)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved Comparative Judgement DIF lacked authenticated",
                  "fitted judge-role provenance and was omitted; refit DIF",
                  "before reporting inference"),
            call. = FALSE)
  }
  if (!legacy && isTRUE(old_contrasts)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved planned DIF contrasts predate the current",
                  "complete-cell support rules and were omitted; recompute",
                  "them before reporting estimates or probabilities"),
            call. = FALSE)
  }
  if (!legacy && isTRUE(old_resolution)) {
    attr(project, "rasch_project_legacy_dropped") <- unique(dropped)
    warning(paste("the saved automatic DIF resolution predates the current",
                  "factor model and reporting rules and was omitted; fitted",
                  "models and analysis history are unchanged; rerun resolution",
                  "from the pre-resolution fit with the intended factors and effects"),
            call. = FALSE)
  }
  project
}
