# rasch :: export
# ===========================================================================
# save_outputs() writes the complete analysis to disk: every table as CSV,
# every plot as PNG (and optionally PDF), and a plain-text summary. The
# The Shiny app zips the resulting fitted-model folder for its model-results
# archive.
# ===========================================================================

# A file stem drops the characters a file system cannot carry, which can map
# two different item names onto one file: the later plot would overwrite the
# earlier and the export would still look complete. Disambiguate after
# sanitising so every plot keeps its own file.
# Device dimensions are opened before anything is drawn: an invalid size
# fails every plot in a batch, which is a caller error worth stating once
# rather than reporting as many drawing failures.
.check_pos_num <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      !is.null(dim(x)) || !is.null(oldClass(x)) || !is.finite(x) || x <= 0)
    stop("`", name, "` must be one positive finite value", call. = FALSE)
  invisible(x)
}

.check_out_path <- function(x, name) {
  if (length(x) != 1L || !is.character(x) || !is.null(dim(x)) ||
      !is.null(oldClass(x)) || is.na(x) || !nzchar(trimws(x)))
    stop("`", name, "` must be one non-missing path", call. = FALSE)
  invisible(x)
}

.check_export_dir <- function(x) {
  .check_out_path(x, "dir")
  if (dir.exists(x)) {
    entries <- list.files(x, all.files = TRUE, no.. = TRUE)
    if (length(entries))
      stop("`dir` must be a new or empty directory; refusing to overwrite ",
           "existing export files", call. = FALSE)
  } else if (file.exists(x)) {
    stop("`dir` must be a new path or an existing empty directory; the ",
         "supplied path is an existing file", call. = FALSE)
  }
  invisible(x)
}

.check_export_fit <- function(fit, btl = TRUE) {
  ok <- inherits(fit, "rasch") || (isTRUE(btl) && inherits(fit, "rasch_btl"))
  if (!ok)
    stop("`fit` must be a fitted ",
         if (isTRUE(btl)) "Rasch or paired-comparison model"
         else "person-by-item Rasch model", call. = FALSE)
  if (inherits(fit, "rasch_btl")) {
    if (!isTRUE(fit$converged))
      stop("the paired-comparison calibration did not converge; a complete export is unavailable",
           call. = FALSE)
  } else {
    if (!isTRUE(fit$est$converged))
      stop("the fitted calibration did not converge; a complete export is unavailable",
           call. = FALSE)
    if (!.efrm_link_converged(fit))
      stop("the fitted set-unit link did not converge; a complete export is unavailable",
           call. = FALSE)
  }
  invisible(fit)
}

.check_device_size <- function(width, height, dpi) {
  .check_pos_num(width, "width")
  .check_pos_num(height, "height")
  .check_pos_num(dpi, "dpi")
  invisible(TRUE)
}

.rr_safe_stem <- function(x) {
  s <- gsub("[^A-Za-z0-9_.-]", "_", as.character(x))
  make.unique(s, sep = "_")
}

.rr_device <- function(path, fmt, width, height, dpi) {
  if (fmt == "png") png(path, width = width, height = height, units = "in", res = dpi)
  else pdf(path, width = width, height = height)
}

.rr_save_plot <- function(expr, stem, dir, formats, width, height, dpi) {
  .check_device_size(width, height, dpi)
  files <- character(0)
  for (fmt in formats) {
    path <- file.path(dir, paste0(stem, ".", fmt))
    # close only a device this call opened: if the device itself failed to
    # open, dev.off() would close the caller's device instead
    before <- grDevices::dev.cur()
    ok <- tryCatch({
      .rr_device(path, fmt, width, height, dpi)
      force(expr())
      dev.off()
      TRUE
    }, error = function(e) {
      if (!identical(grDevices::dev.cur(), before))
        try(dev.off(), silent = TRUE)
      # an omitted display must not pass silently for an export the user
      # will treat as complete
      warning("plot '", stem, "' could not be drawn: ",
              conditionMessage(e), call. = FALSE)
      FALSE
    })
    if (ok) files <- c(files, path)
  }
  files
}

# One plot per element, to a multi-page PDF or a ZIP of PNGs by extension.
.rr_plot_batch <- function(thunks, names, stem, file, width, height, dpi) {
  .check_out_path(file, "file")
  ext <- tolower(tools::file_ext(file))
  if (!grepl("^(/|[A-Za-z]:)", file)) file <- file.path(getwd(), file)
  .check_device_size(width, height, dpi)
  failed <- character(0)
  if (ext == "pdf") {
    pdf(file, width = width, height = height, onefile = TRUE)
    drawn <- 0L
    tryCatch({
      for (j in seq_along(thunks))
        tryCatch({ thunks[[j]](); drawn <- drawn + 1L },
                 error = function(e) failed <<- c(failed, names[j]))
    }, finally = dev.off())
    # a file the caller would treat as the export must not be returned when
    # nothing was drawn into it
    if (!drawn) {
      unlink(file)
      stop("no plots could be drawn, so no file was created: ",
           paste(failed, collapse = ", "), call. = FALSE)
    }
  } else if (ext == "zip") {
    dir <- tempfile("rasch_plots_"); dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    paths <- character(0)
    safe <- .rr_safe_stem(names)
    for (j in seq_along(thunks)) {
      p <- file.path(dir, paste0(stem, "_", safe[j], ".png"))
      before <- grDevices::dev.cur()
      ok <- tryCatch({
        png(p, width = width, height = height, units = "in", res = dpi)
        thunks[[j]](); dev.off(); TRUE
      }, error = function(e) {
        if (!identical(grDevices::dev.cur(), before))
          try(dev.off(), silent = TRUE)
        failed <<- c(failed, names[j]); FALSE })
      if (ok) paths <- c(paths, p)
    }
    if (!length(paths))
      stop("no plots could be drawn, so no archive was created: ",
           paste(failed, collapse = ", "), call. = FALSE)
    # zip() ADDS to an existing archive: without this the returned path
    # would carry plots from an earlier call and read as this call's set
    unlink(file)
    wd <- getwd(); setwd(dir); on.exit(setwd(wd), add = TRUE)
    status <- utils::zip(file, files = basename(paths), flags = "-q9X")
    setwd(wd)
    if (!isTRUE(as.integer(status) == 0L) || !file.exists(file))
      stop("the plot archive could not be created at ", file, call. = FALSE)
  } else stop("`file` must end in .pdf or .zip")
  if (length(failed))
    warning("plot(s) could not be drawn and were omitted: ",
            paste(failed, collapse = ", "), call. = FALSE)
  invisible(file)
}

#' Save a plot for every item
#'
#' Writes one plot per item -- the item characteristic curve, category
#' probability curves, threshold probability curves, or category
#' frequencies -- to a single multi-page PDF or a ZIP archive of PNGs,
#' chosen by the extension of \code{file}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param what Which plot: \code{"icc"}, \code{"ccc"}, \code{"tpc"}, or
#'   \code{"cfreq"}.
#' @param file Output path ending in \code{.pdf} (one page per item) or
#'   \code{.zip} (one PNG per item).
#' @param items Item names or indices; all items by default.
#' @param n_groups Class intervals for observed overlays; \code{NULL} retains
#'   the fit's allocation for each item.
#' @param grid Logit grid for the curves.
#' @param observed Overlay observed proportions on the category and
#'   threshold probability curves.
#' @param width,height,dpi Device size in inches and PNG resolution.
#' @return Invisibly, the output path.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' f <- rasch(X)
#' save_item_plots(f, "icc", file.path(tempdir(), "icc_all.pdf"))
#' @export
save_item_plots <- function(fit, what = c("icc", "ccc", "tpc", "cfreq"),
                            file, items = NULL, n_groups = NULL,
                            grid = NULL, observed = TRUE,
                            width = 8, height = 5.5, dpi = 300) {
  .check_export_fit(fit, btl = FALSE)
  what <- match.arg(what)
  its <- if (is.null(items)) fit$items$item else {
    if (!is.atomic(items) || !is.null(dim(items)) || !length(items))
      stop("`items` must be a non-empty ordinary vector of names or indices",
           call. = FALSE)
    if (is.numeric(items)) {
      if (any(!is.finite(items)) || any(items != floor(items)))
        stop("numeric `items` must contain whole indices", call. = FALSE)
      canonical <- as.character(items)
    } else {
      supplied <- as.character(items)
      if (anyNA(supplied) || any(!nzchar(trimws(supplied))))
        stop("`items` must contain non-missing, non-empty item names",
             call. = FALSE)
      ii <- match(supplied, fit$items$item)
      fallback <- is.na(ii)
      if (any(fallback))
        ii[fallback] <- match(.role_text_values(supplied[fallback]),
                              fit$items$item)
      canonical <- supplied
      canonical[!is.na(ii)] <- fit$items$item[ii[!is.na(ii)]]
    }
    if (anyDuplicated(canonical))
      stop("`items` must not name the same item more than once", call. = FALSE)
    if (is.numeric(items)) items else canonical
  }
  draw <- function(it) switch(what,
    icc   = plot_icc(fit, it, n_groups = n_groups, grid = grid,
                     observed = observed),
    ccc   = plot_ccc(fit, it, grid = grid, observed = observed,
                     n_groups = n_groups),
    tpc   = plot_threshold_prob(fit, it, grid = grid, observed = observed,
                                n_groups = n_groups),
    cfreq = plot_catfreq(fit, it))
  thunks <- lapply(its, function(it) function() draw(it))
  .rr_plot_batch(thunks, as.character(its), what, file, width, height, dpi)
}

#' Save a kidmap for every person
#'
#' Writes one kidmap (\code{\link{plot_kidmap}}) per person to a single
#' multi-page PDF or a ZIP archive of PNGs, chosen by the extension of
#' \code{file}. Persons without a location estimate are skipped.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param file Output path ending in \code{.pdf} (one page per person) or
#'   \code{.zip} (one PNG per person).
#' @param persons Row numbers or IDs; all estimated persons by default.
#' @param level Confidence level of the band marking unexpected responses.
#' @param width,height,dpi Device size in inches and PNG resolution.
#' @return Invisibly, the output path.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(60 * 6, 1, plogis(outer(rnorm(60), d, "-"))), 60, 6)
#' colnames(X) <- paste0("I", 1:6)
#' f <- rasch(X)
#' save_person_plots(f, file.path(tempdir(), "kidmaps.pdf"), persons = 1:5)
#' @export
save_person_plots <- function(fit, file, persons = NULL, level = 0.95,
                              width = 8, height = 6, dpi = 300) {
  .check_export_fit(fit, btl = FALSE)
  ps <- if (is.null(persons)) which(!is.na(fit$person$theta)) else persons
  if (!is.atomic(ps) || !is.null(dim(ps)) || !length(ps) || anyNA(ps))
    stop("`persons` must be a non-empty ordinary vector of row numbers or IDs",
         call. = FALSE)
  if (is.numeric(ps)) {
    if (any(!is.finite(ps)) || any(ps != floor(ps)) || any(ps < 1L) ||
        any(ps > nrow(fit$person)))
      stop("numeric `persons` must be whole row numbers in 1..",
           nrow(fit$person), call. = FALSE)
    ps <- as.integer(ps)
    if (anyDuplicated(ps))
      stop("`persons` must not name the same row more than once", call. = FALSE)
    ids <- fit$person$id[ps]
  } else {
    ps <- as.character(ps)
    if (any(!nzchar(trimws(ps))))
      stop("person IDs must be non-empty", call. = FALSE)
    hits <- lapply(ps, function(p)
      which(as.character(fit$person$id) == p))
    if (any(!lengths(hits)))
      stop("person ID(s) not found: ",
           paste(ps[!lengths(hits)], collapse = ", "), call. = FALSE)
    repeated <- lengths(hits) > 1L
    if (any(repeated))
      stop("person ID(s) occur on several response rows: ",
           paste(ps[repeated], collapse = ", "),
           "; use row numbers instead", call. = FALSE)
    if (anyDuplicated(ps))
      stop("`persons` must not name the same ID more than once", call. = FALSE)
    ids <- ps
  }
  thunks <- lapply(ps, function(p)
    function() plot_kidmap(fit, p, level = level))
  .rr_plot_batch(thunks, as.character(ids), "kidmap", file,
                 width, height, dpi)
}

.save_btl_outputs <- function(fit, dir, formats, width, height, dpi,
                              object_plots, bootstrap = NULL, dif = NULL,
                              dif_bootstrap = NULL, dimensionality = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tdir <- file.path(dir, "tables"); pdir <- file.path(dir, "plots")
  odir <- file.path(pdir, "objects")
  dir.create(tdir, showWarnings = FALSE); dir.create(pdir, showWarnings = FALSE)
  if (object_plots) dir.create(odir, showWarnings = FALSE)
  files <- character(0)
  wtab <- function(d, name) {
    if (is.null(d)) return()
    path <- file.path(tdir, paste0(name, ".csv"))
    .write_csv_plain(d, path)
    files <<- c(files, path)
  }
  wtab(fit_summary_table(fit), "fit_summary")
  wtab(fit$objects, "object_estimates")
  wtab(fit$pairs, "pair_fit")
  wtab(fit$judges, "judge_fit")
  wtab(fit$comparisons, "comparisons")
  if (!is.null(bootstrap)) {
    wtab(bootstrap$objects, "bootstrap_object_fit")
    wtab(bootstrap$pairs, "bootstrap_pair_fit")
    wtab(bootstrap$judges, "bootstrap_judge_fit")
    wtab(data.frame(quantity = names(unlist(bootstrap$total)),
                    value = unlist(bootstrap$total)), "bootstrap_total")
  }
  if (!is.null(dif)) {
    wtab(dif$summary, "dif_anova")
    wtab(dif$terms, "dif_anova_terms")
    wtab(dif$levels, "dif_resolved_locations")
    wtab(dif$sizes, "dif_magnitudes")
  }
  if (!is.null(dif_bootstrap)) {
    wtab(dif_bootstrap$summary, "dif_bootstrap")
    wtab(dif_bootstrap$terms, "dif_bootstrap_terms")
    wtab(data.frame(
      quantity = c("requested", "used", "failed", "nonconverged",
                   "errors", "family_tests"),
      value = c(dif_bootstrap$B, dif_bootstrap$B_used,
                dif_bootstrap$B_failed, dif_bootstrap$B_nonconverged,
                dif_bootstrap$B_errors, dif_bootstrap$family_n)),
      "dif_bootstrap_accounting")
  }
  if (inherits(fit, "rasch_btl_explanatory")) {
    wtab(explanatory_test(fit), "explanatory_model_comparison")
    wtab(fit$object_coefficients, "explanatory_predictor_effects")
    wtab(explanatory_diagnostics(fit), "explanatory_fixed_departure_diagnostics")
    if (nrow(fit$explanatory$relaxations))
      wtab(fit$explanatory$relaxations, "explanatory_fixed_departures")
  }
  info <- tryCatch(btl_information(fit), error = function(e) NULL)
  if (!is.null(info)) {
    wtab(info$objects, "object_information")
    wtab(info$pairs, "pair_information")
  }
  tr <- tryCatch(btl_transitivity(fit), error = function(e) NULL)
  if (!is.null(tr)) {
    wtab(tr$objects, "transitivity_objects")
    wtab(tr$judges, "transitivity_judges")
  }
  bd <- if (!is.null(dimensionality)) dimensionality else
    tryCatch(btl_dimensionality(fit, seed = 1L), error = function(e) e)
  if (inherits(bd, "error")) {
    wtab(data.frame(note = paste("Dimensionality was not available:",
                                 conditionMessage(bd))),
         "residual_dimensionality")
  } else {
    wtab(bd$bimensions, "residual_bimensions")
    wtab(bd$coords, "residual_bimension_coordinates")
    ref <- bd$reference
    wtab(data.frame(
      mean = ref$mean, upper_5_percent = ref$p95, p = ref$p,
      p_adj = ref$p_adj, requested = ref$reps, used = ref$n_used,
      alpha = ref$alpha, seed = ref$seed %||% NA_integer_,
      method = ref$method,
      inference_available = .btl_dimensionality_has_inference(bd),
      independent_comparisons = ref$independent_comparisons %||% NA,
      note = paste(bd$notes, collapse = " ")),
      "residual_dimensionality_reference")
  }
  if (inherits(fit, "rasch_btl_efrm")) {
    wtab(fit$phi_table, "panel_units_phi")
    wtab(fit$alpha_table, "set_units_alpha")
    wtab(fit$kappa_table, "set_origins_kappa")
    wtab(fit$unit_omnibus, "unit_omnibus_tests")
    wtab(fit$frames, "frame_fit")
    wtab(data.frame(
      model = c("Equal units", "Frame-dependent units"),
      loglik = c(fit$equal_unit$loglik_single,
                 fit$equal_unit$loglik_frames),
      parameters = c(fit$equal_unit$parameters_single,
                     fit$equal_unit$parameters_frames),
      two_delta_loglik = c(NA_real_, fit$equal_unit$two_delta_ll)),
      "frame_model_comparison")
  }
  spath <- file.path(dir, "summary.txt")
  # a blank row in the exported DIF table is a test the design could not
  # estimate, not an absence of DIF: the summary carries the reason
  writeLines(c(utils::capture.output(print(fit)), "", fit$notes,
               if (length(dif$notes))
                 paste("DIF notes:", paste(dif$notes, collapse = "; "))),
             spath)
  files <- c(files, spath)
  sp <- function(f, stem) files <<- c(files,
    .rr_save_plot(f, stem, pdir, formats, width, height, dpi))
  sp(function() plot_btl(fit), "object_locations")
  sp(function() plot_btl_targeting(fit), "design_information")
  if (!is.null(tr))
    sp(function() plot_btl_transitivity(tr), "transitivity")
  if (!inherits(bd, "error")) {
    sp(function() plot_btl_scree(bd), "residual_bimensions")
    sp(function() plot_btl_dim_map(bd), "residual_bimension_map")
  }
  if (inherits(fit, "rasch_btl_efrm"))
    sp(function() plot_btl_units(fit), "frame_units")
  if (object_plots) {
    obs <- fit$objects$object
    osafe <- .rr_safe_stem(obs)
    for (j in seq_along(obs)) local({
      object <- obs[j]; s_ <- osafe[j]
      files <<- c(files, .rr_save_plot(
        function() plot_btl_icc(fit, object),
        paste0(s_, "_icc"), odir, formats, width, height, dpi))
    })
  }
  invisible(files)
}

.dimensionality_or_note <- function(fit) {
  tryCatch(dimensionality_test(fit), error = function(e)
    list(note = conditionMessage(e), multidimensional = NA))
}

# A deliberate design refusal belongs in the report. Unexpected calculation
# errors still stop the export rather than being disguised as missing results.
.distractors_or_note <- function(fit) {
  if (is.null(fit$mc)) return(NULL)
  tryCatch(list(table = distractor_analysis(fit), note = NULL),
    rasch_refusal = function(e)
      list(table = NULL, note = conditionMessage(e)))
}

#' Save the outputs of a Rasch analysis
#'
#' Writes the summary, estimates, diagnostic tables, person measures, and
#' model-specific results as CSV. Plots are written as PNG and, optionally,
#' PDF, together with a plain-text analysis summary. For MFRM and EFRM fits,
#' item estimates and response-cell diagnostics are saved separately.
#' Computed tailored item shifts and externally weighted secondary person
#' measures can be retained with the fitted-model output; the latter include
#' their resolved weights.
#' For keyed fits with repeated person IDs, the export records why distractor
#' analysis is unavailable and retains the other model outputs.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param dir Output directory; created if absent. An existing directory must
#'   be empty; exports refuse to overwrite files from an earlier analysis.
#' @param formats Plot formats, any of \code{"png"} and \code{"pdf"}.
#' @param width,height Plot size in inches.
#' @param dpi PNG resolution.
#' @param dif Optional \code{\link{dif_anova}} result to export as computed
#'   --- an application analysis carries the DIF model the analyst chose,
#'   which a default recomputation would silently replace. \code{NULL}
#'   computes the default when the fit carries person factors.
#' @param bootstrap Optional \code{\link{fit_bootstrap}} result from this fit. Its item and
#'   person tables, or pair, object and judge tables for paired comparisons,
#'   join the export with the whole-test readings.
#' @param dif_bootstrap Optional \code{\link{dif_bootstrap}} sensitivity
#'   analysis from this fit and DIF specification.
#' @param dimensionality Optional computed \code{\link{plot_scree}} or
#'   \code{\link{btl_dimensionality}} result from this fit. Supplying it keeps
#'   the exported table and figure identical to the analysis already run.
#' @param subtest Optional computed \code{\link{dimensionality_test}} result
#'   from this fit. Supplying it keeps the nominated item split and bootstrap
#'   calibration used in the analysis.
#' @param invariance Optional computed \code{\link{frame_invariance}} result
#'   from an EFRM fit.
#' @param tailored Optional computed \code{\link{tailored_analysis}} result
#'   from this ordinary dichotomous fit.
#' @param person_weights Optional table returned by
#'   \code{\link{weighted_person_estimates}}. An explicit table takes
#'   precedence over a compatible result retained by the application.
#' @param item_plots Also write the per-item plot set (one ICC, category curve,
#'   threshold curve, and frequency chart per item).
#' @return Invisibly, the vector of files written.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(150 * 6, 1, plogis(outer(rnorm(150), d, "-"))), 150, 6)
#' colnames(X) <- paste0("I", 1:6)
#' out <- tempfile("rasch-out-")
#' save_outputs(rasch(X), out, formats = "png", item_plots = FALSE, dpi = 96)
#' @export
save_outputs <- function(fit, dir, formats = c("png", "pdf"), width = 9,
                         height = 6, dpi = 300, item_plots = TRUE,
                         dif = NULL, bootstrap = NULL,
                         dif_bootstrap = NULL, dimensionality = NULL,
                         invariance = NULL, subtest = NULL,
                         tailored = NULL, person_weights = NULL) {
  .check_export_fit(fit)
  person_weight_result <- .resolve_report_person_weights(fit, person_weights)
  formats <- match.arg(formats, c("png", "pdf"), several.ok = TRUE)
  .validate_boot_dif_result(dif, fit)
  .validate_fit_bootstrap(bootstrap, fit)
  if (!is.null(dif_bootstrap) && is.null(dif))
    stop("`dif` must accompany `dif_bootstrap` so the primary and sensitivity analyses use the same specification")
  .validate_dif_bootstrap(dif_bootstrap, fit, dif)
  if (inherits(fit, "rasch_btl")) {
    if (!is.null(subtest))
      stop("`subtest` is not available for a paired-comparison fit")
    .validate_btl_dimensionality(dimensionality, fit)
  } else {
    .validate_scree_result(dimensionality, fit)
    .validate_dimensionality_test(subtest, fit)
  }
  .validate_frame_invariance(invariance, fit)
  .validate_tailored_result(tailored, fit)
  # everything is checked before a directory is made or a table written: a
  # bad plot size otherwise leaves a populated folder that reads as a
  # complete export but carries no plots
  .check_export_dir(dir)
  .check_device_size(width, height, dpi)
  .check_flag(item_plots, "item_plots")
  if (inherits(fit, "rasch_btl"))
    return(.save_btl_outputs(fit, dir, formats, width, height, dpi,
                             object_plots = item_plots,
                             bootstrap = bootstrap, dif = dif,
                             dif_bootstrap = dif_bootstrap,
                             dimensionality = dimensionality))
  distractors <- .distractors_or_note(fit)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tdir <- file.path(dir, "tables"); pdir <- file.path(dir, "plots")
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  idir <- file.path(pdir, "items")
  cdir <- if (structural) file.path(pdir, "response_cells") else idir
  dir.create(tdir, showWarnings = FALSE)
  dir.create(pdir, showWarnings = FALSE)
  if (item_plots) {
    dir.create(cdir, showWarnings = FALSE)
    if (inherits(fit, "rasch_efrm")) dir.create(idir, showWarnings = FALSE)
  }
  files <- character(0)
  wtab <- function(d, name) {
    path <- file.path(tdir, paste0(name, ".csv"))
    .write_csv_plain(d, path)
    files <<- c(files, path)
  }

  dt <- if (!is.null(subtest)) subtest else .dimensionality_or_note(fit)
  if (is.null(dt$note)) {
    wtab(data.frame(
      split = dt$split, items_positive = paste(dt$items_positive, collapse = ";"),
      items_negative = paste(dt$items_negative, collapse = ";"),
      prop_significant = dt$prop_significant, ci_lower = dt$ci[1],
      ci_upper = dt$ci[2], n = dt$n,
      n_excluded_extreme = dt$n_excluded_extreme,
      multidimensional = dt$multidimensional,
      binomial_multidimensional = dt$binomial_multidimensional,
      verdict_method = dt$verdict_method,
      p_boot = dt$p_boot %||% NA_real_,
      bootstrap_resolution = dt$bootstrap_resolution %||% NA_real_,
      B_requested = dt$bootstrap$B %||% NA_integer_,
      B_used = dt$bootstrap$B_used %||% NA_integer_,
      B_nonconverged = dt$bootstrap$B_nonconverged %||% NA_integer_,
      B_errors = dt$bootstrap$B_errors %||% NA_integer_,
      stringsAsFactors = FALSE),
      "unidimensionality_t_test")
    if (!is.null(dt$bootstrap))
      wtab(data.frame(replicate = seq_along(dt$bootstrap$null),
                      prop_significant = dt$bootstrap$null),
           "unidimensionality_t_test_bootstrap")
  } else {
    wtab(data.frame(note = dt$note), "unidimensionality_t_test")
  }

  # --- tables ---------------------------------------------------------------
  wtab(fit_summary_table(fit), "fit_summary")
  wtab(targeting_table(fit), "targeting_and_reliability")
  wtab(fit$items, if (structural) "response_cell_statistics" else
    "item_statistics")
  wtab(fit$item_anova, if (structural) "response_cell_anova_fit" else
    "item_anova_fit")
  thr <- fit$thresholds
  thr$item <- fit$items$item[thr$item]
  wtab(thr[, c("item", "k", "tau", "se")], if (structural)
    "response_cell_thresholds" else "thresholds")
  if (!is.null(fit$est$components)) wtab(fit$est$components, "principal_components")
  if (inherits(fit, "rasch_explanatory")) {
    wtab(explanatory_test(fit), "explanatory_model_comparison")
    wtab(fit$est$coefficients, "explanatory_predictor_effects")
    wtab(explanatory_diagnostics(fit), "explanatory_fixed_departure_diagnostics")
    if (nrow(fit$explanatory$relaxations))
      wtab(fit$explanatory$relaxations, "explanatory_fixed_departures")
  }
  wtab(fit$person, "person_estimates")
  if (!is.null(person_weight_result)) {
    wtab(person_weight_result$table, "person_estimates_externally_weighted")
    wtab(attr(person_weight_result$table, "weighting", exact = TRUE),
         "person_estimate_external_weights")
  }
  if (!is.null(fit$score_table)) wtab(score_table(fit), "score_to_measure")
  ctt <- tryCatch(ctt_table(fit), error = function(e) NULL)
  if (!is.null(ctt)) wtab(ctt$table, "traditional_statistics")
  cd_all <- do.call(rbind, lapply(fit$items$item, function(it) {
    cd <- tryCatch(chisq_detail(fit, it), error = function(e) NULL)
    if (is.null(cd)) return(NULL)
    cbind(item = cd$item, cd$intervals)
  }))
  wtab(cd_all, if (structural) "response_cell_chisq_class_interval_detail"
       else "chisq_class_interval_detail")
  rc <- residual_correlations(fit)
  wtab(data.frame(item = rownames(rc$matrix), round(rc$matrix, 4),
                  check.names = FALSE), if (structural)
                    "response_cell_residual_correlations" else
                      "residual_correlations")
  wtab(rc$pairs, if (structural) "response_cell_q3_statistics" else
    "q3_statistics")
  if (nrow(rc$flagged)) wtab(rc$flagged, if (structural)
    "response_cell_local_dependence_flagged" else "local_dependence_flagged")
  # residual PCA refuses structurally disjoint designs (extended-frame
  # groups, facet cells) -- record the reason instead of failing the export
  pc <- tryCatch(residual_pca(fit), error = function(e)
    structure(list(msg = conditionMessage(e)), class = "rr_pca_refusal"))
  if (inherits(pc, "rr_pca_refusal")) {
    wtab(data.frame(note = pc$msg), "pca_loadings")
  } else {
    wtab(pc$loadings_matrix, "pca_loadings")
    # a reference that cannot be completed leaves the observed eigenvalues
    # in the table and the scree plot out of the export; an omitted display
    # must not pass silently for an export the user will treat as complete
    scree <- if (!is.null(dimensionality)) dimensionality else
      tryCatch(.scree_analysis(fit, parallel = !structural, seed = 1L),
               error = function(e) {
                 warning("plot 'scree' could not be drawn: ",
                         conditionMessage(e), call. = FALSE)
                 structure(list(msg = conditionMessage(e)),
                           class = "rr_pca_refusal")
               })
    wtab(if (inherits(scree, "rr_pca_refusal")) pc$eigen_table else scree,
         "residual_eigenvalues")
  }
  cf <- do.call(rbind, lapply(fit$thresholds_diag, function(d)
    data.frame(item = d$item, category = seq_along(d$category_counts) - 1L,
               count = d$category_counts)))
  wtab(cf, if (structural) "response_cell_category_frequencies" else
    "category_frequencies")
  whole_item_design <- .classical_design_applicable(fit)
  if (whole_item_design && all(fit$m == 1L)) {
    gt <- guttman_table(fit)
    wtab(data.frame(id = rownames(gt$matrix), gt$matrix, check.names = FALSE),
         "guttman_ordered_responses")
  }
  if (!is.null(distractors))
    wtab(if (is.null(distractors$note)) distractors$table else
      data.frame(note = distractors$note), "distractor_analysis")
  if (inherits(fit, "rasch_mfrm")) {
    wtab(fit$item_effects, "item_effects")
    wtab(fit$item_thresholds, "item_structural_thresholds")
    fsafe <- .rr_safe_stem(fit$facet_spec)
    for (j in seq_along(fit$facet_spec))
      wtab(fit$facet_effects[[fit$facet_spec[j]]],
           paste0("facet_", fsafe[j]))
    if (!is.null(fit$interaction_effects)) {
      wtab(fit$interaction_test, "interaction_omnibus_test")
      wtab(fit$interaction_effects, "item_by_facet_interactions")
    }
  }
  if (inherits(fit, "rasch_efrm")) {
    wtab(fit$frames, "frames")
    x <- fit$efrm_vs_rasch
    wtab(data.frame(
      model = c("Equal group units", "Group-dependent units"),
      loglik = c(x$ll_equal, x$ll_efrm),
      unit_parameters = c(0L, x$extra_parameters),
      two_delta_loglik = c(NA_real_, x$two_delta_ll)),
      "frame_model_comparison")
    wtab(x$unit_omnibus, "unit_omnibus_tests")
    wtab(x$unit_tests, "unit_contrasts")
    # the fit holds each item at one location across frames, so the units are
    # only as good as that assumption; save the test of it beside them
    inv <- if (!is.null(invariance)) invariance else
      tryCatch(frame_invariance(fit), error = function(e) e)
    if (inherits(inv, "error")) {
      wtab(data.frame(note = paste("Frame invariance was not available:",
                                   conditionMessage(inv))),
           "frame_invariance_summary")
    } else {
      wtab(inv$summary, "frame_invariance_summary")
      wtab(inv$locations, "frame_invariance_locations")
      wtab(inv$discrimination, "frame_invariance_discrimination")
      if (identical(inv$se_method, "bootstrap"))
        wtab(data.frame(
          family_n = inv$family_n,
          requested = inv$boot_reps,
          usable = inv$boot_reps_used,
          nonconverged = inv$boot_reps_nonconverged,
          other_failures = inv$boot_reps_errors,
          minimum_usable = inv$boot_minimum_usable),
          "frame_invariance_bootstrap_accounting")
      if (!is.null(inv$excluded) && nrow(inv$excluded))
        wtab(inv$excluded, "frame_invariance_excluded")
    }
    wtab(fit$phi_table, "group_units_phi")
    wtab(fit$alpha_table, "set_units_alpha")
    wtab(fit$set_table, "set_locations")
    wtab(fit$item_arbitrary, "items_common_unit")
    wtab(fit$thresholds_arbitrary, "thresholds_common_unit")
    wtab(fit$score_curves, "score_curves")
  }
  dif_notes <- character(0)
  if (!is.null(dif) || !is.null(fit$factors)) {
    # an analysis exported from the application carries the DIF model the
    # analyst actually chose; recomputing at defaults would export a
    # different analysis than the one on screen
    da <- if (!is.null(dif)) dif
          else tryCatch(dif_anova(fit), error = function(e) NULL)
    if (!is.null(da)) {
      wtab(da$summary, "dif_anova")
      wtab(da$terms, "dif_anova_terms")
      # a blank row in the exported table is a test the design could not
      # estimate, not an absence of DIF: the summary carries the reason
      dif_notes <- as.character(da$notes)
    }
  }
  if (!is.null(bootstrap)) {
    wtab(bootstrap$items, "bootstrap_fit")
    if (!is.null(bootstrap$persons))
      wtab(bootstrap$persons, "bootstrap_person_fit")
    wtab(data.frame(quantity = names(unlist(bootstrap$total)),
                    value = unlist(bootstrap$total)), "bootstrap_total")
  }
  if (!is.null(dif_bootstrap)) {
    wtab(dif_bootstrap$summary, "dif_conditional_bootstrap")
    wtab(dif_bootstrap$terms, "dif_conditional_bootstrap_terms")
    wtab(data.frame(
      quantity = c("requested", "usable", "failed", "nonconverged",
                   "other_errors", "family_size"),
      value = c(dif_bootstrap$B, dif_bootstrap$B_used,
                dif_bootstrap$B_failed, dif_bootstrap$B_nonconverged,
                dif_bootstrap$B_errors, dif_bootstrap$family_n)),
      "dif_conditional_bootstrap_accounting")
  }
  if (!is.null(tailored)) {
    wtab(tailored$table, "tailored_item_shifts")
    wtab(data.frame(
      algorithm = tailored$algorithm,
      chance = tailored$chance,
      responses_removed = tailored$n_removed,
      anchor_items = paste(tailored$anchor_items, collapse = ";"),
      anchor_selection = if (is.null(tailored$anchor_items_requested))
        "automatic" else "supplied",
      anchor_items_requested = if (is.null(tailored$anchor_items_requested))
        NA_character_ else paste(tailored$anchor_items_requested,
                                 collapse = ";"),
      se_method = tailored$se_method,
      seed = tailored$seed %||% NA_integer_,
      requested = tailored$boot_reps,
      usable = tailored$boot_reps_used,
      nonconverged = tailored$boot_reps_nonconverged,
      other_failures = tailored$boot_reps_errors,
      minimum_usable = tailored$boot_minimum_usable,
      fit_fingerprint = tailored$fit_signature$fingerprint,
      result_signature = tailored$result_signature,
      stringsAsFactors = FALSE), "tailored_analysis")
  }
  if (any(fit$person$extreme)) {
    pe <- tryCatch(person_extrapolated(fit), error = function(e) NULL)
    if (!is.null(pe)) wtab(pe, "person_estimates_extrapolated")
  }

  # --- summary ---------------------------------------------------------------
  spath <- file.path(dir, "summary.txt")
  con <- file(spath, "w")
  sink(con); on.exit({ sink(); close(con) }, add = TRUE)
  summary(fit)
  if (is.null(dt$note)) {
    resolution_limited <- grepl(
      "resolution insufficient", dt$verdict_method %||% "", fixed = TRUE)
    verdict <- if (isTRUE(dt$multidimensional)) "evidence against unidimensionality" else
      if (identical(dt$multidimensional, FALSE)) "consistent with one dimension" else
        if (resolution_limited)
          "inferential verdict withheld because bootstrap resolution is insufficient" else
          "inferential verdict withheld for the data-driven split"
    cat(sprintf("\nUnidimensionality t-test: %.1f%% significant (Clopper-Pearson 95%% CI %.1f%% to %.1f%%), %s\n",
                100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2],
                verdict))
    if (!is.null(dt$p_boot))
      cat(sprintf("Parametric-bootstrap p = %s (%d of %d replicates used)\n",
                  .fmt_p(dt$p_boot), dt$bootstrap$B_used, dt$bootstrap$B))
    if (!is.null(dt$caution)) cat("Caution:", dt$caution, "\n")
    if (is.na(dt$multidimensional)) {
      if (resolution_limited)
        cat(sprintf(paste0(
          "Note: the smallest attainable bootstrap p is %.3f; increase B ",
          "for an inferential verdict at alpha %.3f\n"),
          dt$bootstrap_resolution, dt$alpha))
      else
        cat("Note: the item split was chosen from the residuals; use a content split or a bootstrap calibration for an inferential verdict\n")
    }
  } else cat("\nUnidimensionality t-test:", dt$note, "\n")
  cat(sprintf("Average residual correlation: %.3f; binary Q3 flags withheld (no universal critical value)\n",
              rc$average))
  if (!is.null(ctt))
    cat(sprintf("Traditional statistics (complete cases n = %d): raw mean %.2f, SD %.2f, alpha %.3f, SEM %.2f\n",
                ctt$n, ctt$mean, ctt$sd, ctt$alpha, ctt$sem))
  if (!is.null(distractors$note))
    cat("Distractor analysis unavailable:", distractors$note, "\n")
  if (length(dif_notes))
    cat("DIF notes:", paste(dif_notes, collapse = "; "), "\n")
  if (!is.null(tailored)) {
    cat(sprintf(paste0("Tailored analysis: %d responses below chance %.2f ",
                       "set to missing; origin anchored by %s"),
                tailored$n_removed, tailored$chance,
                paste(tailored$anchor_items, collapse = ", ")))
    if (identical(tailored$se_method, "bootstrap"))
      cat(sprintf("; %d of %d person-bootstrap replicates usable\n",
                  tailored$boot_reps_used, tailored$boot_reps))
    else cat("; item shifts are descriptive\n")
  }
  sink(); close(con); on.exit()
  files <- c(files, spath)

  # --- test-level plots --------------------------------------------------------
  sp <- function(f, stem) files <<- c(files,
    .rr_save_plot(f, stem, pdir, formats, width, height, dpi))
  sp(function() plot_pimap(fit), if (structural)
    "person_calibration_distribution" else "person_item_distribution")
  sp(function() plot_wright(fit), "wright_map")
  sp(function() plot_threshold_map(fit), if (structural)
    "calibration_threshold_map" else "threshold_map")
  sp(function() plot_tcc(fit), "test_characteristic_curve")
  sp(function() plot_tif(fit), "test_information")
  sp(function() plot_item_map(fit), if (structural)
    "response_cell_fit_map" else "item_fit_map")
  sp(function() plot_person_fit(fit), "person_fit")
  sp(function() plot_resid_cor(fit), "residual_correlations")
  sp(function() plot_pca(fit), "pca_loadings")
  # A structural fit has an observed residual decomposition when its cells
  # overlap, but no valid full-refit parallel reference over the mutually
  # exclusive virtual design. Export the observed scree deliberately instead
  # of attempting the unavailable reference and warning from a normal export.
  if (exists("scree", inherits = FALSE) &&
      !inherits(scree, "rr_pca_refusal"))
    sp(function() plot_scree(fit, result = scree), "scree")
  if (whole_item_design && all(fit$m == 1L))
    sp(function() plot_guttman(fit), "guttman_scalogram")
  sp(function() plot_resid_dist(fit, "items"), if (structural)
    "response_cell_residual_distribution" else "item_residual_distribution")
  sp(function() plot_resid_dist(fit, "persons"), "person_residual_distribution")
  if (inherits(fit, "rasch_mfrm")) {
    fsafe <- .rr_safe_stem(fit$facet_spec)
    for (j in seq_along(fit$facet_spec)) local({
      f_ <- fit$facet_spec[j]; s_ <- fsafe[j]
      sp(function() plot_facets(fit, f_),
         paste0("facet_severities_", s_))
    })
  }
  if (inherits(fit, "rasch_efrm")) {
    sp(function() plot_frames(fit), "frame_units")
    if (item_plots) {
      vit <- unique(fit$virtual_map$item)
      vsafe <- .rr_safe_stem(vit)
      for (j in seq_along(vit)) local({
        it_ <- vit[j]; s_ <- vsafe[j]
        files <<- c(files, .rr_save_plot(function() plot_icc_frames(fit, it_),
          paste0(s_, "_icc_frames"), idir, formats, width, height, dpi))
      })
    }
  }

  # --- per-item plots ------------------------------------------------------------
  if (item_plots && !is.null(distractors$table)) {
    mcit <- colnames(fit$mc$raw)
    mcsafe <- .rr_safe_stem(mcit)
    for (j in seq_along(mcit)) local({
      it_ <- mcit[j]; s_ <- mcsafe[j]
      files <<- c(files, .rr_save_plot(function() plot_distractors(fit, it_),
        paste0(s_, "_options"), idir, formats, width, height, dpi))
    })
  }
  if (item_plots) {
    safe_items <- .rr_safe_stem(fit$items$item)
    for (j in seq_along(fit$items$item)) {
      it <- fit$items$item[j]
      safe <- safe_items[j]
      files <- c(files,
        .rr_save_plot(function() plot_icc(fit, it),
                      paste0(safe, "_icc"), cdir, formats, width, height, dpi),
        .rr_save_plot(function() plot_ccc(fit, it),
                      paste0(safe, "_categories"), cdir, formats, width, height, dpi),
        .rr_save_plot(function() plot_threshold_prob(fit, it),
                      paste0(safe, "_thresholds"), cdir, formats, width, height, dpi),
        .rr_save_plot(function() plot_catfreq(fit, it),
                      paste0(safe, "_frequencies"), cdir, formats, width, height, dpi))
    }
  }
  invisible(files)
}

# base-R base64 encoder (RFC 4648), used to embed plot images in the
# self-contained HTML report without adding dependencies
.b64 <- function(path) {
  raw <- readBin(path, "raw", file.info(path)$size)
  alphabet <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", "")[[1]]
  n <- length(raw)
  pad <- (3 - n %% 3) %% 3
  raw <- c(raw, as.raw(rep(0, pad)))
  m <- matrix(as.integer(raw), nrow = 3)
  b1 <- m[1, ] %/% 4
  b2 <- (m[1, ] %% 4) * 16 + m[2, ] %/% 16
  b3 <- (m[2, ] %% 16) * 4 + m[3, ] %/% 64
  b4 <- m[3, ] %% 64
  out <- alphabet[rbind(b1, b2, b3, b4) + 1L]
  if (pad > 0) out[(length(out) - pad + 1):length(out)] <- "="
  paste(out, collapse = "")
}

.report_css <- "
  body { font-family: -apple-system, 'Segoe UI', Roboto, Helvetica, Arial,
         sans-serif; color: #0f172a; margin: 0; background: #f8fafc; }
  .wrap { max-width: 980px; margin: 0 auto; padding: 2.5rem 1.5rem 4rem; }
  h1 { font-size: 1.6rem; margin: 0 0 .25rem; }
  h2 { font-size: 1.15rem; margin: 2.2rem 0 .6rem; padding-bottom: .3rem;
       border-bottom: 2px solid #e2e8f0; }
  .meta { color: #64748b; font-size: .85rem; margin-bottom: 1.5rem; }
  .note { color: #64748b; font-size: .82rem; margin: .3rem 0 .8rem; }
  .table-wrap { width: 100%; overflow-x: auto; margin: .4rem 0 1rem;
                -webkit-overflow-scrolling: touch; }
  table { border-collapse: collapse; width: 100%; font-size: .82rem;
          background: #fff; }
  th { text-align: left; font-weight: 600; border-bottom: 2px solid #cbd5e1;
       padding: .35rem .55rem; white-space: nowrap; }
  td { border-bottom: 1px solid #eef2f7; padding: .3rem .55rem; }
  td.num { text-align: right; font-variant-numeric: tabular-nums; }
  tr:hover td { background: #f1f5f9; }
  img { max-width: 100%; background: #fff; border: 1px solid #e2e8f0;
        border-radius: 8px; margin: .4rem 0 1rem; }
  .flag { color: #dc2626; font-weight: 600; }
  .chip { display: inline-block; background: #eff6ff; color: #1d4ed8;
          border-radius: 999px; padding: .1rem .6rem; font-size: .78rem;
          margin-right: .35rem; }
  @media (max-width: 620px) {
    .wrap { padding: 1.25rem .75rem 2.5rem; }
    h1 { font-size: 1.35rem; }
    table { min-width: 620px; }
  }
"

.html_escape <- function(x)
  gsub(">", "&gt;", gsub("<", "&lt;", gsub("&", "&amp;", as.character(x))))

.html_table <- function(d, digits = 3, max_rows = 500) {
  if (is.null(d) || !nrow(d)) return("")
  trunc_note <- ""
  if (nrow(d) > max_rows) {
    trunc_note <- sprintf("<p class='note'>Showing the first %d of %d rows.</p>",
                          max_rows, nrow(d))
    d <- d[seq_len(max_rows), , drop = FALSE]
  }
  esc <- .html_escape
  # drop all-FALSE logical flag columns and constant 'max' columns
  drop <- vapply(seq_along(d), function(j)
    (is.logical(d[[j]]) && !any(d[[j]], na.rm = TRUE)) ||
    (names(d)[j] == "max" && length(unique(d[[j]])) == 1L), TRUE)
  d <- d[, !drop, drop = FALSE]
  if (!ncol(d)) return("")
  fd <- .fmt_df(d, digits)
  num <- vapply(d, is.numeric, TRUE)
  cells <- vapply(seq_len(ncol(d)), function(j) {
    v <- fd[[j]]
    if (num[j]) gsub("<", "&lt;", v) else esc(v)
  }, character(nrow(d)))
  if (is.null(dim(cells))) cells <- matrix(cells, nrow = 1)
  head_html <- paste0("<th>", esc(names(d)), "</th>", collapse = "")
  body_html <- paste(vapply(seq_len(nrow(d)), function(i) {
    paste0("<tr>", paste0("<td", ifelse(num, " class='num'", ""), ">",
                          cells[i, ], "</td>", collapse = ""), "</tr>")
  }, ""), collapse = "\n")
  paste0(trunc_note, "<div class='table-wrap'><table><thead><tr>",
         head_html, "</tr></thead><tbody>", body_html,
         "</tbody></table></div>")
}

#' Write a self-contained HTML report of a Rasch analysis
#'
#' Writes one HTML file containing the summary statistics, diagnostic tables,
#' and test-level plots. Images and styles are embedded in the file. Computed
#' tailored item shifts and externally weighted secondary person measures can
#' be included with the fitted-model results.
#' For keyed fits with repeated person IDs, the report explains why distractor
#' analysis is unavailable; the other model outputs remain available.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param file Path of the HTML file to write.
#' @param title Report title.
#' @param dpi Resolution of the embedded plots.
#' @param dif,bootstrap Optional computed \code{\link{dif_anova}} and
#'   \code{\link{fit_bootstrap}} results from this fit, exported as run; the DIF table is
#'   otherwise recomputed at defaults when the fit carries person factors.
#' @param dif_bootstrap Optional \code{\link{dif_bootstrap}} sensitivity
#'   analysis from this fit and DIF specification.
#' @param dimensionality Optional computed \code{\link{plot_scree}} result
#'   from this fit. Supplying it keeps the table and figure identical to the
#'   analysis already run.
#' @param subtest Optional computed \code{\link{dimensionality_test}} result
#'   from this fit. Supplying it keeps the item split and bootstrap calibration
#'   used in the analysis.
#' @param invariance Optional computed \code{\link{frame_invariance}} result
#'   from an EFRM fit.
#' @param tailored Optional computed \code{\link{tailored_analysis}} result
#'   from this ordinary dichotomous fit.
#' @param person_weights Optional table returned by
#'   \code{\link{weighted_person_estimates}}. An explicit table takes
#'   precedence over a compatible result retained by the application.
#' @return Invisibly, \code{file}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(150 * 6, 1, plogis(outer(rnorm(150), d, "-"))), 150, 6)
#' colnames(X) <- paste0("I", 1:6)
#' out <- file.path(tempdir(), "report.html")
#' report_html(rasch(X), out)
#' @export
report_html <- function(fit, file, title = "Rasch measurement analysis",
                        dpi = 150, dif = NULL, bootstrap = NULL,
                        dif_bootstrap = NULL, dimensionality = NULL,
                        invariance = NULL, subtest = NULL,
                        tailored = NULL, person_weights = NULL) {
  # a vector title would be pasted into as many documents as it has entries,
  # and a non-positive dpi can take the graphics device down with the
  # session rather than raising a catchable error
  .check_export_fit(fit, btl = FALSE)
  person_weight_result <- .resolve_report_person_weights(fit, person_weights)
  .check_out_path(file, "file")
  .validate_boot_dif_result(dif, fit)
  .validate_fit_bootstrap(bootstrap, fit)
  if (!is.null(dif_bootstrap) && is.null(dif))
    stop("`dif` must accompany `dif_bootstrap` so the primary and sensitivity analyses use the same specification")
  .validate_dif_bootstrap(dif_bootstrap, fit, dif)
  .validate_scree_result(dimensionality, fit)
  .validate_dimensionality_test(subtest, fit)
  .validate_frame_invariance(invariance, fit)
  .validate_tailored_result(tailored, fit)
  if (length(title) != 1L || !is.character(title) || !is.null(dim(title)) ||
      !is.null(oldClass(title)) || is.na(title))
    stop("`title` must be one non-missing title")
  .check_pos_num(dpi, "dpi")
  distractors <- .distractors_or_note(fit)
  tmp <- tempfile("rmtplots"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  shot <- function(f, name, w = 9, h = 5.4) {
    path <- file.path(tmp, paste0(name, ".png"))
    before <- grDevices::dev.cur()
    ok <- tryCatch({
      png(path, width = w, height = h, units = "in", res = dpi)
      f(); dev.off(); TRUE
    }, error = function(e) {
      # close only a device this call opened, never the caller's
      if (!identical(grDevices::dev.cur(), before))
        try(dev.off(), silent = TRUE)
      warning("report plot '", name, "' could not be drawn: ",
              conditionMessage(e), call. = FALSE)
      FALSE })
    if (!ok) return("")
    sprintf("<img src='data:image/png;base64,%s' alt='%s'/>", .b64(path), name)
  }
  s <- function(...) paste0(...)
  # data-derived text (title, notes, item names) must be escaped before it
  # is pasted into markup: an item named "A<b>" must render as text
  esc <- function(x) gsub(">", "&gt;", gsub("<", "&lt;",
                          gsub("&", "&amp;", as.character(x))))
  title <- esc(title)
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  # a refusal (structurally disjoint columns) or a failed reference leaves
  # the scree out of the report; say why, in the document and as a warning,
  # rather than drop the section without a word
  scree_refusal <- NULL
  scree <- if (!is.null(dimensionality)) dimensionality else
    tryCatch(.scree_analysis(fit, parallel = !structural, seed = 1L),
             error = function(e) {
               scree_refusal <<- conditionMessage(e)
               warning("report plot 'scree' could not be drawn: ",
                       conditionMessage(e), call. = FALSE)
               NULL
             })
  alpha_design <- .classical_design_applicable(fit)
  item_count <- if (inherits(fit, "rasch_mfrm")) nrow(fit$item_effects)
    else if (inherits(fit, "rasch_efrm")) nrow(fit$item_arbitrary)
    else ncol(fit$X)
  model_name <- if (inherits(fit, "rasch_explanatory"))
    fit$explanatory_model else fit$model
  chips <- s("<span class='chip'>", esc(model_name), "</span>",
             "<span class='chip'>", nrow(fit$X), " persons</span>",
             "<span class='chip'>", item_count, " items</span>",
             if (structural) s("<span class='chip'>", ncol(fit$X),
                               " response cells</span>") else "",
             sprintf("<span class='chip'>PSI %.3f</span>", fit$psi$PSI),
             if (alpha_design && is.finite(fit$alpha$alpha))
               sprintf("<span class='chip'>alpha %.3f</span>",
                       fit$alpha$alpha) else "")
  calibration_unit <- if (structural) "response-cell" else "item"
  fit_unit <- if (structural) "Response-cell" else "Item"
  separation_unit <- if (structural) "response-cell" else "item"
  separation_quality <- fit$separation_quality %||% fit$power_of_fit %||%
    .separation_quality(fit$psi$PSI)
  total_p <- if (length(fit$total_chisq_p) == 1L &&
                 is.finite(fit$total_chisq_p))
    .fmt_p(fit$total_chisq_p) else "unavailable"
  summ <- s(
    sprintf("<p>%s %s in %d iterations. ",
            .html_escape(.estimation_label(fit)),
            if (isTRUE(fit$est$converged)) "converged" else "did <b>not</b> converge",
            fit$est$iterations),
    sprintf("Approximate asymptotic total %s-trait chi-square %.2f on %d df (p = %s). ",
            calibration_unit, fit$total_chisq, fit$total_df,
            total_p),
    sprintf("%s fit residual mean %.2f, SD %.2f; person fit residual mean %.2f, SD %.2f. ",
            fit_unit, fit$item_fit_summary$mean, fit$item_fit_summary$sd,
            fit$person_fit_summary$mean, fit$person_fit_summary$sd),
    sprintf("PSI %.3f (%.3f without extremes); %s separation reliability %.3f; person separation quality: %s. Fit-test sensitivity also depends on sample size, targeting and the departure being tested.</p>",
            fit$psi$PSI, fit$psi_noext$PSI, separation_unit,
            fit$isi$PSI, separation_quality),
    if (length(fit$notes))
      s("<p class='note'>Notes: ", esc(paste(fit$notes, collapse = "; ")), "</p>")
    else "")
  rc <- residual_correlations(fit)
  dt <- if (!is.null(subtest)) subtest else .dimensionality_or_note(fit)
  dim_html <- if (is.null(dt$note)) {
    resolution_limited <- grepl(
      "resolution insufficient", dt$verdict_method %||% "", fixed = TRUE)
    verdict <- if (isTRUE(dt$multidimensional))
      "<span class='flag'>evidence against unidimensionality</span>" else
      if (identical(dt$multidimensional, FALSE)) "consistent with one dimension" else
        if (resolution_limited)
          "inferential verdict withheld because bootstrap resolution is insufficient" else
          "inferential verdict withheld for the data-driven split"
    paste0(sprintf("<p>%.1f%% of person subset t-tests significant (95%% CI %.1f-%.1f%%): %s.</p>",
            100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2], verdict),
           if (!is.null(dt$p_boot)) sprintf(
             "<p>Parametric-bootstrap p = %s (%d of %d replicates used).</p>",
             .fmt_p(dt$p_boot), dt$bootstrap$B_used, dt$bootstrap$B) else "",
           if (!is.null(dt$caution))
             sprintf("<p class='note'>%s</p>", esc(dt$caution)) else "",
           if (is.na(dt$multidimensional) && resolution_limited)
             sprintf(paste0(
               "<p class='note'>The smallest attainable bootstrap p is %.3f. ",
               "Increase B for an inferential verdict at alpha %.3f.</p>"),
               dt$bootstrap_resolution, dt$alpha) else
           if (is.na(dt$multidimensional))
             "<p class='note'>The item split was chosen from the residuals. Use a content split or a bootstrap calibration for an inferential verdict.</p>" else "")
  }
  else sprintf("<p class='note'>%s</p>", esc(dt$note))
  ctt <- tryCatch(ctt_table(fit), error = function(e) NULL)
  item_tab <- if (inherits(fit, "rasch_mfrm")) fit$item_effects else
    if (inherits(fit, "rasch_efrm")) fit$item_arbitrary else fit$items
  if (.has_repeated_residual_units(fit))
    for (nm in intersect(c("p", "p_adj", "p_bonf", "p_anova",
                           "p_anova_adj", "p_anova_bonf"), names(item_tab)))
      item_tab[[nm]][] <- NA_real_
  item_cols <- intersect(c("item", "set", "max", "location", "se", "n",
                           "fit_resid", "fit_resid_pooled", "infit_ms",
                           "outfit_ms", "chisq", "df", "p_adj", "weak"),
                         names(item_tab))
  common_thresholds <- if (inherits(fit, "rasch_mfrm")) fit$item_thresholds
    else if (inherits(fit, "rasch_efrm")) fit$thresholds_arbitrary
    else {
      th <- fit$thresholds
      th$item <- fit$items$item[th$item]
      th[, c("item", "k", "tau", "se")]
    }
  if ("delta" %in% names(common_thresholds) &&
      !"tau" %in% names(common_thresholds))
    names(common_thresholds)[names(common_thresholds) == "delta"] <- "tau"
  common_thresholds <- common_thresholds[, intersect(
    c("item", "set", "k", "tau", "se", "weak"),
    names(common_thresholds)), drop = FALSE]
  cell_thresholds <- NULL
  if (structural) {
    cell_thresholds <- fit$thresholds
    cell_thresholds$item <- fit$items$item[cell_thresholds$item]
    cell_thresholds <- cell_thresholds[, c("item", "k", "tau", "se"),
                                       drop = FALSE]
  }

  html <- s(
    "<!DOCTYPE html><html><head><meta charset='utf-8'/>",
    "<meta name='viewport' content='width=device-width, initial-scale=1'/>",
    "<title>", title, "</title><style>", .report_css, "</style></head><body>",
    "<div class='wrap'>",
    "<h1>", title, "</h1>",
    "<p class='meta'>", format(Sys.time(), "%Y-%m-%d %H:%M"),
    " &middot; rasch ", as.character(utils::packageVersion("rasch")), "</p>",
    "<p>", chips, "</p>",
    "<h2>Summary</h2>", summ,
    if (inherits(fit, "rasch_explanatory")) s(
      "<h2>Explanatory model</h2>",
      "<p class='note'>Formula: ", esc(fit$explanatory$formula_text), "</p>",
      "<h3>Comparison with free calibration</h3>",
      .html_table(explanatory_test(fit)),
      "<h3>Predictor effects</h3>",
      .html_table(fit$est$coefficients),
      "<h3>Fixed-departure diagnostics</h3>",
      .html_table(explanatory_diagnostics(fit)),
      if (nrow(fit$explanatory$relaxations)) s(
        "<h3>Fixed departures</h3>",
        .html_table(fit$explanatory$relaxations)) else "") else "",
    "<h2>Targeting</h2>",
    shot(function() plot_pimap(fit), "targeting"),
    shot(function() plot_wright(fit), "wright_map"),
    "<h2>", if (structural) "Common-scale item estimates" else
      "Item statistics", "</h2>",
    if ("p_adj" %in% item_cols)
      if (.has_repeated_residual_units(fit))
        "<p class='note'>Item-fit probabilities are withheld because person IDs repeat and the available references assume independent response rows. The fit statistics remain descriptive.</p>"
      else
        "<p class='note'>The p_adj column is an approximate asymptotic Holm probability. It treats estimated person locations as known; use the bootstrap fit statistics for calibrated inference.</p>"
    else "",
    if (structural)
      "<p class='note'>Item estimates on the common measurement scale.</p>"
    else "",
    .html_table(item_tab[, item_cols, drop = FALSE]),
    shot(function() plot_item_map(fit), "item_map"),
    "<h2>", if (structural) "Common-scale threshold estimates" else
      "Thresholds", "</h2>",
    .html_table(common_thresholds),
    if (structural) s("<h2>Response-cell fit</h2>",
      "<p class='note'>Observed item-by-frame or item-by-facet cells used in estimation. The p_adj column is an approximate asymptotic Holm probability.</p>",
      .html_table(fit$items[, intersect(
        c("item", "max", "location", "se", "fit_resid", "infit_ms",
          "outfit_ms", "chisq", "df", "p_adj"), names(fit$items)),
        drop = FALSE]),
      "<h2>Response-cell thresholds</h2>", .html_table(cell_thresholds)) else "",
    { dis <- names(which(vapply(fit$thresholds_diag, function(dd)
        !dd$ordered && length(dd$thresholds) > 1L, TRUE)))
      if (length(dis)) sprintf("<p class='flag'>Disordered %sthresholds: %s.</p>",
                               if (structural) "response-cell " else "",
                               esc(paste(dis, collapse = ", ")))
      else if (structural)
        "<p class='note'>All polytomous response cells have ordered thresholds.</p>"
      else "<p class='note'>All polytomous items have ordered thresholds.</p>" },
    shot(function() plot_threshold_map(fit), "threshold_map"),
    "<h2>Test characteristic and information</h2>",
    shot(function() plot_tcc(fit), "tcc"),
    shot(function() plot_tif(fit), "tif"),
    if (!is.null(fit$score_table)) s("<h2>Score to measure</h2>",
      .html_table(score_table(fit))) else "",
    "<h2>Fit residual distributions</h2>",
    shot(function() plot_resid_dist(fit, "items"), "resid_items"),
    shot(function() plot_resid_dist(fit, "persons"), "resid_persons"),
    "<h2>Dimensionality</h2>", dim_html,
    if (!is.null(scree_refusal))
      sprintf("<p class='note'>Scree plot unavailable: %s</p>",
              esc(scree_refusal)) else "",
    if (!is.null(scree)) .html_table(scree) else "",
    if (!is.null(scree))
      shot(function() plot_scree(fit, result = scree), "scree") else "",
    "<h2>Local dependence</h2>",
    sprintf(paste0("<p>Average residual correlation %.3f; binary Q3 flags ",
                   "withheld because there is no universal critical value.</p>"),
            rc$average),
    if (nrow(rc$flagged)) .html_table(rc$flagged) else "",
    shot(function() plot_resid_cor(fit), "residcor"),
    if (!is.null(ctt)) s("<h2>Classical companions</h2>",
      sprintf("<p class='note'>Complete cases n = %d; raw mean %.2f, SD %.2f; alpha %.3f; classical SEM %.2f.</p>",
              ctt$n, ctt$mean, ctt$sd, ctt$alpha, ctt$sem),
      .html_table(ctt$table)) else "",
    if (!is.null(bootstrap)) s("<h2>Bootstrap fit statistics</h2>",
      sprintf(paste0("<p class='note'>Parametric bootstrap null (%s scheme): ",
                     "%d of %d replicates; %d did not converge and %d otherwise ",
                     "failed. Adjustment is separate for each statistic and ",
                     "refers to the fitted global null.</p>"),
              bootstrap$theta, bootstrap$B_used, bootstrap$B,
              bootstrap$B_nonconverged, bootstrap$B_errors),
      .html_table(bootstrap$items[, intersect(c("item", "chisq",
                                     "chisq_p_boot_adj", "n_boot_chisq",
                                     "fit_resid", "fit_resid_p_boot_adj",
                                     "n_boot_fit_resid"),
                                   names(bootstrap$items))]),
      if (!is.null(bootstrap$persons))
        .html_table(bootstrap$persons[, intersect(c("id", "raw", "theta",
          "fit_resid", "fit_resid_p_boot_adj", "n_boot_fit_resid"),
          names(bootstrap$persons)), drop = FALSE]) else "") else "",
    if (!is.null(tailored)) s(
      "<h2>Tailored analysis for guessing</h2>",
      sprintf(paste0("<p class='note'>%d responses with fitted probability ",
                     "below chance %.2f were set to missing. The common ",
                     "origin uses: %s.</p>"), tailored$n_removed,
              tailored$chance, esc(paste(tailored$anchor_items,
                                         collapse = ", "))),
      if (identical(tailored$se_method, "bootstrap"))
        sprintf(paste0("<p class='note'>Person bootstrap: %d of %d ",
                       "replicates usable; probabilities are Holm-adjusted ",
                       "over items.</p>"), tailored$boot_reps_used,
                tailored$boot_reps)
      else
        "<p class='note'>Item shifts are descriptive; no bootstrap uncertainty was requested.</p>",
      .html_table(tailored$table)) else "",
    if (!is.null(dif) || !is.null(fit$factors)) {
      da <- if (!is.null(dif)) dif
            else tryCatch(dif_anova(fit), error = function(e) NULL)
      if (!is.null(da)) s("<h2>Differential item functioning</h2>",
        # a test the design could not estimate leaves a blank row in the
        # table and still counts in the adjusted family: say so beside it
        if (length(da$notes)) s("<p class='note'>Notes: ",
          esc(paste(da$notes, collapse = "; ")), "</p>") else "",
        .html_table(da$summary[, intersect(c("item", "term", "F_uniform",
                                     "p_uniform_adj", "eta2_uniform",
                                     "F_nonuniform", "p_nonuniform_adj",
                                     "eta2_nonuniform", "uniform_DIF",
                                     "nonuniform_DIF"), names(da$summary))]),
        if (!is.null(da$sizes) && nrow(da$sizes)) s(
          "<h3>DIF magnitude</h3>",
          .html_table(da$sizes[, intersect(
            c("item", "term", "level_a", "level_b", "difference", "se",
              "t", "z", "df", "p_adj", "practical"), names(da$sizes)),
            drop = FALSE]))
        else "") else ""
    } else "",
    if (!is.null(dif_bootstrap)) s(
      "<h3>Bootstrap sensitivity analysis</h3>",
      sprintf(paste0("<p class='note'>%d of %d complete-family replicates; ",
                     "%d did not converge and %d otherwise failed. ",
                     "Bootstrap decisions use single-step minimum-p ",
                     "familywise probabilities under the fitted global ",
                     "invariant null; they do not provide strong control when ",
                     "another member has DIF. The adjusted residual ANOVA ",
                     "above remains the primary DIF analysis. Null generator: ",
                     "%s.</p>"),
              dif_bootstrap$B_used, dif_bootstrap$B,
              dif_bootstrap$B_nonconverged, dif_bootstrap$B_errors,
              dif_bootstrap$null_method %||% "conditional"),
      .html_table(dif_bootstrap$summary[, intersect(
        c("item", "term", "p_uniform_boot_adj", "uniform_DIF_boot",
          "p_nonuniform_boot_adj", "nonuniform_DIF_boot"),
        names(dif_bootstrap$summary)), drop = FALSE])) else "",
    if (!is.null(distractors)) s("<h2>Distractor analysis</h2>",
      if (!is.null(distractors$note)) s("<p class='note'>Not available: ",
        esc(distractors$note), "</p>") else s(
        "<p class='note'>Locations use the rest measure; a distractor whose takers are abler than the keyed option's flags a possible miskey.</p>",
        .html_table(distractors$table))) else "",
    if (inherits(fit, "rasch_mfrm")) s("<h2>Facet severities</h2>",
      paste(vapply(fit$facet_spec, function(f) s("<h3>", esc(f), "</h3>",
        .html_table(fit$facet_effects[[f]][, intersect(c("level", "severity",
          "se", "n", "fit_resid", "fit_resid_pooled", "infit_ms",
          "outfit_ms"), names(fit$facet_effects[[f]])), drop = FALSE])), ""),
        collapse = ""),
      if (!is.null(fit$interaction_test)) s(
        "<h2>Item-by-facet interaction</h2>",
        "<p class='note'>The omnibus test assesses the complete interaction. Cell comparisons are Holm-adjusted follow-ups.</p>",
        .html_table(fit$interaction_test),
        .html_table(fit$interaction_effects[, intersect(
          c("item", "level", "gamma", "se", "t", "df", "z", "p_adj",
            "significant"),
          names(fit$interaction_effects)), drop = FALSE])) else "") else "",
    if (inherits(fit, "rasch_efrm")) {
      x <- fit$efrm_vs_rasch
      s("<h2>Frame model comparison</h2>",
        "<p class='note'>Within-frame thresholds and group units use pairwise conditional calibration. Item-set units use common persons and a finite-grid semiparametric link with a separate nuisance distribution for each observed person group.</p>",
        "<p class='note'>The likelihood difference concerns the group-unit stage. The available Wald tests assess the group- and set-unit families.</p>",
        .html_table(data.frame(
          model = c("Equal group units", "Group-dependent units"),
          loglik = c(x$ll_equal, x$ll_efrm),
          unit_parameters = c(0L, x$extra_parameters),
          two_delta_loglik = c(NA_real_, x$two_delta_ll))),
        "<h2>Unit tests</h2>",
        "<h3>Omnibus tests</h3>", .html_table(x$unit_omnibus),
        "<h3>Unit contrasts</h3>", .html_table(x$unit_tests),
        "<h2>Frames and units</h2>",
        "<h3>Frames</h3>",
        .html_table(fit$frames[, intersect(c("set", "group", "rho",
          "se_log_rho", "origin", "fit_resid", "n_responses"),
          names(fit$frames)), drop = FALSE]),
        "<h3>Group units</h3>", .html_table(fit$phi_table),
        "<h3>Item-set units</h3>", .html_table(fit$alpha_table),
        "<h3>Item-set locations</h3>", .html_table(fit$set_table))
    } else "",
    if (inherits(fit, "rasch_efrm"))
      .html_frame_invariance(fit, invariance) else "",
    "<h2>Person estimates</h2>",
    .html_table(fit$person[, intersect(c("id", names(fit$factors), "raw",
                                         "max_raw", "theta", "se", "extreme",
                                         "fit_resid"),
                                       names(fit$person))]),
    if (!is.null(person_weight_result)) s(
      "<h2>Externally weighted secondary person measures</h2>",
      "<p class='note'>These measures use the supplied relative item or item-set weights. They do not replace the ordinary person estimates used elsewhere in the analysis.</p>",
      .html_table(person_weight_result$table),
      "<h3>Resolved external weights</h3>",
      .html_table(attr(person_weight_result$table, "weighting", exact = TRUE))) else "",
    "</div></body></html>")
  writeLines(html, file, useBytes = TRUE)
  invisible(file)
}

#' Write an editable or print-ready analysis report
#'
#' Renders the active Rasch or paired-comparison fit as a self-contained HTML
#' document, an editable Word document, or a PDF. The report contains the
#' principal estimates, model-specific tables, diagnostic figures, and
#' software provenance. Complete machine-readable results remain available
#' from \code{\link{save_outputs}}. Reports downloaded from the application
#' retain compatible tailored item shifts and externally weighted secondary
#' person measures.
#' For keyed fits with repeated person IDs, the report explains why distractor
#' analysis is unavailable and retains the other model outputs.
#'
#' @param fit A fitted object from \code{\link{rasch}}, \code{\link{rasch_mfrm}},
#'   \code{\link{rasch_efrm}}, \code{\link{btl}}, or \code{\link{btl_efrm}}.
#' @param file Output path ending in \code{.html}, \code{.docx}, or \code{.pdf}.
#' @param format Output format. By default it is inferred from \code{file}.
#' @param dif,bootstrap Optional computed \code{\link{dif_anova}} and
#'   \code{\link{fit_bootstrap}} results from this fit, rendered as run rather than
#'   recomputed at defaults.
#' @param dif_bootstrap Optional \code{\link{dif_bootstrap}} sensitivity
#'   analysis from this fit and DIF specification.
#' @param dimensionality Optional computed \code{\link{plot_scree}} or
#'   \code{\link{btl_dimensionality}} result from this fit.
#' @param subtest Optional computed \code{\link{dimensionality_test}} result
#'   from this fit.
#' @param invariance Optional computed \code{\link{frame_invariance}} result
#'   from an EFRM fit.
#' @param tailored Optional computed \code{\link{tailored_analysis}} result
#'   from this ordinary dichotomous fit.
#' @param person_weights Optional table returned by
#'   \code{\link{weighted_person_estimates}}. An explicit table takes
#'   precedence over a compatible result retained by the application.
#' @param title Report title.
#' @return Invisibly, the output path.
#' @details Word and HTML output require Pandoc, supplied with RStudio and
#'   available through \pkg{rmarkdown}. PDF output also requires a LaTeX
#'   installation such as TinyTeX.
#' @examples
#' \dontrun{
#' fit <- rasch(matrix(rbinom(3000, 1, .5), 300, 10))
#' report_document(fit, file.path(tempdir(), "analysis.docx"))
#' }
#' @export
report_document <- function(fit, file,
                            format = c("auto", "html", "docx", "pdf"),
                            title = "Rasch measurement analysis",
                            dif = NULL, bootstrap = NULL,
                            dif_bootstrap = NULL, dimensionality = NULL,
                            invariance = NULL, subtest = NULL,
                            tailored = NULL, person_weights = NULL) {
  .check_export_fit(fit)
  person_weight_result <- .resolve_report_person_weights(fit, person_weights)
  .validate_boot_dif_result(dif, fit)
  .validate_fit_bootstrap(bootstrap, fit)
  if (!is.null(dif_bootstrap) && is.null(dif))
    stop("`dif` must accompany `dif_bootstrap` so the primary and sensitivity analyses use the same specification")
  .validate_dif_bootstrap(dif_bootstrap, fit, dif)
  if (inherits(fit, "rasch_btl")) {
    if (!is.null(subtest))
      stop("`subtest` is not available for a paired-comparison fit")
    .validate_btl_dimensionality(dimensionality, fit)
  } else {
    .validate_scree_result(dimensionality, fit)
    .validate_dimensionality_test(subtest, fit)
  }
  .validate_frame_invariance(invariance, fit)
  .validate_tailored_result(tailored, fit)
  # computed results travel as attributes on the serialised fit, so the
  # template renders the analysis as run rather than a default recomputation
  if (!is.null(dif)) attr(fit, "report_dif") <- dif
  if (!is.null(bootstrap)) attr(fit, "report_bootstrap") <- bootstrap
  if (!is.null(dif_bootstrap))
    attr(fit, "report_dif_bootstrap") <- dif_bootstrap
  if (!is.null(dimensionality))
    attr(fit, "report_dimensionality") <- dimensionality
  if (!is.null(subtest)) attr(fit, "report_subtest") <- subtest
  if (!is.null(invariance)) attr(fit, "report_invariance") <- invariance
  if (!is.null(tailored)) attr(fit, "report_tailored") <- tailored
  if (!is.null(person_weights))
    attr(fit, "report_person_weights") <- person_weight_result
  .check_out_path(file, "file")
  if (length(title) != 1L || !is.character(title) || !is.null(dim(title)) ||
      !is.null(oldClass(title)) || is.na(title))
    stop("`title` must be one non-missing title")
  format <- match.arg(format)
  ext <- tolower(tools::file_ext(file))
  if (format == "auto") {
    format <- switch(ext, html = "html", htm = "html",
                     docx = "docx", pdf = "pdf", NA_character_)
    if (is.na(format))
      stop("infer the report format from a .html, .docx, or .pdf filename")
  }
  wanted <- c(html = "html", docx = "docx", pdf = "pdf")[[format]]
  if (!ext %in% c(wanted, if (format == "html") "htm"))
    stop("the filename extension does not match the requested format")
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("report_document() needs the suggested package rmarkdown")
  if (!rmarkdown::pandoc_available())
    stop("Pandoc is unavailable; install it or use RStudio's bundled Pandoc")

  template <- system.file("rmarkdown", "rasch-report.Rmd", package = "rasch")
  if (!nzchar(template)) {
    candidate <- file.path("inst", "rmarkdown", "rasch-report.Rmd")
    if (file.exists(candidate)) template <- candidate
  }
  if (!nzchar(template) || !file.exists(template))
    stop("the analysis report template is missing")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  fit_file <- tempfile("rasch-report-fit-", fileext = ".rds")
  saveRDS(fit, fit_file, version = 3)
  on.exit(unlink(fit_file), add = TRUE)
  meta <- c("--metadata", paste0("title=", title))
  output_format <- switch(
    format,
    html = rmarkdown::html_document(self_contained = TRUE,
                                    pandoc_args = meta),
    docx = rmarkdown::word_document(pandoc_args = meta),
    pdf = rmarkdown::pdf_document(latex_engine = "xelatex",
                                  pandoc_args = meta))
  env <- new.env(parent = asNamespace("rasch"))
  rendered <- rmarkdown::render(
    input = template, output_format = output_format,
    output_file = basename(file), output_dir = dirname(file),
    params = list(title = title, fit_file = fit_file),
    envir = env, quiet = TRUE)
  if (!file.exists(rendered))
    stop("the report renderer did not create the requested file")
  invisible(normalizePath(rendered, mustWork = TRUE))
}

# The frame model holds each item at one location across frames and scales it
# by the frame unit, so it cannot test that assumption from its own fit. A
# report that shows the units without the test invites the reader to trust
# them further than the analysis warrants, so the test travels with them.
.html_frame_invariance <- function(fit, invariance = NULL) {
  inv <- if (!is.null(invariance)) invariance else
    tryCatch(frame_invariance(fit), error = function(e) e)
  if (inherits(inv, "error"))
    return(paste0("<h2>Item invariance across frames</h2>",
      "<p class='note'>Frame invariance was not available: ",
      .html_escape(conditionMessage(inv)), "</p>"))
  bootstrap <- identical(inv$se_method, "bootstrap")
  fl <- inv$locations[inv$locations$flagged %in% TRUE, , drop = FALSE]
  fd <- inv$discrimination[inv$discrimination$flagged %in% TRUE, ,
                           drop = FALSE]
  nb <- sum(inv$discrimination$disc_boundary %in% TRUE)
  paste0("<h2>Item invariance across frames</h2>",
    if (bootstrap) paste0(
      "<p class='note'>Each frame is calibrated separately and compared on the",
      " common scale. Whole-person bootstrap uncertainty within person group includes the",
      " fitted frame units; Holm adjustment covers ", inv$family_n,
      " location and discrimination comparisons. ", inv$boot_reps_used, " of ",
      inv$boot_reps, " refits were usable; ",
      inv$boot_reps_nonconverged, " did not converge and ",
      inv$boot_reps_errors, " otherwise failed.</p>") else paste0(
      "<p class='note'>Each frame is calibrated separately and compared on the",
      " common scale. Conditional location tests treat the fitted frame units",
      " as fixed and are Holm-adjusted over ", inv$family_n,
      " location comparisons. Discrimination",
      " comparisons are descriptive; bootstrap uncertainty is required for",
      " discrimination tests.</p>"),
    .html_table(as.data.frame(inv$summary)),
    if (!is.null(inv$excluded) && nrow(inv$excluded))
      paste0("<h3>Excluded comparisons</h3>",
             .html_table(as.data.frame(inv$excluded))) else "",
    if (nrow(fl)) paste0("<h3>Locations differing across frames</h3>",
      .html_table(as.data.frame(fl[, intersect(
        c("set", "frame_1", "frame_2", "item", "location_1", "location_2",
          "difference", "se", "statistic", "p_adj"), names(fl))])))
    else "<p class='note'>No available item-location comparison differs across frames.</p>",
    if (!bootstrap) paste0("<h3>Descriptive discrimination comparisons</h3>",
      .html_table(as.data.frame(inv$discrimination[, intersect(
        c("set", "frame_1", "frame_2", "item", "infit_1", "infit_2",
          "infit_z", "disc_1", "disc_2", "disc_ratio", "disc_boundary"),
        names(inv$discrimination)), drop = FALSE])))
    else paste0(
      if (nb) sprintf(paste0(
        "<p class='note'>%d discrimination probability/probabilities were ",
        "withheld because a fitted slope is on its boundary.</p>"), nb) else "",
      if (nrow(fd)) paste0(
      "<h3>Discrimination differing across frames</h3>",
      .html_table(as.data.frame(fd[, intersect(
        c("set", "frame_1", "frame_2", "item", "log_disc_ratio",
          "se_log_disc_ratio", "statistic", "p_adj", "disc_1", "disc_2",
          "disc_ratio", "disc_boundary"), names(fd)), drop = FALSE])))
      else "<p class='note'>No available discrimination comparison differs across frames.</p>"))
}
