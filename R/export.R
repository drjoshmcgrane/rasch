# RaschR :: export
# ===========================================================================
# save_outputs() writes the complete analysis to disk: every table as CSV,
# every plot as PNG (and optionally PDF), and a plain-text summary. The
# Shiny app zips the resulting folder for its "download everything" button.
# ===========================================================================

.rr_device <- function(path, fmt, width, height, dpi) {
  if (fmt == "png") png(path, width = width, height = height, units = "in", res = dpi)
  else pdf(path, width = width, height = height)
}

.rr_save_plot <- function(expr, stem, dir, formats, width, height, dpi) {
  files <- character(0)
  for (fmt in formats) {
    path <- file.path(dir, paste0(stem, ".", fmt))
    ok <- tryCatch({
      .rr_device(path, fmt, width, height, dpi)
      force(expr())
      dev.off()
      TRUE
    }, error = function(e) { try(dev.off(), silent = TRUE); FALSE })
    if (ok) files <- c(files, path)
  }
  files
}

#' Save every output of a Rasch analysis to a folder
#'
#' Writes all tables (item statistics, thresholds with standard errors, person
#' estimates including ID and factors, the score-to-measure table, residual
#' correlations, principal-component loadings, category frequencies, and DIF
#' results for every nominated factor) as CSV; every plot, including the
#' per-item characteristic, category, threshold, and frequency plots, as PNG
#' and optionally PDF; and a plain-text analysis summary.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param dir Output directory; created if absent.
#' @param formats Plot formats, any of \code{"png"} and \code{"pdf"}.
#' @param width,height Plot size in inches.
#' @param dpi PNG resolution; the default 300 is publication quality.
#' @param item_plots Also write the per-item plot set (one ICC, category curve,
#'   threshold curve, and frequency chart per item).
#' @return Invisibly, the vector of files written.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' out <- file.path(tempdir(), "rasch-out")
#' save_outputs(rasch(X), out, formats = "png", item_plots = FALSE)
#' @export
save_outputs <- function(fit, dir, formats = c("png", "pdf"), width = 9,
                         height = 6, dpi = 300, item_plots = TRUE) {
  formats <- match.arg(formats, c("png", "pdf"), several.ok = TRUE)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tdir <- file.path(dir, "tables"); pdir <- file.path(dir, "plots")
  idir <- file.path(pdir, "items")
  dir.create(tdir, showWarnings = FALSE)
  dir.create(pdir, showWarnings = FALSE)
  if (item_plots) dir.create(idir, showWarnings = FALSE)
  files <- character(0)
  wtab <- function(d, name) {
    path <- file.path(tdir, paste0(name, ".csv"))
    utils::write.csv(d, path, row.names = FALSE)
    files <<- c(files, path)
  }

  # --- tables ---------------------------------------------------------------
  wtab(fit$items, "item_statistics")
  thr <- fit$thresholds
  thr$item <- fit$items$item[thr$item]
  wtab(thr[, c("item", "k", "tau", "se")], "thresholds")
  wtab(fit$person, "person_estimates")
  wtab(fit$score_table, "score_to_measure")
  rc <- residual_correlations(fit)
  wtab(data.frame(item = rownames(rc$matrix), round(rc$matrix, 4),
                  check.names = FALSE), "residual_correlations")
  if (nrow(rc$flagged)) wtab(rc$flagged, "local_dependence_flagged")
  pc <- residual_pca(fit)
  wtab(pc$loadings_matrix, "pca_loadings")
  wtab(pc$eigen_table, "residual_eigenvalues")
  cf <- do.call(rbind, lapply(fit$thresholds_diag, function(d)
    data.frame(item = d$item, category = seq_along(d$category_counts) - 1L,
               count = d$category_counts)))
  wtab(cf, "category_frequencies")
  gt <- guttman_table(fit)
  wtab(data.frame(id = rownames(gt$matrix), gt$matrix, check.names = FALSE),
       "guttman_ordered_responses")
  if (!is.null(fit$mc)) wtab(distractor_analysis(fit), "distractor_analysis")
  if (inherits(fit, "rasch_mfrm")) {
    wtab(fit$item_effects, "item_effects")
    wtab(fit$item_thresholds, "item_structural_thresholds")
    for (f in fit$facet_spec)
      wtab(fit$facet_effects[[f]], paste0("facet_", gsub("[^A-Za-z0-9_.-]", "_", f)))
    if (!is.null(fit$interaction_effects))
      wtab(fit$interaction_effects, "item_by_facet_interactions")
  }
  if (inherits(fit, "rasch_efrm")) {
    wtab(fit$frames, "frames")
    wtab(fit$phi_table, "group_units_phi")
    wtab(fit$alpha_table, "set_units_alpha")
    wtab(fit$set_table, "set_locations")
    wtab(fit$item_arbitrary, "items_common_unit")
    wtab(fit$thresholds_arbitrary, "thresholds_common_unit")
    wtab(fit$score_curves, "score_curves")
  }
  if (!is.null(fit$factors)) {
    dif <- tryCatch(dif_anova(fit), error = function(e) NULL)
    if (!is.null(dif)) wtab(dif, "dif_anova")
  }

  # --- summary ---------------------------------------------------------------
  spath <- file.path(dir, "summary.txt")
  con <- file(spath, "w")
  sink(con); on.exit({ sink(); close(con) }, add = TRUE)
  summary(fit)
  dt <- dimensionality_test(fit)
  if (is.null(dt$note)) {
    cat(sprintf("\nUnidimensionality t-test: %.1f%% significant (exact 95%% CI %.1f%% to %.1f%%), %s\n",
                100 * dt$prop_significant, 100 * dt$ci[1], 100 * dt$ci[2],
                if (dt$multidimensional) "MULTIDIMENSIONAL" else "consistent with one dimension"))
  }
  cat(sprintf("Average residual correlation: %.3f; %d flagged dependent pair(s)\n",
              rc$average, nrow(rc$flagged)))
  sink(); close(con); on.exit()
  files <- c(files, spath)

  # --- test-level plots --------------------------------------------------------
  sp <- function(f, stem) files <<- c(files,
    .rr_save_plot(f, stem, pdir, formats, width, height, dpi))
  sp(function() plot_pimap(fit), "person_item_distribution")
  sp(function() plot_threshold_map(fit), "threshold_map")
  sp(function() plot_tcc(fit), "test_characteristic_curve")
  sp(function() plot_tif(fit), "test_information")
  sp(function() plot_item_map(fit), "item_fit_map")
  sp(function() plot_person_fit(fit), "person_fit")
  sp(function() plot_resid_cor(fit), "residual_correlations")
  sp(function() plot_pca(fit), "pca_loadings")
  sp(function() plot_scree(fit), "scree")
  sp(function() plot_guttman(fit), "guttman_scalogram")
  if (inherits(fit, "rasch_mfrm")) {
    for (f in fit$facet_spec) local({
      f_ <- f
      sp(function() plot_facets(fit, f_),
         paste0("facet_severities_", gsub("[^A-Za-z0-9_.-]", "_", f_)))
    })
  }
  if (inherits(fit, "rasch_efrm")) {
    sp(function() plot_frames(fit), "frame_units")
    if (item_plots) for (it in unique(fit$virtual_map$item)) local({
      it_ <- it
      files <<- c(files, .rr_save_plot(function() plot_icc_frames(fit, it_),
        paste0(gsub("[^A-Za-z0-9_.-]", "_", it_), "_icc_frames"),
        idir, formats, width, height, dpi))
    })
  }

  # --- per-item plots ------------------------------------------------------------
  if (item_plots && !is.null(fit$mc)) {
    for (it in colnames(fit$mc$raw)) local({
      it_ <- it
      files <<- c(files, .rr_save_plot(function() plot_distractors(fit, it_),
        paste0(gsub("[^A-Za-z0-9_.-]", "_", it_), "_options"),
        idir, formats, width, height, dpi))
    })
  }
  if (item_plots) {
    for (it in fit$items$item) {
      safe <- gsub("[^A-Za-z0-9_.-]", "_", it)
      files <- c(files,
        .rr_save_plot(function() plot_icc(fit, it),
                      paste0(safe, "_icc"), idir, formats, width, height, dpi),
        .rr_save_plot(function() plot_ccc(fit, it),
                      paste0(safe, "_categories"), idir, formats, width, height, dpi),
        .rr_save_plot(function() plot_threshold_prob(fit, it),
                      paste0(safe, "_thresholds"), idir, formats, width, height, dpi),
        .rr_save_plot(function() plot_catfreq(fit, it),
                      paste0(safe, "_frequencies"), idir, formats, width, height, dpi))
    }
  }
  invisible(files)
}
