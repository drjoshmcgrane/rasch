# rasch :: summary tables
# ===========================================================================
# The test-of-fit and targeting/reliability summaries as tidy two-column
# tables, so the headline statistics of an analysis can be saved and
# reported rather than read off a text panel.
# ===========================================================================

.inference_count <- function(p, alpha = 0.05) {
  total <- length(p)
  tested <- sum(is.finite(p))
  flagged <- if (tested) sum(p[is.finite(p)] < alpha) else NA_integer_
  unavailable <- total - tested
  text <- if (!tested) "unavailable"
    else if (!unavailable) sprintf("%d of %d", flagged, tested)
    else sprintf("%d of %d tested (%d unavailable)",
                 flagged, tested, unavailable)
  list(flagged = flagged, tested = tested, unavailable = unavailable,
       total = total, text = text)
}

# The estimator is a property of the fit, not a constant: a structural or
# explanatory model does not use the ordinary pairwise conditional routine,
# and naming the wrong one in an exported table contradicts the panel the
# download sits beside.
.estimation_label <- function(fit) {
  if (inherits(fit, "rasch") &&
      (isTRUE(fit$refit_spec$fixed_calibration) ||
       isTRUE(fit$est$n_parameters == 0L)))
    return("fixed item calibration; person scoring")
  if (!is.null(fit$estimation) && nzchar(fit$estimation)) return(fit$estimation)
  if (inherits(fit, "rasch_efrm") && length(unique(fit$set_of)) > 1L)
    "pairwise conditional calibration + semiparametric set linking"
  else if (inherits(fit, c("rasch_efrm", "rasch_mfrm")))
    "pairwise conditional ML over response cells"
  else if (inherits(fit, "rasch_explanatory"))
    "pairwise conditional ML over the explanatory design"
  else "pairwise conditional ML"
}

#' Test-of-fit summary as a table
#'
#' Returns the model, estimation method, trait chi-square, calibration and
#' person fit-residual moments, fit-location correlations, the approximate
#' asymptotic chi-square flag count, and disordered-threshold count as a
#' two-column table. MFRM and EFRM
#' summaries label their fitted item-by-facet or item-by-frame columns as
#' response cells. A wholly unavailable probability family is labelled as
#' unavailable; a partial family reports the tested and unavailable counts.
#'
#' @param fit A fitted object from \code{\link{rasch}}, or a
#'   paired-comparison fit from \code{\link{btl}} (which reports its own
#'   headline set: convergence, pairwise chi-square, object separation,
#'   thresholds structure, and dependence effects when estimated).
#' @return A data frame with columns \code{statistic} and \code{value}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit_summary_table(rasch(X))
#' @export
fit_summary_table <- function(fit) {
  if (inherits(fit, "rasch_btl")) return(.btl_summary_table(fit))
  if (!inherits(fit, "rasch"))
    stop("`fit` must be a fitted model from rasch(), rasch_mfrm(), ",
         "rasch_efrm(), or btl()", call. = FALSE)
  ss <- fit$summary_stats
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  unit <- if (structural) "Response cell" else "Item"
  units <- if (structural) "Response cells" else "Items"
  trait_unit <- if (structural) "response-cell" else "item"
  dis <- names(which(vapply(fit$thresholds_diag, function(d)
    !d$ordered && length(d$thresholds) > 1, TRUE)))
  repeated_ids <- .has_repeated_residual_units(fit)
  item_p <- if (repeated_ids) rep(NA_real_, nrow(fit$items))
    else fit$items$p_adj
  item_inference <- .inference_count(item_p)
  total_p <- if (repeated_ids) NA_real_ else fit$total_chisq_p
  converged <- isTRUE(fit$est$converged) && .efrm_link_converged(fit)
  num <- function(x, d = 3) formatC(x, digits = d, format = "f")
  out <- data.frame(statistic = c(
    "Model", "Estimation", "Converged", "Iterations",
    paste("Approximate asymptotic total",
          paste0(trait_unit, "-trait chi-square")),
    "Degrees of freedom",
    paste0("Approximate asymptotic ",
           if (structural) "response-cell" else "item",
           "-trait probability"), "Class intervals",
    paste(unit, "fit residual mean"), paste(unit, "fit residual SD"),
    paste(unit, "fit residual skewness"),
    paste(unit, "fit residual kurtosis"),
    "Person fit residual mean", "Person fit residual SD",
    "Person fit residual skewness", "Person fit residual kurtosis",
    paste0("Fit-location correlation (", tolower(units), ")"),
    "Fit-location correlation (persons)",
    paste(units, "with approximate asymptotic Holm p < .05"),
    if (structural) "Disordered response-cell thresholds" else
      "Disordered thresholds"),
    value = c(
      if (inherits(fit, "rasch_explanatory")) fit$explanatory_model else fit$model,
      .estimation_label(fit),
      ifelse(converged, "yes", "NO"),
      as.character(fit$est$iterations),
      num(fit$total_chisq), as.character(fit$total_df),
      if (is.finite(total_p)) .fmt_p(total_p) else "unavailable",
      as.character(fit$n_groups),
      num(fit$item_fit_summary$mean, 2), num(fit$item_fit_summary$sd, 2),
      num(fit$item_fit_summary$skewness, 2),
      num(fit$item_fit_summary$kurtosis, 2),
      num(fit$person_fit_summary$mean, 2), num(fit$person_fit_summary$sd, 2),
      num(fit$person_fit_summary$skewness, 2),
      num(fit$person_fit_summary$kurtosis, 2),
      num(ss$cor_item_fit_location), num(ss$cor_person_fit_location),
      item_inference$text,
      if (length(dis)) paste(dis, collapse = ", ") else "none"))
  rownames(out) <- NULL
  out
}

#' Targeting and reliability summary as a table
#'
#' The person and calibration location moments, threshold range and coverage,
#' and the applicable separation and reliability indices as a two-column
#' table suitable for saving and reporting. MFRM and EFRM fits describe
#' item-by-facet or item-by-frame response cells rather than additional items.
#' Coefficient alpha is not applicable when an item has several response
#' cells; it is retained for the one-cell-per-item reduction. Item separation
#' excludes wholly fixed item locations in an anchored calibration.
#' A targeting summary requires a converged calibration and, for EFRM, a
#' converged set-unit link.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return A data frame with columns \code{statistic} and \code{value}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
#' colnames(X) <- paste0("I", 1:6)
#' targeting_table(rasch(X))
#' @export
targeting_table <- function(fit) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be a response-data fit from rasch(), rasch_mfrm(), ",
         "or rasch_efrm()", call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; the targeting summary is unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; the targeting summary is unavailable",
         call. = FALSE)
  ss <- fit$summary_stats; t <- fit$targeting
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  alpha_design <- .classical_design_applicable(fit)
  cal <- if (structural) "Response-cell" else "Item"
  threshold <- if (structural) "Calibration threshold" else "Threshold"
  num <- function(x, d = 3) ifelse(is.finite(x),
                                   formatC(x, digits = d, format = "f"), "NA")
  out <- data.frame(statistic = c(
    "Person location mean", "Person location SD",
    "Person location skewness", "Person location kurtosis",
    "Person location mean (no extremes)", "Person location SD (no extremes)",
    paste(cal, "location mean (constrained)"), paste(cal, "location SD"),
    paste(cal, "location skewness"), paste(cal, "location kurtosis"),
    paste(threshold, "minimum"), paste(threshold, "maximum"),
    paste("Persons below", tolower(threshold), "range (%)"),
    paste("Persons above", tolower(threshold), "range (%)"),
    "PSI", "Separation", "Person strata", "PSI without extremes",
    "n without extremes",
    if (structural) "Response-cell separation reliability" else
      "Item separation reliability",
    "Coefficient alpha",
    "n complete (alpha)"),
    value = c(
      num(ss$person_location$mean), num(ss$person_location$sd),
      num(ss$person_location$skewness, 2), num(ss$person_location$kurtosis, 2),
      num(ss$person_location_noext$mean), num(ss$person_location_noext$sd),
      num(ss$item_location$mean), num(ss$item_location$sd),
      num(ss$item_location$skewness, 2), num(ss$item_location$kurtosis, 2),
      num(t$threshold_range[1]), num(t$threshold_range[2]),
      num(100 * t$prop_below, 1), num(100 * t$prop_above, 1),
      num(fit$psi$PSI), num(fit$psi$separation, 2),
      num(fit$psi$strata, 1),
      num(fit$psi_noext$PSI), as.character(fit$psi_noext$n),
      num(fit$isi$PSI),
      if (!alpha_design) "not applicable" else num(fit$alpha$alpha),
      if (!alpha_design) "not applicable" else if (is.null(fit$alpha$n)) "NA"
        else as.character(fit$alpha$n)))
  rownames(out) <- NULL
  out
}

# The paired-comparison headline set, same two-column shape.
.btl_summary_table <- function(fit) {
  num <- function(x, d = 3) ifelse(is.finite(x),
                                   formatC(x, digits = d, format = "f"), "NA")
  framed <- inherits(fit, "rasch_btl_efrm")
  rows <- list(
    c("Model", if (framed) "Paired comparisons with frame-dependent units"
      else if (inherits(fit, "rasch_btl_explanatory"))
        "Explanatory comparative judgement"
      else if (!is.null(fit$m) && fit$m > 1L)
      sprintf("Graded paired comparisons (%d categories)", fit$m + 1L)
      else "Paired comparisons (BTL)"),
    c("Estimation", if (framed) "two-stage maximum likelihood"
      else "maximum likelihood"),
    c("Converged", ifelse(isTRUE(fit$converged), "yes", "NO")),
    c("Objects", as.character(nrow(fit$objects))),
    c("Comparisons", formatC(fit$n_comparisons, format = "d")),
    c("Judges", if (!is.null(fit$judges)) as.character(nrow(fit$judges))
      else "not identified"),
    c("Standard errors", if (framed)
      switch(fit$se_method,
             judge_bootstrap = "judge-resampling bootstrap",
             bootstrap = "parametric bootstrap",
             conditional = "conditional analytic")
      else if (isTRUE(fit$clustered)) "sandwich, clustered by judge"
      else "sandwich"),
    c("Pairwise chi-square", num(fit$total_chisq, 2)),
    c("Degrees of freedom", as.character(fit$total_df)),
    c("Pairwise fit probability", if (is.finite(fit$total_p))
      .fmt_p(fit$total_p) else if (framed || isTRUE(fit$clustered))
        "unavailable (judge clustering)" else "unavailable"),
    c("Object separation index", num(fit$osi$PSI)),
    c("Separation", num(fit$osi$separation)),
    c("Log-likelihood", num(fit$loglik, 2)))
  if (!framed)
    rows <- append(rows, list(c("Iterations", as.character(fit$iterations))),
                   after = 3L)
  else {
    rows[[length(rows) + 1L]] <- c("Object sets", length(fit$sets))
    rows[[length(rows) + 1L]] <- c("Judge panels", length(fit$panels))
  }
  if (!is.null(fit$thr_structure) && !is.null(fit$m) && fit$m > 1L)
    rows[[length(rows) + 1L]] <- c("Threshold structure",
                                   if (fit$thr_structure == "pc")
                                     "principal components (spread)"
                                   else "free symmetric")
  if (!is.null(fit$dependence)) {
    for (r in seq_len(nrow(fit$dependence))) {
      shown_p <- if (is.numeric(fit$dependence$p_adj) &&
                     length(fit$dependence$p_adj) >= r &&
                     is.finite(fit$dependence$p_adj[r]))
        paste0("Holm p = ", .fmt_p(fit$dependence$p_adj[r])) else
          "Holm p unavailable"
      effect_label <- if (identical(fit$dependence$effect[r], "position"))
        "First-position advantage" else
          sprintf("Within-judge %s",
                  gsub("_", "-", fit$dependence$effect[r]))
      rows[[length(rows) + 1L]] <- c(
        sprintf("%s (logits)", effect_label),
        sprintf("%s (SE %s, %s)", num(fit$dependence$estimate[r]),
                num(fit$dependence$se[r]), shown_p))
    }
  }
  out <- data.frame(statistic = vapply(rows, `[`, "", 1),
                    value = vapply(rows, `[`, "", 2))
  rownames(out) <- NULL
  out
}
