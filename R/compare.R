# RaschR :: model comparison
# ===========================================================================
# Side-by-side comparison of fitted models. Two kinds of evidence are
# reported. (1) The pairwise conditional log-likelihood with the number of
# structural parameters and, for fits of the SAME response data (identical
# items, categories, and persons, hence identical conditional information),
# twice the log-likelihood difference from the reference fit. Because the
# likelihood is a composite (pairwise) one, the difference is descriptive
# and is not chi-square calibrated; it is most meaningful for nested
# structures (for example the rating scale model inside the partial credit
# model, or equal units inside the extended frame of reference model).
# (2) Calibration-free fit descriptors that remain comparable across
# different data preparations: the total item-trait chi-square against its
# degrees of freedom, the spread of the item and person fit residuals
# (ideal SD 1), the person separation index, and Cronbach's alpha.
# ===========================================================================

#' Compare fitted Rasch models
#'
#' Builds a comparison table for two or more fits from \code{\link{rasch}},
#' \code{\link{rasch_mfrm}}, or \code{\link{rasch_efrm}}. For fits of the
#' same response data (identical item columns, maximum scores, and number of
#' persons) the pairwise conditional log-likelihoods share their conditional
#' information, and twice the difference from the reference fit is reported
#' with the difference in parameter counts; this is descriptive (composite
#' likelihood), and most meaningful for nested structures such as RSM inside
#' PCM. Across different data preparations (subtests, splits, facet or frame
#' structures) the likelihoods are not comparable and the calibration-free
#' columns carry the comparison: total item-trait chi-square per degree of
#' freedom, item and person fit residual SDs (ideal 1), PSI, and alpha.
#'
#' @param ... Two or more fitted objects, ideally named
#'   (\code{compare_fits(PCM = f1, RSM = f2)}).
#' @param reference Index or name of the reference fit for the
#'   log-likelihood difference; defaults to the first.
#' @return A data frame with one row per fit: label, model, persons, items,
#'   parameters, log-likelihood, comparability with the reference,
#'   \code{two_delta_ll} and \code{delta_parameters} (same-data fits only),
#'   chi-square per df, fit residual SDs, PSI, and alpha.
#' @examples
#' set.seed(1)
#' simP <- function(th, tau) { x <- 0:length(tau); p <- exp(x * th - c(0, cumsum(tau))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(seq(-1, 1, length.out = 6), function(b)
#'   sapply(th, function(t) sample(0:3, 1, prob = simP(t, b + c(-0.8, 0, 0.8)))))
#' colnames(X) <- paste0("R", 1:6)
#' compare_fits(PCM = rasch(X, model = "PCM"), RSM = rasch(X, model = "RSM"))
#' @export
compare_fits <- function(..., reference = 1) {
  fits <- list(...)
  if (length(fits) < 2) stop("supply at least two fits to compare")
  bad <- !vapply(fits, inherits, TRUE, what = "rasch")
  if (any(bad)) stop("argument(s) ", paste(which(bad), collapse = ", "),
                     " are not rasch fits")
  labs <- names(fits)
  if (is.null(labs)) labs <- rep("", length(fits))
  labs[labs == ""] <- paste0("fit", seq_along(fits))[labs == ""]
  if (is.character(reference)) reference <- match(reference, labs)
  if (is.na(reference) || reference < 1 || reference > length(fits))
    stop("no such reference fit")

  sig <- function(f) list(items = colnames(f$X), m = unname(f$m),
                          n = nrow(f$X))
  ref_sig <- sig(fits[[reference]])
  rows <- lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    s <- sig(f)
    same <- identical(s, ref_sig)
    data.frame(
      label = labs[i], model = f$model,
      persons = nrow(f$X), items = ncol(f$X),
      parameters = if (is.null(f$est$n_parameters)) NA_integer_
                   else f$est$n_parameters,
      loglik = f$est$loglik,
      same_data = same,
      two_delta_ll = NA_real_, delta_parameters = NA_integer_,
      chisq_per_df = f$total_chisq / f$total_df,
      item_fit_sd = f$item_fit_summary$sd,
      person_fit_sd = f$person_fit_summary$sd,
      PSI = f$psi$PSI, alpha = f$alpha$alpha)
  })
  out <- do.call(rbind, rows)
  ref <- out[reference, ]
  cmp <- out$same_data & seq_len(nrow(out)) != reference
  out$two_delta_ll[cmp] <- 2 * (out$loglik[cmp] - ref$loglik)
  out$delta_parameters[cmp] <- out$parameters[cmp] - ref$parameters
  rownames(out) <- NULL
  attr(out, "reference") <- labs[reference]
  attr(out, "note") <- paste(
    "two_delta_ll is a composite (pairwise) likelihood difference against",
    "the reference fit, reported only for fits of the same response data;",
    "it is descriptive, not chi-square calibrated. Across different data",
    "preparations compare chisq_per_df, the fit residual SDs (ideal 1),",
    "PSI, and alpha.")
  class(out) <- c("rasch_compare", "data.frame")
  out
}

#' @export
print.rasch_compare <- function(x, ...) {
  cat("Model comparison (reference:", attr(x, "reference"), ")\n\n")
  y <- as.data.frame(x)
  num <- vapply(y, is.numeric, TRUE)
  y[num] <- lapply(y[num], round, 3)
  print(y, row.names = FALSE)
  cat("\n", attr(x, "note"), "\n", sep = "")
  invisible(x)
}
