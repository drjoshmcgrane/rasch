# RaschR :: differential item functioning
# ===========================================================================
# DIF by analysis of variance of the standardised residuals (Hagquist &
# Andrich 2017). For each item, residuals are analysed by person factor and
# trait class interval: a factor main effect indicates uniform DIF and a
# factor-by-interval interaction indicates non-uniform DIF. With several
# person factors the analysis can be run factor-at-a-time (dif_anova) or as
# one factorial model per item (dif_anova_factorial) with factor-by-factor
# interactions, Tukey HSD comparisons on the significant group terms, and
# the convention that a significant interaction supersedes the main effects
# of the factors involved. Multiplicity across items is handled by
# Benjamini-Hochberg false-discovery-rate adjustment.
# ===========================================================================

.dif_factors <- function(fit, factors) {
  if (is.null(factors)) factors <- fit$factors
  if (is.null(factors)) stop("no person factors supplied or stored in the fit")
  if (is.character(factors) && !is.null(fit$factors) &&
      all(factors %in% names(fit$factors)))
    factors <- fit$factors[, factors, drop = FALSE]
  if (!is.data.frame(factors)) factors <- data.frame(group = factors)
  factors
}

.dif_class_intervals <- function(fit, n_groups) {
  ci <- fit$person$class_interval
  if (is.null(ci) || !identical(n_groups, fit$n_groups))
    ci <- .class_intervals(fit$person$theta, fit$person$extreme, n_groups)
  factor(ci)
}

#' Differential item functioning by two-way residual ANOVA
#'
#' For each item and each person factor separately, analyses the
#' standardised residuals by factor group and trait class interval. The
#' group main effect indicates uniform DIF and the group-by-interval
#' interaction indicates non-uniform DIF. Probabilities are adjusted across
#' items within each factor by the Benjamini-Hochberg false-discovery-rate
#' procedure (or any \code{\link[stats]{p.adjust}} method). With several
#' factors, consider \code{\link{dif_anova_factorial}}, which models them
#' jointly.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors A vector (one factor), a data frame of person factors, or a
#'   character vector naming factor columns already nominated in the fit (via
#'   \code{rasch(..., factors = )}). Defaults to every factor stored in the
#'   fit.
#' @param n_groups Number of trait class intervals; defaults to the value used
#'   in the fit.
#' @param p_adjust Multiplicity adjustment method passed to
#'   \code{\link[stats]{p.adjust}}; default \code{"BH"}.
#' @param alpha Significance level applied to the adjusted probabilities.
#' @return A data frame of uniform and non-uniform DIF statistics per item and
#'   factor, with F statistics, raw and adjusted probabilities, and flags.
#' @examples
#' set.seed(1); n <- 600
#' d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
#' sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1
#' X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
#' colnames(X) <- paste0("I", 1:8)
#' dif_anova(rasch(X), factors = data.frame(group = g))
#' @export
dif_anova <- function(fit, factors = NULL, n_groups = NULL, p_adjust = "BH",
                      alpha = 0.05) {
  Z <- fit$residuals; L <- ncol(Z)
  factors <- .dif_factors(fit, factors)
  if (is.null(n_groups)) n_groups <- fit$n_groups
  ci <- .dif_class_intervals(fit, n_groups)

  res <- list()
  for (fname in names(factors)) {
    grp <- factor(factors[[fname]])
    out <- data.frame(factor = fname, item = colnames(Z),
                      F_uniform = NA_real_, p_uniform = NA_real_,
                      F_nonuniform = NA_real_, p_nonuniform = NA_real_)
    for (i in seq_len(L)) {
      d <- data.frame(z = Z[, i], g = grp, ci = ci)
      d <- d[stats::complete.cases(d), ]
      if (nrow(d) < 10 || length(unique(d$g)) < 2) next
      a <- tryCatch(stats::anova(stats::lm(z ~ g * ci, data = d)),
                    error = function(e) NULL)
      if (is.null(a)) next
      rn <- rownames(a)
      if ("g" %in% rn) {
        out$F_uniform[i] <- a["g", "F value"]; out$p_uniform[i] <- a["g", "Pr(>F)"]
      }
      if ("g:ci" %in% rn) {
        out$F_nonuniform[i] <- a["g:ci", "F value"]
        out$p_nonuniform[i] <- a["g:ci", "Pr(>F)"]
      }
    }
    out$p_uniform_adj <- p.adjust(out$p_uniform, method = p_adjust)
    out$p_nonuniform_adj <- p.adjust(out$p_nonuniform, method = p_adjust)
    out$uniform_DIF <- !is.na(out$p_uniform_adj) & out$p_uniform_adj < alpha
    out$nonuniform_DIF <- !is.na(out$p_nonuniform_adj) &
      out$p_nonuniform_adj < alpha
    res[[fname]] <- out
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

# variables of an ANOVA term label, e.g. "g1:ci" -> c("g1", "ci")
.term_vars <- function(term) strsplit(term, ":", fixed = TRUE)[[1]]

#' Factorial DIF analysis with Tukey comparisons
#'
#' Models all nominated person factors jointly: for each item the
#' standardised residuals are analysed by the full factorial of the person
#' factors crossed with the trait class interval,
#' \code{z ~ (f1 * f2 * ...) * ci}. Terms not involving the class interval
#' are uniform DIF effects (main effects and factor-by-factor interactions);
#' terms involving it are non-uniform. Probabilities are adjusted across
#' items within each term (Benjamini-Hochberg by default). A significant
#' interaction supersedes the main effects (and lower-order interactions) of
#' the factors it involves, which is recorded in the \code{superseded}
#' column; interpret the highest-order significant terms. Tukey HSD
#' comparisons are returned for every significant, non-superseded group term
#' (the cell-mean contrasts for interactions), with Tukey's own familywise
#' adjustment within each term.
#'
#' Sums of squares are sequential (factors in the order given, class
#' interval last), as is conventional for this residual diagnostic; with
#' markedly unbalanced groups the term order matters and the factor-at-a-time
#' \code{\link{dif_anova}} is a useful cross-check.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param factors As in \code{\link{dif_anova}}; at least one factor, usually
#'   two or more.
#' @param n_groups Number of trait class intervals; defaults to the value
#'   used in the fit.
#' @param p_adjust Multiplicity adjustment across items within each term;
#'   default \code{"BH"}.
#' @param alpha Significance level applied to the adjusted probabilities.
#' @param effects \code{"factorial"} (default) crosses every person factor
#'   with every other and with the class interval; \code{"main"} fits the
#'   factors additively (each factor's main effect and its interaction with
#'   the class interval, but no factor-by-factor terms).
#' @return A list with \code{terms}, the complete per-item analysis of
#'   variance table (term, df, sum of squares, mean square, F, raw and
#'   adjusted p, significance, supersession, including the residual row),
#'   and \code{tukey} (per item, term, and level comparison: difference,
#'   95 per cent interval, and Tukey-adjusted p), plus the \code{alpha} and
#'   adjustment used. Tukey comparisons are reported for significant,
#'   non-superseded group terms except two-level main effects, where the
#'   F test is already the only comparison.
#' @examples
#' set.seed(1); n <- 800
#' d <- seq(-1.5, 1.5, length.out = 6)
#' g1 <- rep(c("a", "b"), each = n / 2)
#' g2 <- rep(c("x", "y"), times = n / 2)
#' sh <- matrix(0, n, 6); sh[g1 == "b", 2] <- 0.8
#' X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
#' colnames(X) <- paste0("I", 1:6)
#' fit <- rasch(data.frame(X, g1 = g1, g2 = g2), factors = c("g1", "g2"))
#' dif_anova_factorial(fit)$terms
#' @export
dif_anova_factorial <- function(fit, factors = NULL, n_groups = NULL,
                                p_adjust = "BH", alpha = 0.05,
                                effects = c("factorial", "main")) {
  effects <- match.arg(effects)
  Z <- fit$residuals; L <- ncol(Z)
  factors <- .dif_factors(fit, factors)
  if (is.null(n_groups)) n_groups <- fit$n_groups
  ci <- .dif_class_intervals(fit, n_groups)

  fnames <- names(factors)
  safe <- paste0("f", seq_along(fnames))           # syntactic stand-ins
  op <- if (effects == "factorial") " * " else " + "
  form <- stats::as.formula(paste("z ~ (", paste(safe, collapse = op), ") * ci"))

  fits <- vector("list", L); rows <- list()
  for (i in seq_len(L)) {
    d <- data.frame(z = Z[, i], ci = ci)
    for (j in seq_along(fnames)) d[[safe[j]]] <- factor(factors[[fnames[j]]])
    d <- d[stats::complete.cases(d), ]
    if (nrow(d) < 10 || any(vapply(safe, function(s)
      length(unique(d[[s]])) < 2, TRUE))) next
    a <- tryCatch(stats::aov(form, data = d), error = function(e) NULL)
    if (is.null(a)) next
    fits[[i]] <- a
    sm <- summary(a)[[1]]
    rows[[length(rows) + 1L]] <- data.frame(
      item = colnames(Z)[i], term = trimws(rownames(sm)), df = sm$Df,
      sum_sq = sm$`Sum Sq`, mean_sq = sm$`Mean Sq`,
      F_value = sm$`F value`, p = sm$`Pr(>F)`)
  }
  if (!length(rows)) stop("no item yielded an estimable factorial ANOVA")
  terms <- do.call(rbind, rows)
  rownames(terms) <- NULL

  # adjust across items within each term (the residual rows carry no test)
  terms$p_adj <- NA_real_
  for (tt in setdiff(unique(terms$term), "Residuals")) {
    sel <- terms$term == tt
    terms$p_adj[sel] <- p.adjust(terms$p[sel], method = p_adjust)
  }
  terms$significant <- !is.na(terms$p_adj) & terms$p_adj < alpha

  # a significant higher-order interaction supersedes the lower-order terms
  # built from a subset of its variables (within the same item)
  terms$superseded <- FALSE
  for (it in unique(terms$item)) {
    sel <- which(terms$item == it & terms$significant)
    if (length(sel) < 2) next
    vlist <- lapply(terms$term[sel], .term_vars)
    for (a_i in seq_along(sel)) for (b_i in seq_along(sel)) {
      if (a_i == b_i) next
      if (all(vlist[[a_i]] %in% vlist[[b_i]]) &&
          length(vlist[[a_i]]) < length(vlist[[b_i]]))
        terms$superseded[sel[a_i]] <- TRUE
    }
  }

  # Tukey HSD for significant, non-superseded terms that do not involve the
  # class interval (the group structure itself)
  tk <- list()
  for (i in seq_len(L)) {
    a <- fits[[i]]; if (is.null(a)) next
    it <- colnames(Z)[i]
    cand <- terms[terms$item == it & terms$significant & !terms$superseded, ]
    # group terms only; and no comparisons for a two-level main effect,
    # where the F test is already the only contrast
    keep_t <- !vapply(cand$term, function(tt) "ci" %in% .term_vars(tt), TRUE) &
      !(cand$df == 1L & !grepl(":", cand$term, fixed = TRUE))
    cand <- cand$term[keep_t]
    if (!length(cand)) next
    th <- tryCatch(stats::TukeyHSD(a, which = cand), error = function(e) NULL)
    if (is.null(th)) next
    for (tt in names(th)) {
      tb <- as.data.frame(th[[tt]])
      tk[[length(tk) + 1L]] <- data.frame(
        item = it, term = tt, comparison = rownames(tb),
        difference = tb$diff, lower = tb$lwr, upper = tb$upr,
        p_tukey = tb$`p adj`, row.names = NULL)
    }
  }
  tukey <- if (length(tk)) do.call(rbind, tk) else
    data.frame(item = character(), term = character(),
               comparison = character(), difference = numeric(),
               lower = numeric(), upper = numeric(), p_tukey = numeric())

  # map the syntactic stand-ins back to the nominated factor names
  relabel <- function(x) {
    for (j in rev(seq_along(fnames)))
      x <- gsub(paste0("\\bf", j, "\\b"), fnames[j], x)
    x
  }
  terms$term <- relabel(terms$term)
  tukey$term <- relabel(tukey$term)

  list(terms = terms, tukey = tukey, alpha = alpha, p_adjust = p_adjust)
}
