# RaschR :: multiple choice
# ===========================================================================
# Scoring of raw multiple-choice responses against a key, with the raw
# responses retained for distractor analysis: per option, the count,
# proportion, mean person location, and point-biserial with the person
# measure. A distractor chosen by abler persons than the keyed option is the
# classic signature of a miskeyed item and is flagged. Option curves show
# the proportion choosing each option across trait class intervals.
# ===========================================================================

# Resolve a key given as a named vector (item -> correct response) or a data
# frame with columns item and key.
.resolve_key <- function(key) {
  if (is.data.frame(key)) {
    if (!all(c("item", "key") %in% names(key)))
      stop("a key data frame needs columns item, key")
    key <- setNames(as.character(key$key), as.character(key$item))
  }
  if (is.null(names(key))) stop("the key must be named by item")
  setNames(trimws(toupper(as.character(key))), names(key))
}

# Score raw responses 0/1 against the key. Blank, NA, and missing-data codes
# become NA. Matching is case-insensitive after trimming.
.score_mc <- function(X, key) {
  keyed <- intersect(names(key), colnames(X))
  if (!length(keyed)) stop("no key item matches an item column")
  raw <- matrix(trimws(toupper(as.character(X[, keyed]))), nrow(X),
                length(keyed), dimnames = list(NULL, keyed))
  raw[raw %in% c("", "NA", "-1")] <- NA
  scored <- matrix(NA_integer_, nrow(raw), ncol(raw),
                   dimnames = dimnames(raw))
  for (j in seq_along(keyed))
    scored[, j] <- ifelse(is.na(raw[, j]), NA_integer_,
                          as.integer(raw[, j] == key[keyed[j]]))
  list(scored = scored, raw = raw, key = key[keyed])
}

#' Distractor analysis for multiple-choice items
#'
#' For every keyed item and response option: the count and proportion
#' choosing it, the mean location of those persons, and the point-biserial
#' correlation between choosing the option and the person measure. Locations
#' and correlations use the rest measure (the person estimate from the other
#' items), so the analysed item cannot credit its own takers. The keyed
#' option should attract the ablest persons and carry the only positive
#' point-biserial; a distractor whose takers are abler than the keyed
#' option's (with at least \code{min_n} takers) is flagged as a possible
#' miskey.
#'
#' @param fit A fitted object from \code{\link{rasch}} run with a \code{key}.
#' @param items Optional subset of item names; defaults to every keyed item.
#' @param min_n Minimum takers for an option to be eligible for the miskey
#'   flag.
#' @return A data frame with one row per item-option: \code{item},
#'   \code{option}, \code{keyed}, \code{n}, \code{prop},
#'   \code{mean_location}, \code{point_biserial}, and \code{flag}.
#' @examples
#' set.seed(1); Np <- 400
#' th <- rnorm(Np)
#' raw <- sapply(seq(-1, 1, length.out = 6), function(d) {
#'   ok <- rbinom(Np, 1, plogis(th - d))
#'   ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
#' })
#' colnames(raw) <- paste0("M", 1:6)
#' fit <- rasch(raw, key = setNames(rep("A", 6), colnames(raw)))
#' head(distractor_analysis(fit))
#' @export
distractor_analysis <- function(fit, items = NULL, min_n = 10) {
  if (is.null(fit$mc)) stop("the fit has no key: run rasch(..., key = )")
  raw <- fit$mc$raw; key <- fit$mc$key
  if (is.null(items)) items <- colnames(raw)
  items <- intersect(items, colnames(raw))
  out <- list()
  for (it in items) {
    r <- raw[, it]
    idx <- match(it, colnames(fit$X))
    th <- .person_estimates(fit$X[, -idx, drop = FALSE],
                            fit$tau_list[-idx])$theta   # rest measure
    ok <- !is.na(r) & !is.na(th)
    opts <- sort(unique(r[ok]))
    rows <- data.frame(item = it, option = opts, keyed = opts == key[it],
                       n = NA_integer_, prop = NA_real_,
                       mean_location = NA_real_, point_biserial = NA_real_)
    for (j in seq_along(opts)) {
      sel <- ok & r == opts[j]
      rows$n[j] <- sum(sel)
      rows$prop[j] <- sum(sel) / sum(ok)
      rows$mean_location[j] <- mean(th[sel])
      ind <- as.integer(r[ok] == opts[j])
      rows$point_biserial[j] <- if (var(ind) > 0)
        cor(ind, th[ok]) else NA_real_
    }
    key_mean <- rows$mean_location[rows$keyed]
    rows$flag <- if (length(key_mean) == 1 && !is.na(key_mean))
      !rows$keyed & rows$n >= min_n & rows$mean_location > key_mean
    else rep(FALSE, nrow(rows))
    out[[it]] <- rows
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Plot multiple-choice option curves
#'
#' The proportion choosing each response option across class intervals of the
#' rest measure (the person estimate from the other items), with the keyed
#' option drawn solid and bold. The keyed option should rise with the trait
#' and every distractor should fall; a rising distractor is the graphical
#' signature of a miskey or an ambiguous option.
#'
#' @param fit A fitted object from \code{\link{rasch}} run with a \code{key}.
#' @param item Keyed item name.
#' @param n_groups Number of class intervals.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1); Np <- 400
#' th <- rnorm(Np)
#' raw <- sapply(seq(-1, 1, length.out = 6), function(d) {
#'   ok <- rbinom(Np, 1, plogis(th - d))
#'   ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
#' })
#' colnames(raw) <- paste0("M", 1:6)
#' fit <- rasch(raw, key = setNames(rep("A", 6), colnames(raw)))
#' plot_distractors(fit, "M3")
#' @export
plot_distractors <- function(fit, item, n_groups = fit$n_groups) {
  if (is.null(fit$mc)) stop("the fit has no key: run rasch(..., key = )")
  if (!item %in% colnames(fit$mc$raw)) stop("no such keyed item: ", item)
  r <- fit$mc$raw[, item]
  idx <- match(item, colnames(fit$X))
  th <- .person_estimates(fit$X[, -idx, drop = FALSE],
                          fit$tau_list[-idx])$theta     # rest measure
  ok <- !is.na(r) & !is.na(th)
  ng <- min(n_groups, max(2, floor(sum(ok) / 25)))
  ci <- cut(rank(th[ok], ties.method = "first"), ng, labels = FALSE)
  mid <- tapply(th[ok], ci, mean)
  opts <- sort(unique(r[ok]))
  key_opt <- fit$mc$key[item]
  op <- .rr_canvas(range(mid) + c(-0.2, 0.2), c(0, 1),
                   "Person location (logits)", "Proportion choosing option",
                   paste0("Option curves \u2013 ", item,
                          "  (key: ", key_opt, ")"))
  on.exit(par(op))
  for (j in seq_along(opts)) {
    pr <- tapply(r[ok] == opts[j], ci, mean)
    colr <- .rr$pal[(j - 1L) %% length(.rr$pal) + 1L]
    keyed <- opts[j] == key_opt
    lines(mid, pr, lwd = if (keyed) 3.2 else 1.8,
          lty = if (keyed) 1 else 5, col = colr)
    points(mid, pr, pch = 21, bg = colr, col = "white",
           cex = if (keyed) 1.5 else 1.1, lwd = 1)
  }
  .rr_legend("right", paste0(opts, ifelse(opts == key_opt, " (key)", "")),
             lwd = ifelse(opts == key_opt, 3.2, 1.8),
             lty = ifelse(opts == key_opt, 1, 5),
             col = .rr$pal[seq_along(opts)])
  invisible(NULL)
}
