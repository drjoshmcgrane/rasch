# rasch :: multiple choice
# ===========================================================================
# Scoring of raw multiple-choice responses against a key, with the raw
# responses retained for distractor analysis: per option, the count,
# proportion, mean person location, and point-biserial with the person
# measure. A distractor chosen by abler persons than the keyed option is the
# classic signature of a miskeyed item and is flagged. Option curves show
# the proportion choosing each option across trait class intervals.
# Three key forms are supported:
# a single correct option per item (scored 0/1), double keying (several
# options separated by "/" all scoring 1), and full polytomous option
# scoring (a data frame item/option/score assigning partial credit to
# informative distractors; Andrich & Styles 2011). distractor_rescore()
# proposes such a scoring from the rest-measure evidence.
# ===========================================================================

# Resolve a key into a per-item scoring map (named integer vector
# option -> score). Accepted forms: a named vector or data frame
# (item, key) where key is the correct option, possibly several separated
# by "/" (double keying, all scoring 1); or a data frame
# (item, option, score) assigning an integer score to every credited
# option (unlisted options score 0). Matching is case-insensitive.
.resolve_key <- function(key) {
  if (is.data.frame(key) && all(c("item", "option", "score") %in% names(key))) {
    key$item <- as.character(key$item)
    key$option <- trimws(toupper(as.character(key$option)))
    key$score <- as.integer(key$score)
    if (anyNA(key$score) || any(key$score < 0))
      stop("option scores must be non-negative integers")
    map <- lapply(split(key, key$item), function(d) {
      if (anyDuplicated(d$option))
        stop("duplicate option row(s) for item ", d$item[1])
      if (max(d$score) < 1)
        stop("item ", d$item[1], " credits no option")
      setNames(d$score, d$option)
    })
    return(map)
  }
  if (is.data.frame(key)) {
    if (!all(c("item", "key") %in% names(key)))
      stop("a key data frame needs columns item, key ",
           "(or item, option, score for polytomous option scoring)")
    key <- setNames(as.character(key$key), as.character(key$item))
  }
  if (is.null(names(key))) stop("the key must be named by item")
  lapply(setNames(trimws(toupper(as.character(key))), names(key)),
         function(k) {
           opts <- trimws(strsplit(k, "/", fixed = TRUE)[[1]])
           opts <- opts[nzchar(opts)]
           if (!length(opts)) stop("empty key entry")
           setNames(rep(1L, length(opts)), opts)
         })
}

# One-line display form of an item's scoring map (e.g. "C", "A/C", or
# "C=2, B=1").
.key_label <- function(m) {
  m <- sort(m[m > 0], decreasing = TRUE)
  if (all(m == 1L)) paste(names(m), collapse = "/")
  else paste(sprintf("%s=%d", names(m), m), collapse = ", ")
}

# Score raw responses against the scoring maps. Observed options absent
# from an item's map score 0; blank, NA, and missing-data codes become NA.
.score_mc <- function(X, map) {
  keyed <- intersect(names(map), colnames(X))
  if (!length(keyed)) stop("no key item matches an item column")
  raw <- matrix(trimws(toupper(as.character(X[, keyed]))), nrow(X),
                length(keyed), dimnames = list(NULL, keyed))
  raw[raw %in% c("", "NA", "-1")] <- NA
  scored <- matrix(NA_integer_, nrow(raw), ncol(raw),
                   dimnames = dimnames(raw))
  for (j in seq_along(keyed)) {
    m <- map[[keyed[j]]]
    s <- unname(m[raw[, j]])
    s[is.na(s) & !is.na(raw[, j])] <- 0L   # observed but uncredited option
    scored[, j] <- s
  }
  list(scored = scored, raw = raw, map = map[keyed],
       key = vapply(map[keyed], .key_label, ""))
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
#'   \code{option}, its assigned \code{score}, \code{keyed} (full credit),
#'   \code{n}, \code{prop}, \code{mean_location}, \code{point_biserial},
#'   and \code{flag}.
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
  raw <- fit$mc$raw; map <- fit$mc$map
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
    m <- map[[it]]
    sc <- unname(m[opts]); sc[is.na(sc)] <- 0L
    rows <- data.frame(item = it, option = opts, score = sc,
                       keyed = sc == max(m),
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
    key_mean <- max(rows$mean_location[rows$keyed], na.rm = TRUE)
    rows$flag <- if (is.finite(key_mean))
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
  m <- fit$mc$map[[item]]
  sc <- unname(m[opts]); sc[is.na(sc)] <- 0L
  keyed_v <- sc == max(m)
  op <- .rr_canvas(range(mid) + c(-0.2, 0.2), c(0, 1),
                   "Person location (logits)", "Proportion choosing option",
                   paste0(item, "  (key: ", fit$mc$key[item], ")"))
  on.exit(par(op))
  for (j in seq_along(opts)) {
    pr <- tapply(r[ok] == opts[j], ci, mean)
    colr <- .rr$pal[(j - 1L) %% length(.rr$pal) + 1L]
    lines(mid, pr, lwd = if (keyed_v[j]) 3.2 else if (sc[j] > 0) 2.4 else 1.8,
          lty = if (keyed_v[j]) 1 else if (sc[j] > 0) 2 else 5, col = colr)
    points(mid, pr, pch = 21, bg = colr, col = "white",
           cex = if (keyed_v[j]) 1.5 else 1.1, lwd = 1)
  }
  labs <- paste0(opts,
                 ifelse(keyed_v, " (key)",
                        ifelse(sc > 0, sprintf(" (credit %d)", sc), "")))
  .rr_legend("right", labs,
             lwd = ifelse(keyed_v, 3.2, ifelse(sc > 0, 2.4, 1.8)),
             lty = ifelse(keyed_v, 1, ifelse(sc > 0, 2, 5)),
             col = .rr$pal[seq_along(opts)])
  invisible(NULL)
}

#' Propose polytomous option scores from the distractor evidence
#'
#' Multiple-choice items can be rescored polytomously so
#' that a distractor carrying information about the trait receives partial
#' credit (Andrich and Styles 2011). This function proposes such a scoring
#' from the rest-measure distractor analysis: within each keyed item, a
#' distractor qualifies for credit when it attracts at least \code{min_n}
#' takers, its takers' mean rest location exceeds that of the uncredited
#' distractors by more than \code{z} standard errors of the difference,
#' and it remains below the keyed option. Qualifying distractors are
#' ranked by mean location and scored 1, 2, ... below the keyed option's
#' top score. The result is a proposal for substantive review, not an
#' automatic decision: inspect \code{\link{plot_distractors}} and the item
#' content, edit as needed, then refit with
#' \code{rasch(raw_data, key = proposal$option_scores)}.
#'
#' @param fit A fitted object from \code{\link{rasch}} run with a \code{key}.
#' @param items Optional subset of keyed item names.
#' @param min_n Minimum takers for a distractor to be considered.
#' @param z Required separation, in standard errors, between a credited
#'   distractor and the uncredited ones.
#' @return A list of class \code{"rasch_rescore"}: \code{option_scores}, a
#'   data frame (\code{item}, \code{option}, \code{score}) ready for
#'   \code{rasch(key = )} and covering every observed option of the
#'   examined items, and \code{evidence}, the distractor analysis with the
#'   proposed scores and the separation z per option.
#' @references Andrich, D. and Styles, I. (2011). Distractors with
#'   information in multiple choice items: A rationale based on the Rasch
#'   model. Journal of Applied Measurement, 12, 67-95.
#' @examples
#' set.seed(1); Np <- 600
#' th <- rnorm(Np)
#' raw <- sapply(seq(-0.5, 0.5, length.out = 4), function(d) {
#'   x <- vapply(th, function(b) sample(0:2, 1,
#'     prob = item_moments(b, c(d - 0.7, d + 0.7))$P), 0L)
#'   c("D", "B", "A")[x + 1]   # B is an informative distractor
#' })
#' colnames(raw) <- paste0("M", 1:4)
#' fit <- rasch(raw, key = setNames(rep("A", 4), colnames(raw)))
#' pr <- distractor_rescore(fit)
#' pr$option_scores
#' @export
distractor_rescore <- function(fit, items = NULL, min_n = 20, z = 1.96) {
  da <- distractor_analysis(fit, items = items, min_n = min_n)
  raw <- fit$mc$raw
  ev <- list(); os <- list()
  for (it in unique(da$item)) {
    d <- da[da$item == it, ]
    idx <- match(it, colnames(fit$X))
    th <- .person_estimates(fit$X[, -idx, drop = FALSE],
                            fit$tau_list[-idx])$theta
    r <- raw[, it]; ok <- !is.na(r) & !is.na(th)
    # per-option SE of the mean rest location
    d$se_location <- vapply(d$option, function(o) {
      x <- th[ok & r == o]
      if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
    }, 0)
    key_row <- which(d$keyed)[which.max(d$mean_location[d$keyed])]
    cand <- which(!d$keyed & d$n >= min_n)
    # the uncredited baseline: distractors not currently under consideration
    d$z_sep <- NA_real_
    credited <- integer(0)
    for (j in cand) {
      others <- setdiff(cand, j)
      if (!length(others)) next
      base <- th[ok & r %in% d$option[others]]
      xj <- th[ok & r == d$option[j]]
      if (length(base) < 2 || length(xj) < 2) next
      sep <- (mean(xj) - mean(base)) /
        sqrt(var(xj) / length(xj) + var(base) / length(base))
      d$z_sep[j] <- sep
      if (sep > z && mean(xj) < d$mean_location[key_row]) credited <- c(credited, j)
    }
    credited <- credited[order(d$mean_location[credited])]
    d$proposed <- 0L
    if (length(credited)) d$proposed[credited] <- seq_along(credited)
    d$proposed[d$keyed] <- length(credited) + 1L
    ev[[it]] <- d
    os[[it]] <- data.frame(item = it, option = d$option, score = d$proposed)
  }
  out <- list(option_scores = do.call(rbind, os),
              evidence = do.call(rbind, ev))
  rownames(out$option_scores) <- rownames(out$evidence) <- NULL
  class(out) <- "rasch_rescore"
  out
}

#' @export
print.rasch_rescore <- function(x, ...) {
  n_credit <- sum(x$option_scores$score > 0 &
                  x$option_scores$score < ave(x$option_scores$score,
                                              x$option_scores$item,
                                              FUN = max))
  cat(sprintf("Polytomous option-scoring proposal (Andrich & Styles 2011): %d distractor(s) credited\n",
              n_credit))
  ev <- x$evidence
  num <- vapply(ev, is.numeric, TRUE)
  ev[num] <- lapply(ev[num], round, 3)
  print(ev, row.names = FALSE)
  cat("review substantively, edit, then refit: rasch(raw, key = proposal$option_scores)\n")
  invisible(x)
}
