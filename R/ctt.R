# rmt :: traditional (classical test theory) statistics
# ===========================================================================
# Classical test theory companion statistics, computed on complete cases
# only (alpha and its relatives have no missing-data form): item facility,
# item-total and item-rest correlations,
# the upper-lower discrimination index over total-score thirds (Andrich &
# Marais 2019 ch. 5), coefficient alpha with alpha-if-item-deleted, and the
# classical summary (mean, SD, SEM) for comparison with the Rasch results.
# ===========================================================================

#' Traditional (classical test theory) statistics
#'
#' The classical companion table conventionally reported alongside a Rasch
#' analysis (Andrich and Marais 2019, chs. 3-5), on complete cases only: per item the facility (mean score over
#' maximum), the item-total and corrected item-rest correlations, the
#' discrimination index DI = PRU - PRL (mean proportion-of-maximum in the
#' upper third of total scores minus the lower third), and alpha if the
#' item is deleted; plus coefficient alpha, the raw-score mean, SD, and the
#' classical standard error of measurement \eqn{s\sqrt{1 - \alpha}}, which
#' unlike the Rasch SE is one value for all persons.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return A list of class \code{"rmt_ctt"}: the per-item \code{table}
#'   (\code{item}, \code{n}, \code{min}, \code{max}, \code{facility},
#'   \code{item_total}, \code{item_rest}, \code{di}, \code{alpha_drop}), and
#'   the scalars
#'   \code{alpha}, \code{n} (complete cases), \code{mean}, \code{sd}, and
#'   \code{sem}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#' colnames(X) <- paste0("I", 1:8)
#' ctt_table(rasch(X))
#' @export
ctt_table <- function(fit) {
  X <- fit$X
  L <- ncol(X)
  n_i <- colSums(!is.na(X))
  if (all(n_i < 3)) stop("fewer than 3 responses per item: no traditional statistics")
  # available-case statistics: each person's total is their proportion of
  # the maximum over the items they answered, so persons with different
  # answered sets remain comparable; with complete data every statistic
  # reduces exactly to its textbook complete-case form
  Mmat <- matrix(rep(fit$m, each = nrow(X)), nrow(X), L)
  Mmat[is.na(X)] <- NA
  tot_p <- rowSums(X, na.rm = TRUE) / rowSums(Mmat, na.rm = TRUE)
  thirds <- cut(rank(tot_p, ties.method = "first"), 3, labels = FALSE)
  alpha <- fit$alpha$alpha
  # pairwise covariance carries the alpha-if-deleted computation under
  # missing data (equal to the variance form when data are complete)
  C <- suppressWarnings(stats::cov(X, use = "pairwise.complete.obs"))
  min_i <- suppressWarnings(vapply(seq_len(L), function(i)
    if (n_i[i] == 0) NA_real_ else min(X[, i], na.rm = TRUE), 0))
  tab <- data.frame(item = colnames(X), n = n_i, min = min_i, max = fit$m,
                    facility = colMeans(X, na.rm = TRUE) / fit$m,
                    item_total = NA_real_, item_rest = NA_real_,
                    di = NA_real_, alpha_drop = NA_real_)
  for (i in seq_len(L)) {
    x <- X[, i]; ok <- !is.na(x)
    if (sum(ok) >= 3 && stats::sd(x[ok]) > 0) {
      rest_p <- (rowSums(X, na.rm = TRUE) - ifelse(ok, x, 0)) /
        pmax(rowSums(Mmat, na.rm = TRUE) - ifelse(ok, fit$m[i], 0), 1)
      tab$item_total[i] <- .safe_cor(x[ok], tot_p[ok])
      tab$item_rest[i] <- .safe_cor(x[ok], rest_p[ok])
      hi <- ok & thirds == 3; lo <- ok & thirds == 1
      if (sum(hi) >= 2 && sum(lo) >= 2)
        tab$di[i] <- mean(x[hi]) / fit$m[i] - mean(x[lo]) / fit$m[i]
    }
    if (L > 2) {
      Cr <- C[-i, -i, drop = FALSE]
      sr <- sum(Cr, na.rm = TRUE)
      if (is.finite(sr) && sr > 0)
        tab$alpha_drop[i] <- (L - 1) / (L - 2) *
          (1 - sum(diag(Cr), na.rm = TRUE) / sr)
    }
  }
  rownames(tab) <- NULL
  cc <- stats::complete.cases(X)
  tot_cc <- rowSums(X[cc, , drop = FALSE])
  out <- list(table = tab, alpha = alpha, n = sum(cc),
              n_range = range(n_i),
              mean = if (sum(cc) >= 3) mean(tot_cc) else NA_real_,
              sd = if (sum(cc) >= 3) stats::sd(tot_cc) else NA_real_,
              sem = if (is.finite(alpha) && sum(cc) >= 3)
                stats::sd(tot_cc) * sqrt(1 - alpha) else NA_real_)
  class(out) <- "rmt_ctt"
  out
}

#' @export
print.rmt_ctt <- function(x, ...) {
  cat(sprintf("Traditional statistics (available cases; item n %d-%d; %d complete)\n",
              x$n_range[1], x$n_range[2], x$n))
  if (is.finite(x$mean))
    cat(sprintf("Raw score mean %.2f, SD %.2f (complete responders); alpha %.3f; SEM %.2f\n",
                x$mean, x$sd, x$alpha, x$sem))
  else
    cat(sprintf("Too few complete responders for total-score summaries; alpha %.3f\n",
                x$alpha))
  y <- x$table
  num <- vapply(y, is.numeric, TRUE)
  y[num] <- lapply(y[num], round, 3)
  print(y, row.names = FALSE)
  invisible(x)
}

#' Reshape repeated measurements for racked or stacked analysis
#'
#' Repeated measurements (the same persons and items at two or more time
#' points) enter a Rasch analysis in one of two designs (Andrich & Marais
#' 2019, ch. 26). \emph{Racking} keeps one row per person and duplicates
#' the items per time point (columns \code{item@time}), so change over time
#' shows in the item estimates. \emph{Stacking} keeps one column per item
#' and duplicates the persons per time point (rows \code{person@time}), so
#' change shows in the person estimates and DIF of items over time can be
#' examined with \code{time} as a person factor.
#'
#' @param data A long data frame with one measurement per row.
#' @param person,time Names of the person and time-point columns.
#' @param items Character vector naming the item columns.
#' @return \code{rack_data}: a wide data frame with one row per person and
#'   \code{length(items) * n_times} item columns. \code{stack_data}: a data
#'   frame with one row per person-time (\code{id} column), the original
#'   item columns, and \code{time} as a factor column for DIF analysis.
#' @examples
#' d <- data.frame(pid = rep(1:100, 2), t = rep(1:2, each = 100),
#'                 Q1 = rbinom(200, 1, 0.6), Q2 = rbinom(200, 1, 0.5))
#' racked <- rack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
#' names(racked)
#' stacked <- stack_data(d, person = "pid", time = "t", items = c("Q1", "Q2"))
#' head(stacked)
#' @export
rack_data <- function(data, person, time, items) {
  data <- as.data.frame(data)
  for (col in c(person, time)) if (!col %in% names(data))
    stop("column not found: ", col)
  bad <- setdiff(items, names(data))
  if (length(bad)) stop("item column(s) not found: ", paste(bad, collapse = ", "))
  times <- sort(unique(data[[time]]))
  ids <- unique(data[[person]])
  out <- data.frame(id = ids)
  for (tt in times) {
    d_t <- data[data[[time]] == tt, , drop = FALSE]
    if (anyDuplicated(d_t[[person]]))
      stop("more than one row for a person at time ", tt)
    idx <- match(ids, d_t[[person]])
    blk <- d_t[idx, items, drop = FALSE]
    names(blk) <- paste0(items, "@", tt)
    out <- cbind(out, blk)
  }
  rownames(out) <- NULL
  out
}

#' @rdname rack_data
#' @export
stack_data <- function(data, person, time, items) {
  data <- as.data.frame(data)
  for (col in c(person, time)) if (!col %in% names(data))
    stop("column not found: ", col)
  bad <- setdiff(items, names(data))
  if (length(bad)) stop("item column(s) not found: ", paste(bad, collapse = ", "))
  out <- data.frame(id = paste0(data[[person]], "@", data[[time]]),
                    time = factor(data[[time]]),
                    data[, items, drop = FALSE])
  rownames(out) <- NULL
  out
}
