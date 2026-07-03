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
#'   (\code{item}, \code{n}, \code{facility}, \code{item_total},
#'   \code{item_rest}, \code{di}, \code{alpha_drop}), and the scalars
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
  X <- fit$X[stats::complete.cases(fit$X), , drop = FALSE]
  if (nrow(X) < 3) stop("fewer than 3 complete cases: no traditional statistics")
  L <- ncol(X); n <- nrow(X)
  tot <- rowSums(X)
  thirds <- cut(rank(tot, ties.method = "first"), 3, labels = FALSE)
  alpha <- fit$alpha$alpha
  tab <- data.frame(item = colnames(X), n = n,
                    facility = colMeans(X) / fit$m,
                    item_total = NA_real_, item_rest = NA_real_,
                    di = NA_real_, alpha_drop = NA_real_)
  for (i in seq_len(L)) {
    x <- X[, i]
    if (sd(x) > 0) {
      tab$item_total[i] <- cor(x, tot)
      tab$item_rest[i] <- cor(x, tot - x)
      tab$di[i] <- mean(x[thirds == 3]) / fit$m[i] -
                   mean(x[thirds == 1]) / fit$m[i]
    }
    Xr <- X[, -i, drop = FALSE]
    vr <- var(rowSums(Xr))
    if (L > 2 && vr > 0)
      tab$alpha_drop[i] <- (L - 1) / (L - 2) * (1 - sum(apply(Xr, 2, var)) / vr)
  }
  rownames(tab) <- NULL
  out <- list(table = tab, alpha = alpha, n = n, mean = mean(tot),
              sd = sd(tot),
              sem = if (is.finite(alpha)) sd(tot) * sqrt(1 - alpha) else NA_real_)
  class(out) <- "rmt_ctt"
  out
}

#' @export
print.rmt_ctt <- function(x, ...) {
  cat(sprintf("Traditional statistics (complete cases, n = %d)\n", x$n))
  cat(sprintf("Raw score mean %.2f, SD %.2f; alpha %.3f; SEM %.2f (one value for all persons)\n",
              x$mean, x$sd, x$alpha, x$sem))
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
