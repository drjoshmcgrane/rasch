# rmt :: Bradley-Terry-Luce paired comparisons
# ===========================================================================
# The Bradley-Terry-Luce model (Bradley & Terry 1952; Luce 1959) for paired
# comparisons: Pr{a beats b} = exp(beta_a - beta_b) / (1 + exp(...)). This
# is the conditional form of the dichotomous Rasch model (Rasch 1961;
# Andrich 1978): given that one of two items is answered correctly, the
# probability the harder one was is exactly of BTL form, and the package's
# pairwise conditional estimation maximises precisely such a likelihood on
# the pair-conditional counts. BTL therefore sits inside the same family
# with the person parameter replaced by exchangeable comparison
# replications, and it is estimated here by the same conventions as
# everything else: Newton-Raphson on the conditional likelihood, sum-zero
# identification, and Godambe sandwich standard errors, clustered by judge
# when judges are identified. Fit follows the same residual logic as the
# rest of the package (Andrich & Marais 2019, ch. 23): per
# comparison z = (y - P)/sqrt(PQ); objects and judges carry the
# log-of-mean-square fit residual over their comparisons with apportioned
# degrees of freedom, and the classical pairwise chi-square compares the
# observed and expected win proportions of every pair.
# ===========================================================================

# union-find connectivity over the comparison graph
.btl_components <- function(K, ia, ib) {
  parent <- seq_len(K)
  find <- function(x) { while (parent[x] != x) x <- parent[x]; x }
  for (r in seq_along(ia)) {
    ra <- find(ia[r]); rb <- find(ib[r])
    if (ra != rb) parent[ra] <- rb
  }
  vapply(seq_len(K), find, 1L)
}

#' Fit the Bradley-Terry-Luce model to paired comparisons
#'
#' Estimates object locations from paired-comparison data by conditional
#' maximum likelihood. The Bradley-Terry-Luce model is the conditional form
#' of the dichotomous Rasch model -- within an item pair, given one correct
#' response, the Rasch probability that it was the easier item is exactly
#' of BTL form -- so it belongs to the same measurement family and is
#' estimated by the same conventions as the rest of the package:
#' Newton-Raphson on the person-free likelihood, locations identified by
#' the sum-zero constraint, and Godambe sandwich standard errors, clustered
#' by judge when a judge column is given (so repeated comparisons by the
#' same judge need not be independent). Objects that win or lose every
#' comparison have no finite estimate and are removed with a note, exactly
#' as extreme persons are set aside in a Rasch calibration; the comparison
#' graph must remain connected.
#'
#' Fit is reported at three levels, mirroring the Rasch diagnostics.
#' Per object and (when given) per judge: the log-of-mean-square fit
#' residual of Andrich and Marais (2019, ch. 23) over their comparisons, with apportioned degrees of freedom --
#' an erratic judge or an object of inconsistent quality shows exactly as
#' an erratic person or misfitting item does. Per pair: the classical
#' goodness-of-fit table comparing observed and expected win proportions,
#' with the total chi-square on (pairs used) - (objects - 1) degrees of
#' freedom. The object separation index is the analogue of the PSI:
#' the proportion of observed location variance not due to error.
#'
#' @param data A data frame with one comparison per row.
#' @param object_a,object_b Names of the columns holding the two objects
#'   compared.
#' @param winner Name of the column holding the winner of each row: its
#'   value must equal the row's \code{object_a} or \code{object_b} entry
#'   (anything else is a tie).
#' @param judge Optional name of a judge column; enables the judge fit
#'   table and clusters the sandwich standard errors by judge.
#' @param count Optional name of a column of replication counts (a row
#'   standing for several identical comparisons).
#' @param ties How to treat ties: \code{"drop"} (default, removed with a
#'   note), \code{"half"} (half a win each way, a common pragmatic device
#'   -- flagged in the notes because the halves are not independent
#'   Bernoulli trials), or \code{"error"}.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return A list of class \code{"rmt_btl"}: \code{objects} (location, se,
#'   comparisons, wins, outfit mean square, fit residual and its df),
#'   \code{pairs} (per pair: n, observed and expected win proportions,
#'   standardised residual, chi-square component), \code{judges} (when
#'   given: per judge n, outfit, fit residual, df), \code{total_chisq},
#'   \code{total_df}, \code{total_p}, the object separation index
#'   \code{osi}, \code{loglik}, convergence details, and \code{notes}.
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#'   incomplete block designs: I. The method of paired comparisons.
#'   Biometrika, 39, 324-345. Luce, R. D. (1959). Individual Choice
#'   Behavior. Wiley. Andrich, D. (1978). Relationships between the
#'   Thurstone and Rasch approaches to item scaling. Applied Psychological
#'   Measurement, 2, 451-462.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pairs <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pairs[, 1], each = 30),
#'                 b = rep(pairs[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' btl(d, object_a = "a", object_b = "b", winner = "win")
#' @export
btl <- function(data, object_a, object_b, winner, judge = NULL, count = NULL,
                ties = c("drop", "half", "error"), maxit = 60, tol = 1e-8) {
  ties <- match.arg(ties)
  data <- as.data.frame(data)
  for (col in c(object_a, object_b, winner, judge, count))
    if (!col %in% names(data)) stop("column not found: ", col)
  a <- trimws(as.character(data[[object_a]]))
  b <- trimws(as.character(data[[object_b]]))
  wn <- trimws(as.character(data[[winner]]))
  jd <- if (is.null(judge)) NULL else as.character(data[[judge]])
  w <- if (is.null(count)) rep(1, nrow(data)) else as.numeric(data[[count]])
  keep <- !is.na(a) & !is.na(b) & !is.na(wn) & a != b & !is.na(w) & w > 0
  notes <- character(0)
  if (any(!keep)) {
    notes <- c(notes, sprintf("%d row(s) dropped (missing, zero-count, or self-comparison)",
                              sum(!keep)))
    a <- a[keep]; b <- b[keep]; wn <- wn[keep]; w <- w[keep]
    if (!is.null(jd)) jd <- jd[keep]
  }
  if (!length(a)) stop("no usable comparisons")

  # outcome: 1 = a wins, 0 = b wins, NA = tie
  y <- ifelse(wn == a, 1, ifelse(wn == b, 0, NA))
  if (anyNA(y)) {
    n_tie <- sum(is.na(y))
    if (ties == "error") stop(n_tie, " tie(s) present; set ties = 'drop' or 'half'")
    if (ties == "drop") {
      notes <- c(notes, sprintf("%d tie(s) dropped", n_tie))
      sel <- !is.na(y)
      a <- a[sel]; b <- b[sel]; y <- y[sel]; w <- w[sel]
      if (!is.null(jd)) jd <- jd[sel]
    } else {
      notes <- c(notes, sprintf("%d tie(s) scored half a win each way (halves are not independent trials)",
                                n_tie))
      t_i <- which(is.na(y))
      a <- c(a, a[t_i]); b <- c(b, b[t_i])
      y[t_i] <- 1; y <- c(y, rep(0, length(t_i)))
      w[t_i] <- w[t_i] / 2; w <- c(w, w[t_i])
      if (!is.null(jd)) jd <- c(jd, jd[t_i])
    }
  }

  # remove undefeated / winless objects iteratively (no finite estimate),
  # as extreme persons are set aside in a Rasch calibration
  repeat {
    objs <- sort(unique(c(a, b)))
    win_of <- setNames(numeric(length(objs)), objs)
    n_of <- win_of
    for (r in seq_along(a)) {
      win_of[a[r]] <- win_of[a[r]] + w[r] * y[r]
      win_of[b[r]] <- win_of[b[r]] + w[r] * (1 - y[r])
      n_of[a[r]] <- n_of[a[r]] + w[r]; n_of[b[r]] <- n_of[b[r]] + w[r]
    }
    ext <- names(win_of)[win_of == 0 | win_of == n_of]
    if (!length(ext)) break
    notes <- c(notes, sprintf("object(s) with no wins or no losses removed (no finite estimate): %s",
                              paste(ext, collapse = ", ")))
    sel <- !(a %in% ext) & !(b %in% ext)
    a <- a[sel]; b <- b[sel]; y <- y[sel]; w <- w[sel]
    if (!is.null(jd)) jd <- jd[sel]
    if (!length(a)) stop("no comparisons remain after removing extreme objects")
  }
  objs <- sort(unique(c(a, b)))
  K <- length(objs)
  if (K < 3) stop("need at least three comparable objects")
  ia <- match(a, objs); ib <- match(b, objs)

  comp <- .btl_components(K, ia, ib)
  if (length(unique(comp)) > 1) {
    parts <- split(objs, comp)
    stop("the comparison graph is disconnected; components: ",
         paste(vapply(parts, paste, "", collapse = ","), collapse = " | "))
  }

  # Newton-Raphson on the sum-zero design (as pcml identifies item locations)
  B <- rbind(diag(K - 1L), rep(-1, K - 1L))
  beta <- numeric(K)
  for (it in seq_len(maxit)) {
    eta <- beta[ia] - beta[ib]
    P <- plogis(eta); V <- pmax(P * (1 - P), 1e-12)
    res <- w * (y - P)
    g_full <- vapply(seq_len(K), function(k)
      sum(res[ia == k]) - sum(res[ib == k]), 0)
    H_full <- matrix(0, K, K)
    for (r in seq_along(ia)) {
      hv <- w[r] * V[r]
      H_full[ia[r], ia[r]] <- H_full[ia[r], ia[r]] + hv
      H_full[ib[r], ib[r]] <- H_full[ib[r], ib[r]] + hv
      H_full[ia[r], ib[r]] <- H_full[ia[r], ib[r]] - hv
      H_full[ib[r], ia[r]] <- H_full[ib[r], ia[r]] - hv
    }
    g <- drop(crossprod(B, g_full))
    H <- crossprod(B, H_full %*% B)
    step <- solve(H, g)
    beta <- beta + drop(B %*% step)
    if (max(abs(step)) < tol) break
  }
  eta <- beta[ia] - beta[ib]
  P <- plogis(eta); V <- pmax(P * (1 - P), 1e-12)
  loglik <- sum(w * (y * log(pmax(P, 1e-12)) +
                     (1 - y) * log(pmax(1 - P, 1e-12))))
  converged <- max(abs(drop(crossprod(B, vapply(seq_len(K), function(k)
    sum((w * (y - P))[ia == k]) - sum((w * (y - P))[ib == k]), 0))))) < 1e-4

  # Godambe sandwich, clustered by judge when identified
  cl <- if (is.null(jd)) seq_along(ia) else jd
  res <- w * (y - P)
  Gm <- matrix(0, length(unique(cl)), K)
  rownames(Gm) <- as.character(unique(cl))
  for (r in seq_along(ia)) {
    rc <- as.character(cl[r])
    Gm[rc, ia[r]] <- Gm[rc, ia[r]] + res[r]
    Gm[rc, ib[r]] <- Gm[rc, ib[r]] - res[r]
  }
  Gb <- Gm %*% B
  J <- crossprod(Gb)
  Hi <- solve(crossprod(B, {
    H_full <- matrix(0, K, K)
    for (r in seq_along(ia)) {
      hv <- w[r] * V[r]
      H_full[ia[r], ia[r]] <- H_full[ia[r], ia[r]] + hv
      H_full[ib[r], ib[r]] <- H_full[ib[r], ib[r]] + hv
      H_full[ia[r], ib[r]] <- H_full[ia[r], ib[r]] - hv
      H_full[ib[r], ia[r]] <- H_full[ib[r], ia[r]] - hv
    }
    H_full
  } %*% B))
  covb <- Hi %*% J %*% Hi
  cov_beta <- B %*% covb %*% t(B)
  se <- sqrt(pmax(diag(cov_beta), 0))

  # fit: per-comparison z; objects and judges pool their cells
  z <- (y - P) / sqrt(V)
  c4v <- (1 - 3 * V) / V - 1                 # Bernoulli C4/V^2 - 1
  n_rows <- sum(w)
  f_cell <- (n_rows - (K - 1)) / n_rows
  pool <- function(sel) {
    if (sum(w[sel]) < 3)
      return(list(outfit_ms = NA_real_, fit_resid = NA_real_, df = NA_real_,
                  n = sum(w[sel])))
    y2 <- sum(w[sel] * z[sel]^2); f <- f_cell * sum(w[sel])
    v <- sum(w[sel] * c4v[sel])
    fr <- if (v > 1e-8 && y2 > 0) f * (log(y2) - log(f)) / sqrt(v) else NA_real_
    list(outfit_ms = y2 / f, fit_resid = fr, df = f, n = sum(w[sel]))
  }
  ofit <- lapply(seq_len(K), function(k) pool(ia == k | ib == k))
  objects <- data.frame(object = objs, location = beta, se = se,
                        comparisons = vapply(ofit, `[[`, 0, "n"),
                        wins = vapply(seq_len(K), function(k)
                          sum(w[ia == k] * y[ia == k]) +
                          sum(w[ib == k] * (1 - y[ib == k])), 0),
                        outfit_ms = vapply(ofit, `[[`, 0, "outfit_ms"),
                        fit_resid = vapply(ofit, `[[`, 0, "fit_resid"),
                        df_fit = vapply(ofit, `[[`, 0, "df"))
  rownames(objects) <- NULL

  judges <- NULL
  if (!is.null(jd)) {
    ju <- sort(unique(jd))
    jfit <- lapply(ju, function(j) pool(jd == j))
    judges <- data.frame(judge = ju,
                         n = vapply(jfit, `[[`, 0, "n"),
                         outfit_ms = vapply(jfit, `[[`, 0, "outfit_ms"),
                         fit_resid = vapply(jfit, `[[`, 0, "fit_resid"),
                         df_fit = vapply(jfit, `[[`, 0, "df"))
    rownames(judges) <- NULL
  }

  # classical pairwise goodness of fit: observed vs expected win proportions
  key <- ifelse(ia < ib, paste(ia, ib), paste(ib, ia))
  wins_lo <- tapply(w * ifelse(ia < ib, y, 1 - y), key, sum)
  n_pair <- tapply(w, key, sum)
  Ppair <- plogis(vapply(strsplit(names(n_pair), " "), function(s)
    beta[as.integer(s[1])] - beta[as.integer(s[2])], 0))
  zp <- (wins_lo - n_pair * Ppair) / sqrt(pmax(n_pair * Ppair * (1 - Ppair), 1e-12))
  idx <- do.call(rbind, strsplit(names(n_pair), " "))
  pairs <- data.frame(object_a = objs[as.integer(idx[, 1])],
                      object_b = objs[as.integer(idx[, 2])],
                      n = as.numeric(n_pair),
                      obs_prop = as.numeric(wins_lo / n_pair),
                      exp_prop = as.numeric(Ppair),
                      residual = as.numeric(zp),
                      chisq = as.numeric(zp^2))
  rownames(pairs) <- NULL
  used <- pairs$n >= 2
  total_chisq <- sum(pairs$chisq[used])
  total_df <- max(sum(used) - (K - 1), 1L)
  osi <- .psi(objects$location, objects$se)

  out <- list(objects = objects, pairs = pairs, judges = judges,
              total_chisq = total_chisq, total_df = total_df,
              total_p = pchisq(total_chisq, total_df, lower.tail = FALSE),
              osi = osi, loglik = loglik, iterations = it,
              converged = converged, n_comparisons = n_rows,
              clustered = !is.null(jd), cov_beta = cov_beta, notes = notes)
  class(out) <- "rmt_btl"
  out
}

#' @export
print.rmt_btl <- function(x, ...) {
  cat(sprintf("Bradley-Terry-Luce analysis: %d objects, %.0f comparisons%s\n",
              nrow(x$objects), x$n_comparisons,
              if (!is.null(x$judges)) sprintf(", %d judges", nrow(x$judges)) else ""))
  cat(sprintf("Conditional ML: %s in %d iterations; sandwich SEs%s\n",
              if (x$converged) "converged" else "NOT converged", x$iterations,
              if (x$clustered) " clustered by judge" else ""))
  cat(sprintf("Object separation index %.3f; pairwise chi-square %.2f on %d df, p = %s\n",
              x$osi$PSI, x$total_chisq, x$total_df, .fmt_p(x$total_p)))
  print(.fmt_df(x$objects[, intersect(c("object", "location", "se",
                                        "comparisons", "wins", "fit_resid"),
                                      names(x$objects))]), row.names = FALSE)
  if (!is.null(x$judges)) {
    mis <- x$judges[!is.na(x$judges$fit_resid) & abs(x$judges$fit_resid) > 2.5, ]
    cat(sprintf("Judges beyond |fit residual| 2.5: %d%s\n", nrow(mis),
                if (nrow(mis)) paste0(" (", paste(mis$judge, collapse = ", "), ")")
                else ""))
  }
  if (length(x$notes)) cat(sprintf("Notes: %s\n", paste(x$notes, collapse = "; ")))
  invisible(x)
}

#' Plot Bradley-Terry-Luce object locations
#'
#' Caterpillar plot of the object locations with 95 per cent error bars,
#' misfitting objects highlighted, in the package's house style.
#'
#' @param fit An object from \code{\link{btl}}.
#' @param band Absolute fit-residual value beyond which an object is
#'   highlighted.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pairs <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pairs[, 1], each = 30),
#'                 b = rep(pairs[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' plot_btl(btl(d, "a", "b", "win"))
#' @export
plot_btl <- function(fit, band = 2.5) {
  d <- fit$objects[order(fit$objects$location), ]
  k <- nrow(d)
  xlim <- range(c(d$location - 1.96 * d$se, d$location + 1.96 * d$se))
  op <- .rr_canvas(xlim + c(-0.15, 0.15) * diff(xlim), c(0.5, k + 0.5),
                   "Location (logits)", "", grid_y = FALSE, grid_x = TRUE,
                   yaxis = FALSE)
  on.exit(par(op))
  mis <- !is.na(d$fit_resid) & abs(d$fit_resid) > band
  segments(d$location - 1.96 * d$se, seq_len(k),
           d$location + 1.96 * d$se, seq_len(k),
           col = ifelse(mis, .rr$red, .rr$soft), lwd = 2.2)
  points(d$location, seq_len(k), pch = 21, cex = 1.6, lwd = 1.2,
         bg = ifelse(mis, .rr$red, .rr$blue), col = "white")
  text(d$location, seq_len(k), d$object, pos = 3, offset = 0.55, cex = 0.8,
       col = .rr$ink)
  if (any(mis))
    .rr_legend("bottomright", sprintf("|fit residual| > %.1f", band),
               pch = 21, pt.bg = .rr$red, col = "white", pt.cex = 1.4)
  invisible(NULL)
}
