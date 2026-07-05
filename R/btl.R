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
#'   (anything else is a tie). Ignored when \code{response} is given.
#' @param margin Optional name of a column holding the extent of the win
#'   ("a little", "much", ...), as an ordered factor or increasing values;
#'   combined with \code{winner} it assembles the graded response without
#'   any orientation bookkeeping ("B by much" means the same thing
#'   whichever column B sits in). Winner values matching neither object
#'   are ties and form the middle category.
#' @param thresholds \code{"free"} (default) estimates every symmetric
#'   threshold parameter; \code{"pc"} pools them to the spread (linear)
#'   principal component -- the symmetric case of the principal-component
#'   threshold structure, whose even skewness component is structurally
#'   zero here -- so thinly used categories borrow strength from every
#'   response. Both modes report the component decomposition.
#' @param response Optional name of a column holding a graded preference
#'   for \code{object_a} over \code{object_b} -- an ordered factor (worst
#'   to best for \code{object_a}) or integer scores \code{0..m}. Fits the
#'   adjacent-categories ordinal extension of BTL (Tutz 1986; Agresti
#'   1992): a partial-credit structure on the difference of locations with
#'   thresholds constrained symmetric, \code{tau_k = -tau_(m+1-k)}, so the
#'   model is invariant to presentation order. Two categories reproduce
#'   BTL exactly; three give the Davidson (1970) ties model.
#' @param judge Optional name of a judge column; enables the judge fit
#'   table and clusters the sandwich standard errors by judge.
#' @param count Optional name of a column of replication counts (a row
#'   standing for several identical comparisons).
#' @param ties How to treat ties in the dichotomous analysis:
#'   \code{"drop"} (default, removed with a note), \code{"half"} (half a
#'   win each way, a common pragmatic device -- flagged in the notes
#'   because the halves are not independent Bernoulli trials), or
#'   \code{"error"}. With graded responses, code ties as a middle
#'   category instead.
#' @param maxit,tol Newton-Raphson iteration cap and convergence tolerance.
#' @return A list of class \code{"rmt_btl"}: \code{objects} (location, se,
#'   comparisons, wins -- or the graded \code{score} -- outfit mean
#'   square, fit residual and its df),
#'   \code{pairs} (per pair: n, observed and expected win proportions --
#'   or mean graded responses --
#'   standardised residual, chi-square component), \code{judges} (when
#'   given: per judge n, outfit, fit residual, df), \code{total_chisq},
#'   \code{total_df}, \code{total_p}, the object separation index
#'   \code{osi}, \code{loglik}, convergence details, and \code{notes}.
#'   Graded fits add \code{thresholds} (the symmetric threshold estimates
#'   with standard errors), \code{m}, and \code{categories}.
#' @references Bradley, R. A. and Terry, M. E. (1952). Rank analysis of
#'   incomplete block designs: I. The method of paired comparisons.
#'   Biometrika, 39, 324-345. Luce, R. D. (1959). Individual Choice
#'   Behavior. Wiley. Andrich, D. (1978). Relationships between the
#'   Thurstone and Rasch approaches to item scaling. Applied Psychological
#'   Measurement, 2, 451-462.
#'
#'   Tutz, G. (1986). Bradley-Terry-Luce models with an ordered response.
#'   Journal of Mathematical Psychology, 30(3), 306-316. Agresti, A.
#'   (1992). Analysis of ordinal paired comparison data. Journal of the
#'   Royal Statistical Society C, 41(2), 287-297. Davidson, R. R. (1970).
#'   On extending the Bradley-Terry model to accommodate ties in paired
#'   comparison experiments. Journal of the American Statistical
#'   Association, 65(329), 317-328.
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
btl <- function(data, object_a, object_b, winner = NULL, response = NULL,
                margin = NULL, judge = NULL, count = NULL,
                ties = c("drop", "half", "error"),
                thresholds = c("free", "pc"), maxit = 60, tol = 1e-8) {
  ties <- match.arg(ties)
  thresholds <- match.arg(thresholds)
  data <- as.data.frame(data)
  if (is.null(winner) && is.null(response))
    stop("give either `winner` (dichotomous) or `response` (graded)")
  if (!is.null(margin) && is.null(winner))
    stop("`margin` requires `winner`")
  for (col in c(object_a, object_b, winner, response, margin, judge, count))
    if (!col %in% names(data)) stop("column not found: ", col)
  a <- trimws(as.character(data[[object_a]]))
  b <- trimws(as.character(data[[object_b]]))
  jd <- if (is.null(judge)) NULL else as.character(data[[judge]])
  w <- if (is.null(count)) rep(1, nrow(data)) else as.numeric(data[[count]])
  notes <- character(0)

  if (!is.null(response)) {
    xr <- data[[response]]
    if (is.factor(xr)) {
      cats <- levels(xr); x <- as.integer(xr) - 1L
    } else {
      x <- as.integer(round(as.numeric(xr)))
      if (any(x < 0, na.rm = TRUE))
        stop("graded responses must be non-negative integers 0..m")
      cats <- as.character(0:max(x, na.rm = TRUE))
    }
    keep <- !is.na(a) & !is.na(b) & !is.na(x) & a != b & !is.na(w) & w > 0
    if (any(!keep)) {
      notes <- c(notes, sprintf(
        "%d row(s) dropped (missing, zero-count, or self-comparison)",
        sum(!keep)))
      a <- a[keep]; b <- b[keep]; x <- x[keep]; w <- w[keep]
      if (!is.null(jd)) jd <- jd[keep]
    }
    if (!length(a)) stop("no usable comparisons")
    return(.btl_graded(a, b, x, jd, w, cats, maxit, tol, notes,
                       thr = thresholds))
  }

  if (!is.null(margin)) {
    # winner + margin entry: orientation-free by construction. The graded
    # response is assembled from "who won" and "by how much"; a winner value
    # matching neither object is a tie and becomes the middle category.
    mg <- data[[margin]]
    lv <- if (is.factor(mg)) levels(droplevels(mg)) else
      as.character(sort(unique(mg[!is.na(mg)])))
    q <- length(lv)
    if (q < 1L) stop("`margin` has no usable levels")
    mgi <- match(as.character(mg), lv)
    wn <- trimws(as.character(data[[winner]]))
    is_a <- !is.na(wn) & wn == a
    is_b <- !is.na(wn) & wn == b
    tie <- !is.na(wn) & !is_a & !is_b
    ties_present <- any(tie)
    keep <- !is.na(a) & !is.na(b) & !is.na(wn) & a != b & !is.na(w) & w > 0 &
      (tie | !is.na(mgi))
    if (any(!keep)) {
      notes <- c(notes, sprintf(
        "%d row(s) dropped (missing winner or margin, zero-count, or self-comparison)",
        sum(!keep)))
      a <- a[keep]; b <- b[keep]; w <- w[keep]
      mgi <- mgi[keep]; is_a <- is_a[keep]; is_b <- is_b[keep]
      tie <- tie[keep]
      if (!is.null(jd)) jd <- jd[keep]
    }
    if (!length(a)) stop("no usable comparisons")
    base <- q - 1L + as.integer(ties_present)
    x <- ifelse(is_a, base + mgi, ifelse(is_b, q - mgi, q))
    cats <- c(paste0("worse by ", rev(lv)), if (ties_present) "tie",
              paste0("better by ", lv))
    if (ties_present)
      notes <- c(notes, sprintf("%d tie(s) placed in the middle category",
                                sum(tie)))
    return(.btl_graded(a, b, as.integer(x), jd, w, cats, maxit, tol, notes,
                       thr = thresholds))
  }

  wn <- trimws(as.character(data[[winner]]))
  keep <- !is.na(a) & !is.na(b) & !is.na(wn) & a != b & !is.na(w) & w > 0
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
              clustered = !is.null(jd), cov_beta = cov_beta,
              comparisons = data.frame(object_a = a, object_b = b,
                                       response = y, weight = w),
              notes = notes)
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
  if (!is.null(x$thresholds)) {
    cat(sprintf("Graded comparisons in %d categories%s; symmetric thresholds: %s\n",
                x$m + 1L,
                if (!is.null(x$categories) &&
                    !all(x$categories == as.character(0:x$m)))
                  paste0(" (", paste(x$categories, collapse = " < "), ")")
                else "",
                paste(sprintf("%.3f", x$thresholds$tau), collapse = ", ")))
  }
  print(.fmt_df(x$objects[, intersect(c("object", "location", "se",
                                        "comparisons", "wins", "score",
                                        "fit_resid"),
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

# ---------------------------------------------------------------------------
# Graded paired comparisons: the adjacent-categories (Rasch-type) ordinal
# extension of BTL (Tutz 1986; Agresti 1992). The response is one of m+1
# ordered categories from object_a's perspective ("much worse" ... "much
# better"); category probabilities follow a partial-credit structure on the
# difference beta_a - beta_b with thresholds constrained symmetric,
# tau_k = -tau_{m+1-k}, so the model is invariant to presentation order and
# judge tendencies cancel. m = 1 is exactly BTL; m = 2 is the Davidson
# (1970) ties model. Estimation, identification, sandwich errors, and fit
# follow the package conventions established in btl().
# ---------------------------------------------------------------------------
.btl_graded <- function(a, b, x, jd, w, cats, maxit, tol, notes,
                        thr = "free") {
  m <- length(cats) - 1L
  if (m < 1L) stop("graded responses need at least two categories")
  # identifiability: empty EXTREME categories leave no finite spread (the
  # data are evidence of infinite spread, as a zero raw score is of an
  # infinite person location); empty interior categories are unidentified
  # under free thresholds but pooled over by the principal-component
  # structure
  xs <- c(x, m - x)
  emp <- which(tabulate(xs + 1L, m + 1L) == 0) - 1L
  if (any(emp %in% c(0L, m)))
    stop("extreme category never used (in either orientation): ",
         paste(cats[intersect(emp, c(0L, m)) + 1L], collapse = ", "),
         "; no finite threshold estimate exists - collapse categories")
  if (length(emp) && thr == "free")
    stop("interior category never used (in either orientation): ",
         paste(cats[emp + 1L], collapse = ", "),
         "; use thresholds = 'pc' (pooled principal-component structure)",
         " or collapse categories")
  if (length(emp))
    notes <- c(notes, sprintf(
      "interior category unused (%s); thresholds pooled by the principal-component structure",
      paste(cats[emp + 1L], collapse = ", ")))

  # objects whose every response sits at the boundary have no finite
  # estimate, as extreme persons are set aside in a Rasch calibration
  repeat {
    objs <- sort(unique(c(a, b)))
    T_of <- setNames(numeric(length(objs)), objs); N_of <- T_of
    for (r in seq_along(a)) {
      T_of[a[r]] <- T_of[a[r]] + w[r] * x[r]
      T_of[b[r]] <- T_of[b[r]] + w[r] * (m - x[r])
      N_of[a[r]] <- N_of[a[r]] + w[r] * m
      N_of[b[r]] <- N_of[b[r]] + w[r] * m
    }
    ext <- names(T_of)[T_of == 0 | T_of == N_of]
    if (!length(ext)) break
    notes <- c(notes, sprintf(
      "object(s) at a response boundary removed (no finite estimate): %s",
      paste(ext, collapse = ", ")))
    sel <- !(a %in% ext) & !(b %in% ext)
    a <- a[sel]; b <- b[sel]; x <- x[sel]; w <- w[sel]
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

  # symmetric-threshold map: tau = Cmat %*% tfree, tau_k = -tau_{m+1-k}.
  # Under thr = "pc" the free symmetric parameters are further pooled to
  # the spread (linear) component alone - Andrich's principal-component
  # structure with the even (skewness) components structurally zero under
  # symmetry - so sparse categories borrow strength from every response.
  if (thr == "pc" && m >= 2L) {
    v1 <- seq_len(m) - (m + 1) / 2
    Cmat <- cbind(v1 / sqrt(sum(v1^2)))
    q <- 1L
  } else {
    q <- m %/% 2L
    Cmat <- matrix(0, m, q)
    for (k in seq_len(q)) { Cmat[k, k] <- 1; Cmat[m + 1L - k, k] <- -1 }
  }
  Bmat <- rbind(diag(K - 1L), rep(-1, K - 1L))
  np <- (K - 1L) + q
  sc <- 0:m

  # per-row moments for current parameters: probabilities, E, V, mu4,
  # survivor S_k = P(X >= k) and EXc_k = E[X 1(X >= k)], k = 1..m
  moments <- function(beta, tfree) {
    tau <- if (q) drop(Cmat %*% tfree) else numeric(m)
    d <- beta[ia] - beta[ib]
    eta <- outer(d, sc) - matrix(rep(c(0, cumsum(tau)), each = length(d)),
                                 length(d), m + 1L)
    eta <- eta - apply(eta, 1, max)
    P <- exp(eta); P <- P / rowSums(P)
    E <- drop(P %*% sc)
    V <- drop(P %*% sc^2) - E^2
    mu4 <- drop(P %*% sc^4) - 4 * E * drop(P %*% sc^3) +
      6 * E^2 * drop(P %*% sc^2) - 3 * E^4
    S <- t(apply(P, 1, function(p) rev(cumsum(rev(p)))))[, -1L, drop = FALSE]
    EXc <- t(apply(P * rep(sc, each = length(d)), 1,
                   function(p) rev(cumsum(rev(p)))))[, -1L, drop = FALSE]
    list(P = P, E = E, V = pmax(V, 1e-12), mu4 = mu4, S = S, EXc = EXc,
         tau = tau)
  }
  cumInd <- outer(x, seq_len(m), ">=") * 1

  beta <- numeric(K); tfree <- numeric(q)
  for (it in seq_len(maxit)) {
    mo <- moments(beta, tfree)
    resE <- w * (x - mo$E)
    g_beta_full <- vapply(seq_len(K), function(k)
      sum(resE[ia == k]) - sum(resE[ib == k]), 0)
    g <- drop(crossprod(Bmat, g_beta_full))
    if (q) g <- c(g, drop(crossprod(Cmat, colSums(w * (mo$S - cumInd)))))
    H <- matrix(0, np, np)
    Hbb <- matrix(0, K, K)
    for (r in seq_along(ia)) {
      hv <- w[r] * mo$V[r]
      Hbb[ia[r], ia[r]] <- Hbb[ia[r], ia[r]] + hv
      Hbb[ib[r], ib[r]] <- Hbb[ib[r], ib[r]] + hv
      Hbb[ia[r], ib[r]] <- Hbb[ia[r], ib[r]] - hv
      Hbb[ib[r], ia[r]] <- Hbb[ib[r], ia[r]] - hv
    }
    H[1:(K - 1L), 1:(K - 1L)] <- crossprod(Bmat, Hbb %*% Bmat)
    if (q) {
      CovXc <- mo$EXc - mo$E * mo$S              # rows: Cov(X, 1(X>=k))
      Hbt_full <- matrix(0, K, q)
      wc <- (w * CovXc) %*% Cmat
      for (r in seq_along(ia)) {
        Hbt_full[ia[r], ] <- Hbt_full[ia[r], ] - wc[r, ]
        Hbt_full[ib[r], ] <- Hbt_full[ib[r], ] + wc[r, ]
      }
      H[1:(K - 1L), K:np] <- crossprod(Bmat, Hbt_full)
      H[K:np, 1:(K - 1L)] <- t(H[1:(K - 1L), K:np])
      Htt <- matrix(0, q, q)
      for (r in seq_along(ia)) {
        Sk <- mo$S[r, ]
        Covcc <- outer(seq_len(m), seq_len(m),
                       function(i, j) Sk[pmax(i, j)]) - tcrossprod(Sk)
        Htt <- Htt + w[r] * crossprod(Cmat, Covcc %*% Cmat)
      }
      H[K:np, K:np] <- Htt
    }
    step <- solve(H, g)
    beta <- beta + drop(Bmat %*% step[1:(K - 1L)])
    if (q) tfree <- tfree + step[K:np]
    if (max(abs(step)) < tol) break
  }
  mo <- moments(beta, tfree)
  tau <- mo$tau
  loglik <- sum(w * log(pmax(mo$P[cbind(seq_along(x), x + 1L)], 1e-300)))
  resE <- w * (x - mo$E)
  g_chk <- drop(crossprod(Bmat, vapply(seq_len(K), function(k)
    sum(resE[ia == k]) - sum(resE[ib == k]), 0)))
  if (q) g_chk <- c(g_chk, drop(crossprod(Cmat,
                                          colSums(w * (mo$S - cumInd)))))
  converged <- max(abs(g_chk)) < 1e-4

  # Godambe sandwich over the full parameter, clustered by judge
  cl <- if (is.null(jd)) as.character(seq_along(ia)) else jd
  ucl <- unique(cl)
  Gm <- matrix(0, length(ucl), np, dimnames = list(ucl, NULL))
  st_tau <- if (q) (w * (mo$S - cumInd)) %*% Cmat else NULL
  for (r in seq_along(ia)) {
    rc <- cl[r]
    zfull <- numeric(K)
    zfull[ia[r]] <- 1; zfull[ib[r]] <- -1
    Gm[rc, 1:(K - 1L)] <- Gm[rc, 1:(K - 1L)] +
      resE[r] * drop(crossprod(Bmat, zfull))
    if (q) Gm[rc, K:np] <- Gm[rc, K:np] + st_tau[r, ]
  }
  # H from the final iterate
  Hbb <- matrix(0, K, K)
  for (r in seq_along(ia)) {
    hv <- w[r] * mo$V[r]
    Hbb[ia[r], ia[r]] <- Hbb[ia[r], ia[r]] + hv
    Hbb[ib[r], ib[r]] <- Hbb[ib[r], ib[r]] + hv
    Hbb[ia[r], ib[r]] <- Hbb[ia[r], ib[r]] - hv
    Hbb[ib[r], ia[r]] <- Hbb[ib[r], ia[r]] - hv
  }
  H <- matrix(0, np, np)
  H[1:(K - 1L), 1:(K - 1L)] <- crossprod(Bmat, Hbb %*% Bmat)
  if (q) {
    CovXc <- mo$EXc - mo$E * mo$S
    Hbt_full <- matrix(0, K, q)
    wc <- (w * CovXc) %*% Cmat
    for (r in seq_along(ia)) {
      Hbt_full[ia[r], ] <- Hbt_full[ia[r], ] - wc[r, ]
      Hbt_full[ib[r], ] <- Hbt_full[ib[r], ] + wc[r, ]
    }
    H[1:(K - 1L), K:np] <- crossprod(Bmat, Hbt_full)
    H[K:np, 1:(K - 1L)] <- t(H[1:(K - 1L), K:np])
    Htt <- matrix(0, q, q)
    for (r in seq_along(ia)) {
      Sk <- mo$S[r, ]
      Covcc <- outer(seq_len(m), seq_len(m),
                     function(i, j) Sk[pmax(i, j)]) - tcrossprod(Sk)
      Htt <- Htt + w[r] * crossprod(Cmat, Covcc %*% Cmat)
    }
    H[K:np, K:np] <- Htt
  }
  Hi <- solve(H)
  covth <- Hi %*% crossprod(Gm) %*% Hi
  cov_beta <- Bmat %*% covth[1:(K - 1L), 1:(K - 1L), drop = FALSE] %*% t(Bmat)
  se <- sqrt(pmax(diag(cov_beta), 0))
  thresholds <- NULL; components <- NULL
  if (q) {
    cov_tau <- Cmat %*% covth[K:np, K:np, drop = FALSE] %*% t(Cmat)
    thresholds <- data.frame(threshold = seq_len(m), tau = tau,
                             se = sqrt(pmax(diag(cov_tau), 0)))
    # principal-component decomposition of the threshold structure: the
    # odd components (spread; kurtosis from five thresholds up) carry the
    # symmetric structure, the even skewness component is structurally
    # zero under presentation-order symmetry
    v1 <- seq_len(m) - (m + 1) / 2
    v1 <- v1 / sqrt(sum(v1^2))
    comp_rows <- list(data.frame(
      component = "spread", estimate = sum(v1 * tau),
      se = sqrt(pmax(drop(t(v1) %*% cov_tau %*% v1), 0))))
    if (m >= 4L) {
      v3 <- (seq_len(m) - (m + 1) / 2)^3
      v3 <- v3 - sum(v3 * v1) * v1
      v3 <- v3 / sqrt(sum(v3^2))
      comp_rows[[2]] <- data.frame(
        component = "kurtosis", estimate = sum(v3 * tau),
        se = if (thr == "pc") 0 else
          sqrt(pmax(drop(t(v3) %*% cov_tau %*% v3), 0)))
    }
    components <- do.call(rbind, comp_rows)
    rownames(components) <- NULL
  } else if (m > 1L) {
    thresholds <- data.frame(threshold = seq_len(m), tau = tau, se = 0)
  }

  # fit: per-comparison z; objects and judges pool their cells
  z <- (x - mo$E) / sqrt(mo$V)
  c4v <- mo$mu4 / mo$V^2 - 1
  n_rows <- sum(w)
  f_cell <- (n_rows - np) / n_rows
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
  score_of <- vapply(seq_len(K), function(k)
    sum(w[ia == k] * x[ia == k]) + sum(w[ib == k] * (m - x[ib == k])), 0)
  objects <- data.frame(object = objs, location = beta, se = se,
                        comparisons = vapply(ofit, `[[`, 0, "n"),
                        score = score_of,
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

  # pairwise goodness of fit on the oriented mean response
  key <- ifelse(ia < ib, paste(ia, ib), paste(ib, ia))
  x_lo <- ifelse(ia < ib, x, m - x)
  E_lo <- ifelse(ia < ib, mo$E, m - mo$E)
  n_pair <- tapply(w, key, sum)
  obs_m <- tapply(w * x_lo, key, sum) / n_pair
  exp_m <- tapply(w * E_lo, key, sum) / n_pair
  v_pair <- tapply(w * mo$V, key, sum)
  zp <- tapply(w * (x_lo - E_lo), key, sum) / sqrt(pmax(v_pair, 1e-12))
  idx <- do.call(rbind, strsplit(names(n_pair), " "))
  pairs <- data.frame(object_a = objs[as.integer(idx[, 1])],
                      object_b = objs[as.integer(idx[, 2])],
                      n = as.numeric(n_pair),
                      obs_mean = as.numeric(obs_m),
                      exp_mean = as.numeric(exp_m),
                      residual = as.numeric(zp),
                      chisq = as.numeric(zp)^2)
  rownames(pairs) <- NULL
  used <- pairs$n >= 2
  total_chisq <- sum(pairs$chisq[used])
  total_df <- max(sum(used) - (K - 1L) - q, 1L)
  osi <- .psi(objects$location, objects$se)

  out <- list(objects = objects, thresholds = thresholds,
              components = components, thr_structure = thr, pairs = pairs,
              judges = judges, m = m, categories = cats,
              total_chisq = total_chisq, total_df = total_df,
              total_p = pchisq(total_chisq, total_df, lower.tail = FALSE),
              osi = osi, loglik = loglik, iterations = it,
              converged = converged, n_comparisons = n_rows,
              clustered = !is.null(jd), cov_beta = cov_beta,
              comparisons = data.frame(object_a = a, object_b = b,
                                       response = x, weight = w),
              notes = notes)
  class(out) <- "rmt_btl"
  out
}

#' Plot graded-comparison category curves
#'
#' For a graded paired-comparison fit, the probability of each response
#' category as a function of the location difference
#' \code{beta_a - beta_b}, with the symmetric threshold structure marked.
#' The display is the paired-comparison counterpart of the category
#' probability curves of a polytomous item.
#'
#' @param fit A graded fit from \code{\link{btl}} (with \code{response}).
#' @param grid Difference grid, in logits.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 40), b = rep(pr[, 2], each = 40))
#' P <- vapply(seq_len(nrow(d)), function(r)
#'   item_moments(beta[d$a[r]] - beta[d$b[r]], c(-1, 0, 1))$P, numeric(4))
#' d$grade <- apply(P, 2, function(p) sample(0:3, 1, prob = p))
#' plot_btl_categories(btl(d, "a", "b", response = "grade"))
#' @export
plot_btl_categories <- function(fit, grid = seq(-4, 4, 0.05)) {
  if (is.null(fit$m) || fit$m < 2L)
    stop("category curves need a graded fit (three or more categories)")
  tau <- fit$thresholds$tau
  P <- vapply(grid, function(d) item_moments(d, tau)$P, numeric(fit$m + 1L))
  op <- .rr_canvas(range(grid), c(0, 1),
                   "Location difference (logits)", "Category probability")
  on.exit(par(op))
  abline(v = tau, lty = 3, col = .rr$soft)
  labs <- if (!is.null(fit$categories)) fit$categories else
    as.character(0:fit$m)
  for (cat in 0:fit$m)
    lines(grid, P[cat + 1L, ], lwd = 2.6,
          col = .rr$pal[cat %% length(.rr$pal) + 1L])
  .rr_legend("right", labs, lwd = 2.6,
             col = .rr$pal[(0:fit$m) %% length(.rr$pal) + 1L])
  invisible(NULL)
}

#' Plot an object characteristic curve
#'
#' The paired-comparison counterpart of the item characteristic curve: the
#' model expected response for one object as a function of opponent
#' location (the win probability, or the expected graded response), with
#' the observed mean response against each opponent overlaid at that
#' opponent\'s estimated location. Observed points shrink in toward the
#' curve as the model holds; an object of inconsistent quality shows
#' points straying from it, exactly as a misfitting item does.
#'
#' @param fit An object from \code{\link{btl}}.
#' @param object Object name.
#' @param grid Opponent-location grid, in logits.
#' @param min_n Opponents met fewer than this many times are drawn hollow.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
#' pr <- t(combn(names(beta), 2))
#' d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
#' p <- plogis(beta[d$a] - beta[d$b])
#' d$win <- ifelse(runif(nrow(d)) < p, d$a, d$b)
#' plot_btl_icc(btl(d, "a", "b", winner = "win"), "C")
#' @export
plot_btl_icc <- function(fit, object, grid = NULL, min_n = 10) {
  ob <- fit$objects
  if (!object %in% ob$object) stop("no such object: ", object)
  m <- if (is.null(fit$m)) 1L else fit$m
  tau <- if (!is.null(fit$thresholds)) fit$thresholds$tau else numeric(1)
  b_o <- ob$location[ob$object == object]
  if (is.null(grid)) {
    rng <- range(ob$location) + c(-1, 1)
    grid <- seq(rng[1], rng[2], length.out = 201)
  }
  Ecurve <- vapply(grid, function(t) item_moments(b_o - t, tau)$E, 0)
  cm <- fit$comparisons
  sel_a <- cm$object_a == object
  sel_b <- cm$object_b == object
  opp <- c(cm$object_b[sel_a], cm$object_a[sel_b])
  resp <- c(cm$response[sel_a], m - cm$response[sel_b])
  wt <- c(cm$weight[sel_a], cm$weight[sel_b])
  keep <- opp %in% ob$object
  opp <- opp[keep]; resp <- resp[keep]; wt <- wt[keep]
  obs <- data.frame(
    opponent = tapply(opp, opp, `[`, 1),
    loc = ob$location[match(names(tapply(wt, opp, sum)), ob$object)],
    mean = as.numeric(tapply(wt * resp, opp, sum) / tapply(wt, opp, sum)),
    n = as.numeric(tapply(wt, opp, sum)))
  op <- .rr_canvas(range(grid), c(0, m), "Opponent location (logits)",
                   if (m == 1L) "Probability preferred" else
                     "Expected graded response",
                   sprintf("%s  (location %.3f)", object, b_o))
  on.exit(par(op))
  lines(grid, Ecurve, lwd = 3, col = .rr$ink)
  abline(v = b_o, lty = 3, col = .rr$soft)
  solid <- obs$n >= min_n
  points(obs$loc[solid], obs$mean[solid], pch = 21, bg = .rr$blue,
         col = "white", cex = 1.5, lwd = 1.2)
  if (any(!solid))
    points(obs$loc[!solid], obs$mean[!solid], pch = 21, bg = "white",
           col = .rr$blue, cex = 1.3, lwd = 1.4)
  text(obs$loc, obs$mean, obs$opponent, pos = 3, offset = 0.45,
       cex = 0.72, col = .rr$soft)
  .rr_legend("topright",
             c("Model", "Observed (per opponent)",
               if (any(!solid)) sprintf("fewer than %d comparisons", min_n)),
             lwd = c(3, NA, if (any(!solid)) NA),
             pch = c(NA, 21, if (any(!solid)) 21),
             pt.bg = c(NA, .rr$blue, if (any(!solid)) "white"),
             col = c(.rr$ink, "white", if (any(!solid)) .rr$blue),
             pt.cex = 1.3)
  invisible(NULL)
}
