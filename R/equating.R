# rasch :: test equating
# ===========================================================================
# Comparison of common-item locations across two separately analysed
# datasets. Because each analysis fixes its own origin, the comparison
# allows for a scale shift (estimated by precision-weighted mean difference)
# and then tests each common item against the shifted identity line:
#
#   t_i = (b1_i - b2_i - shift) / sqrt(se1_i^2 + se2_i^2)
#
# Items that survive define the equating link; items that fail show item
# drift and should be dropped from the link (or anchored individually).
# ===========================================================================

# Item-location covariance over a set of items: block means of the fit's
# threshold covariance (recentred parameterisation), which carries the
# negative correlations induced by the identification constraint. A bank
# table has no covariance and contributes diag(se^2) -- conservative.
.equate_loc_cov <- function(obj, items) {
  if (inherits(obj, "rasch") && !is.null(obj$est$cov_tau)) {
    thr <- obj$est$thr
    # thr$item holds integer item positions in column order, which is the
    # order of obj$items: the match below is by position, not name
    idx <- match(items, obj$items$item)
    rows <- lapply(idx, function(i) thr$id[thr$item == i])
    S <- matrix(0, length(items), length(items))
    for (i in seq_along(rows)) for (j in seq_along(rows))
      S[i, j] <- mean(obj$est$cov_tau[rows[[i]], rows[[j]]])
    return(S)
  }
  ref <- .equate_ref(obj)
  se <- ref$se[match(items, ref$item)]
  diag(se^2, length(items))
}

.equate_ref <- function(reference) {
  if (inherits(reference, "rasch"))
    return(data.frame(item = reference$items$item,
                      location = reference$items$location,
                      se = reference$items$se))
  reference <- as.data.frame(reference)
  if (!all(c("item", "location") %in% names(reference)))
    stop("reference needs columns item, location (and ideally se)")
  if (is.null(reference$se)) reference$se <- 0
  reference[, c("item", "location", "se")]
}

#' Equate two test calibrations through their common items
#'
#' Compares the item locations of a fit with those of a second fit (or a
#' reference table such as an item bank), matched by item name. A scale
#' shift between the two origins is estimated by the precision-weighted mean
#' difference, and each common item is then tested against the shifted
#' identity line; flagged items show drift and weaken the equating link.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param reference A second \code{\link{rasch}} fit, or a data frame with
#'   columns \code{item}, \code{location}, and optionally \code{se}.
#' @param shift \code{"mean"} (default) allows a scale shift between the two
#'   analyses; \code{"none"} compares raw locations, appropriate when both
#'   analyses are already on a shared (anchored) scale.
#' @return A list with the comparison \code{table} (locations, standard
#'   errors, difference, t, raw and BH-adjusted p, drift flag), the
#'   estimated \code{shift},
#'   the location \code{correlation}, the root mean square difference after
#'   shifting (\code{rmsd}), and the number of common items \code{n}.
#' @examples
#' set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
#' mk <- function() {
#'   X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#'   colnames(X) <- paste0("I", 1:8); rasch(X)
#' }
#' eq <- equate_tests(mk(), mk())
#' eq$table
#' @export
equate_tests <- function(fit, reference, shift = c("mean", "none")) {
  shift <- match.arg(shift)
  ref <- .equate_ref(reference)
  cur <- data.frame(item = fit$items$item, location = fit$items$location,
                    se = fit$items$se)
  common <- intersect(cur$item, ref$item)
  if (length(common) < 2) stop("need at least two common items to equate")
  a <- cur[match(common, cur$item), ]
  b <- ref[match(common, ref$item), ]
  d <- a$location - b$location
  v <- a$se^2 + b$se^2
  # an item whose location or SE is unavailable (for example a weakly
  # determined item whose SE is honestly NA) cannot contribute to the
  # precision-weighted shift or the drift tests, but it must not poison
  # the remaining items either: exclude it, say so, and carry its row in
  # the table with NA test columns
  usable <- is.finite(d) & is.finite(v)
  note <- NULL
  if (any(!usable)) {
    note <- sprintf(
      "common item(s) excluded from the shift and drift tests (location or SE unavailable): %s",
      paste(common[!usable], collapse = ", "))
    if (sum(usable) < 2)
      stop("fewer than two common items with usable locations and standard ",
           "errors remain; ", note)
  }
  # the shift c0 is estimated from the same common items the drift tests
  # then examine, and each calibration's locations are correlated through
  # its identification constraint: the drift denominators use
  # Var(d_i - c0) = [(I - 1u') Sigma (I - u 1')]_ii over the usable items,
  # with Sigma the sum of the two calibrations' item-location covariances
  Sg <- .equate_loc_cov(fit, common) + .equate_loc_cov(reference, common)
  var_d <- rep(NA_real_, length(common))
  if (shift == "mean") {
    w <- 1 / pmax(v[usable], 1e-10)
    c0 <- sum(w * d[usable]) / sum(w)
    u <- w / sum(w)
    Suu <- Sg[usable, usable, drop = FALSE]
    Su <- drop(Suu %*% u)
    var_d[usable] <- pmax(diag(Suu) - 2 * Su + drop(t(u) %*% Su), 1e-10)
  } else {
    c0 <- 0
    var_d[usable] <- pmax(diag(Sg)[usable], 1e-10)
  }
  t <- ifelse(usable, (d - c0) / sqrt(var_d), NA_real_)
  p <- 2 * pnorm(-abs(t))
  n <- sum(usable)
  p_adj <- rep(NA_real_, length(p))
  p_adj[usable] <- p.adjust(p[usable], method = "BH")
  tab <- data.frame(item = common,
                    location_1 = a$location, se_1 = a$se,
                    location_2 = b$location, se_2 = b$se,
                    difference = d, adj_difference = d - c0,
                    t = t, p = p, p_adj = p_adj,
                    drift = p_adj < 0.05)
  rownames(tab) <- NULL
  structure(class = "rasch_equate", list(table = tab, shift = c0,
       correlation = cor(a$location[usable], b$location[usable]),
       rmsd = sqrt(mean((d[usable] - c0)^2)), n = n, note = note))
}

#' Plot a test-equating comparison
#'
#' Scatter of the two calibrations' common-item locations with the shifted
#' identity line and per-item 95 per cent bands; drifting items
#' (BH-adjusted) are highlighted and labelled.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param reference A second \code{\link{rasch}} fit, or a data frame with
#'   columns \code{item}, \code{location}, and optionally \code{se}.
#' @param shift Passed to \code{\link{equate_tests}}.
#' @return Called for its plotting side effect; invisibly the
#'   \code{\link{equate_tests}} result.
#' @examples
#' set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
#' mk <- function() {
#'   X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
#'   colnames(X) <- paste0("I", 1:8); rasch(X)
#' }
#' plot_equate(mk(), mk())
#' @export
plot_equate <- function(fit, reference, shift = c("mean", "none")) {
  eq <- equate_tests(fit, reference, shift)
  tab <- eq$table
  rng <- range(c(tab$location_1, tab$location_2)) + c(-0.4, 0.4)
  op <- .rr_canvas(rng, rng, "Reference location (logits)",
                   "Current location (logits)",
                   sprintf("%d common items, shift %.3f, r = %.3f",
                           eq$n, eq$shift, eq$correlation), grid_x = TRUE)
  on.exit(par(op))
  abline(eq$shift, 1, col = .rr$ink, lwd = 2)
  # excluded items (NA SEs) still appear as points, but the average band
  # and the drift highlighting are computed over the usable rows only
  band <- 1.96 * sqrt(mean(tab$se_1^2 + tab$se_2^2, na.rm = TRUE))
  abline(eq$shift + band, 1, lty = 3, col = .rr$soft)
  abline(eq$shift - band, 1, lty = 3, col = .rr$soft)
  segments(tab$location_2, tab$location_1 - 1.96 * tab$se_1,
           tab$location_2, tab$location_1 + 1.96 * tab$se_1,
           col = paste0(.rr$soft, "88"))
  dr <- tab$drift %in% TRUE
  points(tab$location_2, tab$location_1, pch = 21, cex = 1.6,
         bg = ifelse(dr, .rr$red, .rr$blue), col = "white", lwd = 1.2)
  if (any(dr))
    text(tab$location_2[dr], tab$location_1[dr],
         tab$item[dr], pos = 3, offset = 0.5, cex = 0.75, col = .rr$red)
  invisible(eq)
}


#' @export
print.rasch_equate <- function(x, ...) {
  cat(sprintf("Common-item equating over %d item(s): shift %.3f, correlation %.3f, RMSD %.3f\n",
              x$n, x$shift, x$correlation, x$rmsd))
  core <- c("item", "location_1", "location_2", "adj_difference", "t",
            "p_adj", "drift")
  print(.fmt_df(x$table[, intersect(core, names(x$table))]), row.names = FALSE)
  cat("(standard errors and unadjusted columns on $table)\n")
  if (!is.null(x$note)) cat("Note:", x$note, "\n")
  invisible(x)
}
