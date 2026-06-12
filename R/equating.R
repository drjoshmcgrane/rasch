# RaschR :: test equating
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
  w <- 1 / pmax(v, 1e-10)
  c0 <- if (shift == "mean") sum(w * d) / sum(w) else 0
  t <- (d - c0) / sqrt(pmax(v, 1e-10))
  p <- 2 * pnorm(-abs(t))
  n <- length(common)
  p_adj <- p.adjust(p, method = "BH")
  tab <- data.frame(item = common,
                    location_1 = a$location, se_1 = a$se,
                    location_2 = b$location, se_2 = b$se,
                    difference = d, adj_difference = d - c0,
                    t = t, p = p, p_adj = p_adj,
                    drift = p_adj < 0.05)
  rownames(tab) <- NULL
  list(table = tab, shift = c0,
       correlation = cor(a$location, b$location),
       rmsd = sqrt(mean((d - c0)^2)), n = n)
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
                   sprintf("Test equating: %d common items, shift %.3f, r = %.3f",
                           eq$n, eq$shift, eq$correlation), grid_x = TRUE)
  on.exit(par(op))
  abline(eq$shift, 1, col = .rr$ink, lwd = 2)
  abline(eq$shift + 1.96 * sqrt(mean(tab$se_1^2 + tab$se_2^2)), 1,
         lty = 3, col = .rr$soft)
  abline(eq$shift - 1.96 * sqrt(mean(tab$se_1^2 + tab$se_2^2)), 1,
         lty = 3, col = .rr$soft)
  segments(tab$location_2, tab$location_1 - 1.96 * tab$se_1,
           tab$location_2, tab$location_1 + 1.96 * tab$se_1,
           col = paste0(.rr$soft, "88"))
  points(tab$location_2, tab$location_1, pch = 21, cex = 1.6,
         bg = ifelse(tab$drift, .rr$red, .rr$blue), col = "white", lwd = 1.2)
  if (any(tab$drift))
    text(tab$location_2[tab$drift], tab$location_1[tab$drift],
         tab$item[tab$drift], pos = 3, offset = 0.5, cex = 0.75, col = .rr$red)
  invisible(eq)
}
