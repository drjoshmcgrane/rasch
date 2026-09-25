# rasch :: plot-cj
# ===========================================================================
# The joint calibration map: each item's combined location from responses and
# judgements with its interval, and beside it the location each frame gives
# the item on its own, so the agreement the invariance test summarises can be
# read item by item.
# ===========================================================================

#' Plot a joint calibration of responses and judgements
#'
#' Caterpillar plot of the item locations from \code{\link{rasch_cj}}: the
#' combined location of each item with its 95 per cent interval, and beside
#' it the location each frame (the responses, the comparisons, the rankings,
#' or each test when several are linked) gives the item on its own, expressed
#' on the reference scale. The spread of the frame markers around the
#' combined location is the evidence the invariance test summarises; an item
#' whose location differs between frames at a Holm-adjusted p below .05 in
#' the per-object invariance table is drawn in red. A frame that does not
#' reach an item leaves no marker for it.
#'
#' @param fit An item-mode object from \code{\link{rasch_cj}}.
#' @param frames Logical; draw the separate calibration of each frame beside
#'   the combined location (the default), or the combined locations alone.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' delta <- seq(-1.5, 1.5, length.out = 8)
#' names(delta) <- sprintf("I%02d", 1:8)
#' theta <- rnorm(300)
#' X <- sapply(delta, function(d) as.integer(runif(300) < plogis(theta - d)))
#' pairs <- t(combn(names(delta), 2))[sample(28, 200, replace = TRUE), ]
#' p_a <- plogis(0.6 * (delta[pairs[, 1]] - delta[pairs[, 2]]))
#' cj <- data.frame(a = pairs[, 1], b = pairs[, 2],
#'                  winner = ifelse(runif(200) < p_a, pairs[, 1], pairs[, 2]))
#' fit <- rasch_cj(X, comparisons = cj, object_a = "a", object_b = "b",
#'                 winner = "winner")
#' plot_cj(fit)
#' @seealso \code{\link{rasch_cj}}, \code{\link{plot_pl}}.
#' @export
plot_cj <- function(fit, frames = TRUE) {
  if (!inherits(fit, "rasch_cj"))
    stop("not a joint calibration (rasch_cj) fit", call. = FALSE)
  if (identical(fit$mode, "persons"))
    stop("plot_cj() draws the item calibration; a person-mode rasch_cj fit ",
         "measures the persons and has no item locations to draw",
         call. = FALSE)
  if (!isTRUE(fit$converged))
    stop("the joint calibration did not converge; fitted displays are unavailable",
         call. = FALSE)
  if (!is.logical(frames) || length(frames) != 1L || is.na(frames))
    stop("`frames` must be TRUE or FALSE", call. = FALSE)
  d <- fit$items[order(fit$items$location), ]
  k <- nrow(d)
  sep <- if (frames) grep("^location_", names(d), value = TRUE) else character(0)
  frame_names <- sub("^location_", "", sep)
  # An item moves when any judgement frame places it away from its reference
  # location; the per-object table carries one row per frame and object.
  tab <- fit$invariance$items
  moved <- if (is.null(tab)) rep(FALSE, k) else
    d$item %in% tab$item[!is.na(tab$p_adj) & tab$p_adj < 0.05]
  has_se <- is.finite(d$se) & d$se > 0
  lower <- upper <- rep(NA_real_, k)
  lower[has_se] <- d$location[has_se] - 1.96 * d$se[has_se]
  upper[has_se] <- d$location[has_se] + 1.96 * d$se[has_se]
  xall <- c(lower, upper, d$location, unlist(d[sep], use.names = FALSE))
  xlim <- range(xall[is.finite(xall)])
  if (!diff(xlim)) xlim <- xlim + c(-1, 1)
  op <- .rr_canvas(xlim + c(-0.15, 0.15) * diff(xlim), c(0.5, k + 0.5),
                   "Location (logits)", "", grid_y = FALSE, grid_x = TRUE,
                   yaxis = FALSE)
  on.exit(par(op))
  y <- seq_len(k)
  segments(lower[has_se], y[has_se], upper[has_se], y[has_se],
           col = ifelse(moved[has_se], .rr$red, .rr$soft), lwd = 2.2)
  points(d$location, y, pch = 21, cex = 1.3, lwd = 1.2,
         bg = ifelse(moved, .rr$red, .rr$blue), col = "white")
  # One open marker shape and colour per frame, drawn over the combined
  # point so a frame that agrees with it shows as a ring around the dot;
  # the colours keep clear of the combined blue and the red of a moving item.
  frame_pch <- c(22, 24, 23, 25)[(seq_along(sep) - 1L) %% 4L + 1L]
  frame_col <- .rr$pal[c(4, 3, 5, 6, 7, 8)][(seq_along(sep) - 1L) %% 6L + 1L]
  for (j in seq_along(sep))
    points(d[[sep[j]]], y, pch = frame_pch[j], cex = 1.7, lwd = 1.5,
           col = frame_col[j])
  text(d$location, y, d$item, pos = 3, offset = 0.6, cex = 0.8, col = .rr$ink)
  labs <- c("Combined",
            sprintf("%s alone", paste0(toupper(substr(frame_names, 1, 1)),
                                       substring(frame_names, 2))),
            if (any(moved)) "Differs between frames (Holm p < .05)")
  .rr_legend("bottomright", labs,
             pch = c(21, frame_pch, if (any(moved)) 21),
             pt.bg = c(.rr$blue, rep(NA, length(sep)), if (any(moved)) .rr$red),
             col = c("white", frame_col, if (any(moved)) "white"),
             pt.lwd = 1.5,
             pt.cex = c(1.3, rep(1.4, length(sep)), if (any(moved)) 1.3))
  invisible(NULL)
}
