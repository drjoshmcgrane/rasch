# rasch :: plots
# ===========================================================================
# The Rasch diagnostic plot suite in base graphics with a modern flat style:
# item characteristic curves (with group overlay for DIF), category
# probability curves, threshold probability curves, the person-item
# threshold distribution, the threshold map, the test characteristic curve,
# test information and measurement error, item and person fit maps, the
# residual-correlation heatmap, residual principal-component loadings, and
# category frequencies.
# ===========================================================================

.rr <- list(
  ink    = "#1e293b", grid = "#e2e8f0", soft = "#94a3b8",
  blue   = "#2563eb", red  = "#dc2626", amber = "#f59e0b",
  teal   = "#0f766e", purple = "#7c3aed",
  pal    = c("#2563eb", "#dc2626", "#f59e0b", "#0f766e", "#7c3aed",
             "#db2777", "#65a30d", "#475569")
)

# Modern flat canvas: light horizontal grid, open axes. A title is drawn
# only when `main` carries information (item, person, summary figures);
# plot-type names are left to the surrounding context.
.rr_canvas <- function(xlim, ylim, xlab, ylab, main = "", grid_y = TRUE,
                       grid_x = FALSE, yaxis = TRUE, right = 1.5,
                       xaxis = TRUE) {
  has_main <- !is.null(main) && nzchar(main)
  op <- par(mar = c(4.2, 4.4, if (has_main) 3.2 else 1.6, right),
            mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15, cex.lab = 1.0)
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,
       main = "", axes = FALSE)
  if (has_main) title(main = main, adj = 0, line = 1.4)
  xt <- .rr_ticks(xlim); yt <- .rr_ticks(ylim)
  if (grid_y) abline(h = yt, col = .rr$grid, lwd = 0.8)
  if (grid_x) abline(v = xt, col = .rr$grid, lwd = 0.8)
  if (xaxis) axis(1, at = xt, col = .rr$grid, col.ticks = .rr$soft)
  if (yaxis) axis(2, at = yt, col = .rr$grid, col.ticks = .rr$soft)
  invisible(op)
}


# A drawing grid is a vector of at least two finite locations.
# Display limits are a pair of finite ascending values; a map whose limits
# admit nothing is not a display, so an empty range is reported rather than
# drawn as empty panels or left to emit internal warnings.
.check_xlim <- function(xlim) {
  if (!is.null(xlim) &&
      (length(xlim) != 2L || !is.numeric(xlim) || is.complex(xlim) ||
       !is.null(dim(xlim)) || !is.null(oldClass(xlim)) ||
       any(!is.finite(xlim)) ||
       xlim[1] >= xlim[2]))
    stop("`xlim` must be two finite ascending limits", call. = FALSE)
  invisible(xlim)
}

# A repeated person identifier -- a stacked or racked longitudinal design --
# does not select one row. Taking the first would draw a different occasion
# from the one the caller asked for, with nothing on the display to say so.
.person_row <- function(fit, person) {
  if (missing(person) || !is.atomic(person) || !is.null(dim(person)) ||
      length(person) != 1L || is.na(person))
    stop("`person` must name or index exactly one person", call. = FALSE)
  if (is.numeric(person) && person %in% seq_len(nrow(fit$X)))
    return(as.integer(person))
  hits <- which(as.character(fit$person$id) == as.character(person))
  if (!length(hits)) stop("person not found", call. = FALSE)
  if (length(hits) > 1L)
    stop("person identifier '", person, "' appears in ", length(hits),
         " rows (", paste(utils::head(hits, 5), collapse = ", "),
         if (length(hits) > 5) ", ..." else "",
         "); give the row number instead", call. = FALSE)
  hits
}

.check_cex <- function(x, name = "cex_labels") {
  if (length(x) != 1L || !is.numeric(x) || is.complex(x) ||
      !is.null(dim(x)) || !is.null(oldClass(x)) || !is.finite(x) || x <= 0)
    stop("`", name, "` must be one positive finite size", call. = FALSE)
  invisible(x)
}

.check_map_range <- function(n_thresholds, n_persons = NA_integer_) {
  if (n_thresholds == 0L)
    stop("the display range contains no thresholds; widen `xlim`",
         call. = FALSE)
  if (!is.na(n_persons) && n_persons == 0L)
    stop("the display range contains no person estimates; widen `xlim`",
         call. = FALSE)
  invisible(TRUE)
}

# Targeting displays use fitted person and threshold locations directly.  A
# failed calibration can still contain the last numerical iterate, which is
# drawable but is not an estimate and must not be presented as one.
.check_response_display_fit <- function(fit, display = "fitted displays") {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be a fitted response-data Rasch model", call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; ", display, " are unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; ", display, " are unavailable",
         call. = FALSE)
  invisible(fit)
}

.check_targeting_fit <- function(fit) {
  .check_response_display_fit(fit, "targeting plots")
}

# The observed points must sit on the SAME class intervals as the item-trait
# test they illustrate. With missing data the fit carries a per-item
# allocation (each item's own responders, its own interval count), so the
# display reads that allocation rather than recomputing a whole-sample one --
# unless the caller asked for a specific number of intervals, or the points
# are split by a person factor, which is a different allocation by design.
.fit_class_intervals <- function(fit, i, x, th, ex, n_groups,
                                 n_groups_given = FALSE, group = NULL) {
  if (!n_groups_given && is.null(group) && !is.null(fit$ci_item) &&
      length(fit$ci_item) >= i && !is.null(fit$ci_item[[i]]))
    return(fit$ci_item[[i]])
  .class_intervals(ifelse(is.na(x), NA_real_, th), ex, n_groups)
}

.check_grid <- function(grid) {
  if (!is.numeric(grid) || is.complex(grid) || !is.null(dim(grid)) ||
      !is.null(oldClass(grid)) || length(grid) < 2L ||
      any(!is.finite(grid)) ||
      any(diff(grid) <= 0))
    stop("`grid` must be a strictly increasing vector of at least two finite locations",
         call. = FALSE)
  invisible(grid)
}

# Default curve grids follow the fitted origin. The location origin is
# arbitrary and may be set by external anchors; a fixed grid around zero can
# therefore draw an empty curve for a valid fit.
.default_model_grid <- function(fit, half_width = 6, by = 0.05) {
  tau <- unlist(fit$tau_list, use.names = FALSE)
  tau <- tau[is.finite(tau)]
  if (!length(tau))
    stop("the fitted thresholds are unavailable; model curves cannot be drawn",
         call. = FALSE)
  origin <- mean(tau)
  span <- max(half_width, max(abs(tau - origin)) + 1)
  n <- max(2L, ceiling(2 * span / by) + 1L)
  seq(origin - span, origin + span, length.out = n)
}

.model_grid <- function(fit, grid = NULL, half_width = 6, by = 0.05) {
  if (!inherits(fit, "rasch") || inherits(fit, "rasch_btl"))
    stop("`fit` must be a fitted response-data Rasch model", call. = FALSE)
  if (!isTRUE(fit$est$converged))
    stop("the fitted calibration did not converge; model curves are unavailable",
         call. = FALSE)
  if (!.efrm_link_converged(fit))
    stop("the fitted set-unit link did not converge; model curves are unavailable",
         call. = FALSE)
  if (is.null(grid))
    grid <- .default_model_grid(fit, half_width, by)
  .check_grid(grid)
  grid
}

# A misfit band is a positive finite half-width; a logical flag must be
# stated as TRUE or FALSE, not NA.
.check_band <- function(band) {
  if (length(band) != 1L || !is.numeric(band) || is.complex(band) ||
      !is.null(dim(band)) || !is.null(oldClass(band)) ||
      !is.finite(band) || band <= 0)
    stop("`band` must be one positive finite half-width", call. = FALSE)
  invisible(band)
}
.check_flag <- function(x, name) {
  if (!isTRUE(x) && !isFALSE(x))
    stop("`", name, "` must be TRUE or FALSE", call. = FALSE)
  invisible(x)
}

# Whole-unit tick marks whenever the span allows them, so a logit axis is
# labelled at every logit and its extremes are not left between labels. One
# extra integer on each side reaches into the axis expansion; axis() clips
# marks outside the plot region.
.rr_ticks <- function(lim) {
  w <- max(lim) - min(lim)
  if (w > 2.5 && w <= 14)
    seq(ceiling(min(lim) - 1e-9) - 1, floor(max(lim) + 1e-9) + 1)
  else pretty(lim, n = 7)
}

# House axis for plots drawn outside .rr_canvas: ticks from the realised
# plot window, so every location axis follows the same whole-unit rule.
.rr_axis <- function(side, ...) {
  usr <- par("usr")
  lim <- if (side %in% c(1, 3)) usr[1:2] else usr[3:4]
  axis(side, at = .rr_ticks(lim), col = .rr$grid, col.ticks = .rr$soft, ...)
}

.rr_legend <- function(pos, ..., cex = 0.85)
  legend(pos, ..., bty = "n", text.col = .rr$ink, cex = cex)

.item_idx <- function(fit, item) {
  if (!is.atomic(item) || !is.null(dim(item)))
    stop("item must be an ordinary vector of names or indices",
         call. = FALSE)
  if (!length(item)) stop("item must name or index at least one item")
  if (is.character(item) || is.factor(item)) {
    supplied <- as.character(item)
    if (anyNA(supplied) || any(!nzchar(trimws(supplied))))
      stop("item names must be non-missing and non-empty", call. = FALSE)
    # Keep deliberately spaced column names addressable. Otherwise surrounding
    # whitespace is selector syntax rather than part of an ordinary item name.
    i <- match(supplied, fit$items$item)
    fallback <- is.na(i)
    if (any(fallback))
      i[fallback] <- match(.role_text_values(supplied[fallback]),
                           fit$items$item)
    if (anyNA(i))
      stop("no such item: ",
           paste(supplied[is.na(i)], collapse = ", "), call. = FALSE)
    i
  } else {
    if (!is.numeric(item) || any(!is.finite(item)) ||
        any(item != floor(item)) || any(item < 1) ||
        any(item > nrow(fit$items)))
      stop("item index must be whole numbers between 1 and ",
           nrow(fit$items), call. = FALSE)
    as.integer(item)
  }
}

# Discrimination (frame unit) of column i; 1 unless the fit carries units.
.disc_of <- function(fit, i) if (is.null(fit$disc)) 1 else fit$disc[i]

# A pooled MFRM DIF row names the item, whereas fit$items contains the
# item-by-facet response cells. Put every observed cell back on the item's
# scale by subtracting its fitted additive/interaction shift
# from the person location, then form the same class-interval display used by
# an ordinary ICC. No facet cell is privileged as the graphical reference.
.plot_mfrm_item_icc <- function(fit, item, group, n_groups, grid, observed) {
  item <- as.character(item)
  vm <- fit$virtual_map
  rows <- which(vm$item == item)
  if (!length(rows)) stop("no such item: ", item)
  cols <- match(vm$vkey[rows], colnames(fit$X))
  base_tau <- fit$item_thresholds$tau[fit$item_thresholds$item == item]
  if (!length(base_tau)) stop("no common-scale thresholds for item ", item)
  group_label <- NULL
  if (.role_columns(group,
                    if (is.null(fit$factors)) character(0) else names(fit$factors),
                    nrow(fit$X))) {
    if (is.null(fit$factors) || !all(group %in% names(fit$factors)))
      stop("every named group must be a person factor in the fit")
    group_label <- paste(group, collapse = " x ")
    group <- if (length(group) == 1L) fit$factors[[group]] else
      .factor_cells(fit$factors[group], sep = ":")
  }
  if (!is.null(group) && length(group) != nrow(fit$X))
    stop("group must have one value per person")
  if (is.null(n_groups))
    n_groups <- if (is.null(group)) fit$n_groups else
      .dif_n_groups(fit, group)
  stacked <- lapply(cols, function(j) {
    tau_j <- fit$tau_list[[j]]
    if (length(tau_j) != length(base_tau)) return(NULL)
    shift <- mean(tau_j) - mean(base_tau)
    data.frame(theta = fit$person$theta - shift,
               score = fit$X[, j], extreme = fit$person$extreme,
               group = if (is.null(group)) NA_character_ else
                 .role_text_values(group), stringsAsFactors = FALSE)
  })
  stacked <- do.call(rbind, Filter(Negate(is.null), stacked))
  if (is.null(stacked) || !nrow(stacked))
    stop("no comparable response cells for item ", item)
  mmax <- length(base_tau)
  expected <- vapply(grid, function(th)
    item_moments(th, base_tau)$E, 0)
  op <- .rr_canvas(range(grid), c(0, mmax),
                   "Facet-adjusted person location (logits)",
                   "Expected score",
                   if (is.null(group_label)) item else
                     sprintf("%s by %s", item, group_label))
  on.exit(par(op))
  lines(grid, expected, lwd = 3, col = .rr$ink)
  if (!isTRUE(observed)) return(invisible(NULL))
  ci <- .class_intervals(
    ifelse(is.na(stacked$score), NA_real_, stacked$theta),
    stacked$extreme, n_groups)
  ok <- !is.na(ci) & is.finite(stacked$score)
  if (is.null(group)) {
    points(tapply(stacked$theta[ok], ci[ok], mean),
           tapply(stacked$score[ok], ci[ok], mean),
           pch = 21, bg = .rr$blue, col = "white", cex = 1.5, lwd = 1.2)
    .rr_legend("topleft", c("Model", "Observed"),
               lwd = c(3, NA), pch = c(NA, 21),
               pt.bg = c(NA, .rr$blue), col = c(.rr$ink, "white"),
               pt.cex = 1.4)
  } else {
    g <- droplevels(factor(stacked$group[ok]))
    cols_g <- rep(.rr$pal, length.out = nlevels(g))
    for (k in seq_len(nlevels(g))) {
      take <- g == levels(g)[k]
      if (sum(take) < 2L) next
      xx <- tapply(stacked$theta[ok][take], ci[ok][take], mean)
      yy <- tapply(stacked$score[ok][take], ci[ok][take], mean)
      lines(xx, yy, col = cols_g[k], lwd = 1.4, lty = 3)
      points(xx, yy, pch = 21, bg = cols_g[k], col = "white",
             cex = 1.4, lwd = 1.1)
    }
    .rr_legend("topleft", levels(g), lwd = 1.4, lty = 3, pch = 21,
               pt.bg = cols_g, col = cols_g, pt.cex = 1.25)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Item characteristic curve, with optional group overlay (the graphical DIF
# display).
# ---------------------------------------------------------------------------
#' Plot an item characteristic curve
#'
#' Draws the model expected-score curve with observed class-interval means
#' overlaid. Several items may be drawn together; their expected scores are
#' then expressed as proportions of their maximum scores. With \code{group}
#' supplied, observed means are drawn separately per group, the conventional
#' graphical DIF display. For an MFRM fit, a single item may be named; its
#' observed item-by-facet response cells are aligned by their fitted facet and
#' interaction shifts before the class-interval means are formed.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item One or more item names or column indices. Up to eight items may
#'   be overlaid. A group overlay requires a single item.
#' @param group Optional person grouping vector, or one or more names of
#'   factors nominated in the fit, for a DIF overlay; several names give
#'   the factor-combination cells (the factorial display).
#' @param n_groups Number of class intervals for the observed means; by
#'   default the fit's own count, or -- with a \code{group} overlay -- a
#'   count adapted to keep the smallest group's interval cells adequately
#'   filled.
#' @param grid Logit grid over which to draw the model curve.
#' @param observed Whether to add observed class-interval means.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- sprintf("I%02d", 1:6)
#' plot_icc(rasch(X), "I03")
#' @export
plot_icc <- function(fit, item, group = NULL, n_groups = NULL,
                     grid = NULL, observed = TRUE) {
  .check_flag(observed, "observed")
  grid <- .model_grid(fit, grid, half_width = 5)
  n_groups_given <- !is.null(n_groups)
  if (!is.null(n_groups)) n_groups <- .check_whole(n_groups, "n_groups", 2)
  if (!is.null(group)) {
    if (.role_columns(group,
                      if (is.null(fit$factors)) character(0) else names(fit$factors),
                      nrow(fit$X))) {
      unknown <- if (is.null(fit$factors)) group else
        setdiff(group, names(fit$factors))
      if (length(unknown))
        stop("`group` names no fitted person factor: ",
             paste(unknown, collapse = ", "),
             if (!is.null(fit$factors)) paste0("; fitted factor(s): ",
               paste(names(fit$factors), collapse = ", ")) else
               "; the fit carries no person factors")
    } else if (length(group) != nrow(fit$X))
      stop("`group` must name fitted factor(s) or give one value per ",
           "person (", nrow(fit$X), ")")
  }
  if (inherits(fit, "rasch_mfrm") &&
      (is.character(item) || is.factor(item)) &&
      length(item) == 1L && !item %in% fit$items$item &&
      as.character(item) %in% fit$virtual_map$item)
    return(.plot_mfrm_item_icc(fit, as.character(item), group, n_groups,
                               grid, observed))
  i <- .item_idx(fit, item)
  if (anyDuplicated(i))
    stop("An item characteristic curve cannot name the same item more than once",
         call. = FALSE)
  if (!length(i) || anyNA(i) || any(i < 1L | i > nrow(fit$items)))
    stop("Every item must name a fitted item", call. = FALSE)
  if (length(i) > 8L)
    stop("At most eight item characteristic curves may be overlaid", call. = FALSE)
  if (length(i) > 1L && !is.null(group))
    stop("A group overlay requires a single item", call. = FALSE)
  if (length(i) > 1L) {
    if (is.null(n_groups)) n_groups <- fit$n_groups
    th <- fit$person$theta
    ex <- if (!is.null(fit$person$extreme)) fit$person$extreme else
      rep(FALSE, length(th))
    cols <- .rr$pal[seq_along(i)]
    op <- .rr_canvas(range(grid), c(0, 1), "Person location (logits)",
                     "Expected score (proportion of maximum)",
                     "Item characteristic curves", right = 2)
    on.exit(par(op))
    for (j in seq_along(i)) {
      ij <- i[j]
      tau_j <- fit$tau_list[[ij]]
      max_j <- length(tau_j)
      curve <- vapply(grid, function(theta)
        item_moments(theta, tau_j, disc = .disc_of(fit, ij))$E, 0) / max_j
      lines(grid, curve, lwd = 2.6, col = cols[j])
      if (isTRUE(observed)) {
        x <- fit$X[, ij]
        ci_full <- .fit_class_intervals(
          fit, ij, x, th, ex, n_groups, n_groups_given)
        ok <- !is.na(ci_full)
        ci <- ci_full[ok]
        obs_th <- tapply(th[ok], ci, mean)
        obs_x <- tapply(x[ok] / max_j, ci, mean)
        points(obs_th, obs_x, pch = 21, bg = cols[j], col = "white",
               cex = 1.15, lwd = 1)
      }
    }
    labs <- fit$items$item[i]
    .rr_legend("topleft", labs, lwd = 2.6, col = cols,
               cex = if (length(labs) > 6L) 0.72 else 0.8)
    if (isTRUE(observed))
      .rr_legend("bottomright", c("Model", "Observed"),
                 lwd = c(2.6, NA), pch = c(NA, 21),
                 pt.bg = c(NA, .rr$ink), col = c(.rr$ink, "white"),
                 pt.cex = 1.15)
    return(invisible(NULL))
  }
  tau_i <- fit$tau_list[[i]]; mmax <- length(tau_i)
  if (.role_columns(group,
                    if (is.null(fit$factors)) character(0) else names(fit$factors),
                    nrow(fit$X)) &&
      !is.null(fit$factors) && all(group %in% names(fit$factors)))
    group <- if (length(group) == 1L) fit$factors[[group]] else
      .factor_cells(fit$factors[group], sep = ":")
  if (is.null(n_groups))
    n_groups <- if (is.null(group)) fit$n_groups else
      .dif_n_groups(fit, group)
  Ecurve <- vapply(grid, function(th)
    item_moments(th, tau_i, disc = .disc_of(fit, i))$E, 0)
  th <- fit$person$theta; x <- fit$X[, i]
  # the observed points use the SAME class intervals as the fit's item-trait
  # test: extreme-score persons excluded, and tied locations kept together
  # so the allocation is order-invariant (a plain rank/cut split ties by row
  # order and would move points when the data are merely reordered)
  ex <- if (!is.null(fit$person$extreme)) fit$person$extreme else
    rep(FALSE, length(th))
  ci_full <- .fit_class_intervals(fit, i, x, th, ex, n_groups, n_groups_given,
                                  group)
  ok <- !is.na(ci_full)
  ci <- ci_full[ok]
  op <- .rr_canvas(range(grid), c(0, mmax), "Person location (logits)",
                   "Expected score",
                   sprintf("%s  (location %.3f)", fit$items$item[i],
                           fit$items$location[i]))
  on.exit(par(op))
  lines(grid, Ecurve, lwd = 3, col = .rr$ink)
  if (is.null(group) && isTRUE(observed)) {
    obsTh <- tapply(th[ok], ci, mean); obsX <- tapply(x[ok], ci, mean)
    points(obsTh, obsX, pch = 21, bg = .rr$blue, col = "white", cex = 1.5, lwd = 1.2)
    .rr_legend("topleft", c("Model", "Observed"),
               lwd = c(3, NA), pch = c(NA, 21), pt.bg = c(NA, .rr$blue),
               col = c(.rr$ink, "white"), pt.cex = 1.4)
  } else if (!is.null(group) && isTRUE(observed)) {
    g <- factor(group)[ok]
    levs <- levels(droplevels(g))
    for (li in seq_along(levs)) {
      sel <- g == levs[li]
      obsTh <- tapply(th[ok][sel], ci[sel], mean)
      obsX <- tapply(x[ok][sel], ci[sel], mean)
      colr <- .rr$pal[(li - 1L) %% length(.rr$pal) + 1L]
      lines(obsTh, obsX, col = colr, lwd = 1.4, lty = 3)
      points(obsTh, obsX, pch = 21, bg = colr, col = "white", cex = 1.5, lwd = 1.2)
    }
    legend_cols <- .rr$pal[(seq_along(levs) - 1L) %% length(.rr$pal) + 1L]
    .rr_legend("topleft", levs, lwd = 1.4, lty = 3, pch = 21,
               pt.bg = legend_cols, col = legend_cols, pt.cex = 1.3)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Category probability curves.
# ---------------------------------------------------------------------------
#' Plot category probability curves
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @param grid Logit grid over which to draw the curves.
#' @param observed Whether to add category proportions by class interval
#'   (Andrich and Marais 2019, ch. 20).
#' @param n_groups Class intervals for the observed points; by default the
#'   fit's own allocation for this item, so the points match the item-trait
#'   test they illustrate.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(1:4, function(i)
#'   sapply(th, function(t)
#'     sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
#' colnames(X) <- sprintf("P%02d", 1:4)
#' plot_ccc(rasch(X), "P01", observed = TRUE)
#' @export
plot_ccc <- function(fit, item, grid = NULL, observed = FALSE,
                     n_groups = NULL) {
  if (length(item) != 1L) stop("`item` must name exactly one item")
  .check_flag(observed, "observed")
  grid <- .model_grid(fit, grid)
  n_groups_given <- !is.null(n_groups)
  n_groups <- if (n_groups_given) .check_whole(n_groups, "n_groups", 2)
              else fit$n_groups
  i <- .item_idx(fit, item); tau_i <- fit$tau_list[[i]]; mmax <- length(tau_i)
  P <- vapply(grid, function(th)
    item_moments(th, tau_i, disc = .disc_of(fit, i))$P, numeric(mmax + 1))
  ordered <- all(diff(tau_i) > 0) || mmax == 1L
  op <- .rr_canvas(range(grid), c(0, 1), "Person location (logits)",
                   "Category probability",
                   fit$items$item[i])
  on.exit(par(op))
  abline(v = tau_i, lty = 3, col = .rr$soft)
  for (cat in 0:mmax)
    lines(grid, P[cat + 1, ], lwd = 2.6,
          col = .rr$pal[cat %% length(.rr$pal) + 1L])
  if (observed) {
    th <- fit$person$theta; x <- fit$X[, i]
    ex <- fit$person$extreme %||% rep(FALSE, length(th))
    ci <- .fit_class_intervals(fit, i, x, th, ex, n_groups, n_groups_given)
    ok <- !is.na(ci)
    obsTh <- tapply(th[ok], ci[ok], mean)
    for (cat in 0:mmax) {
      obsP <- tapply(x[ok] == cat, ci[ok], mean)
      points(obsTh, obsP, pch = 21, cex = 1.2, lwd = 1.1, col = "white",
             bg = .rr$pal[cat %% length(.rr$pal) + 1L])
    }
  }
  mtext(if (ordered) "thresholds ordered" else "THRESHOLDS DISORDERED",
        side = 3, line = 0.2, adj = 0, cex = 0.8,
        col = if (ordered) .rr$teal else .rr$red)
  .rr_legend("right", paste0("Category ", 0:mmax), lwd = 2.6,
             col = .rr$pal[(0:mmax) %% length(.rr$pal) + 1L])
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Threshold probability curves: conditional adjacent-category ogives, each
# crossing 0.5 at its threshold.
# ---------------------------------------------------------------------------
#' Plot threshold probability curves
#'
#' Conditional probability of success at each threshold,
#' \code{P(X = k | X = k - 1 or k)}, a logistic ogive crossing 0.5 at the
#' threshold location. Disordered thresholds are immediately visible as
#' out-of-sequence ogives. With \code{observed = TRUE} the observed
#' conditional proportions per class interval are overlaid, the direct
#' check on whether each threshold discriminates (and hence whether
#' collapsing categories could ever be justified; Andrich and Marais 2019,
#' ch. 22).
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @param grid Logit grid over which to draw the curves.
#' @param observed Overlay the observed conditional threshold proportions
#'   per class interval.
#' @param n_groups Class intervals for the observed points; by default the
#'   fit's own allocation for this item, so the points match the item-trait
#'   test they illustrate.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(1:4, function(i)
#'   sapply(th, function(t)
#'     sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
#' colnames(X) <- sprintf("P%02d", 1:4)
#' plot_threshold_prob(rasch(X), "P01")
#' @export
plot_threshold_prob <- function(fit, item, grid = NULL,
                                observed = FALSE, n_groups = NULL) {
  if (length(item) != 1L) stop("`item` must name exactly one item")
  grid <- .model_grid(fit, grid)
  .check_flag(observed, "observed")
  n_groups_given <- !is.null(n_groups)
  n_groups <- if (n_groups_given) .check_whole(n_groups, "n_groups", 2)
              else fit$n_groups
  i <- .item_idx(fit, item); tau_i <- fit$tau_list[[i]]
  op <- .rr_canvas(range(grid), c(0, 1), "Person location (logits)",
                   "Threshold probability",
                   fit$items$item[i])
  on.exit(par(op))
  abline(h = 0.5, lty = 2, col = .rr$soft)
  for (k in seq_along(tau_i)) {
    colr <- .rr$pal[(k - 1L) %% length(.rr$pal) + 1L]
    lines(grid, plogis(.disc_of(fit, i) * (grid - tau_i[k])), lwd = 2.6, col = colr)
    points(tau_i[k], 0.5, pch = 21, bg = colr, col = "white", cex = 1.4)
  }
  if (observed) {
    # observed conditional threshold proportions per class interval:
    # among persons responding k - 1 or k, the proportion responding k
    # (Andrich & Marais 2019, ch. 22 rescoring check)
    th <- fit$person$theta; x <- fit$X[, i]
    ex <- fit$person$extreme %||% rep(FALSE, length(th))
    ci <- .fit_class_intervals(fit, i, x, th, ex, n_groups, n_groups_given)
    ok <- !is.na(ci)
    for (k in seq_along(tau_i)) {
      colr <- .rr$pal[(k - 1L) %% length(.rr$pal) + 1L]
      inpair <- ok & (x == k - 1L | x == k)
      cip <- ci[inpair]
      if (!sum(inpair)) next
      obsTh <- tapply(th[inpair], cip, mean)
      obsT <- tapply(x[inpair] == k, cip, mean)
      points(obsTh, obsT, pch = 21, cex = 1.2, lwd = 1.1, col = "white",
             bg = colr)
    }
  }
  .rr_legend(
    "topleft", sprintf("threshold %d (%.3f)", seq_along(tau_i), tau_i),
    lwd = 2.6,
    col = .rr$pal[(seq_along(tau_i) - 1L) %% length(.rr$pal) + 1L])
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Person-item threshold distribution: mirrored histograms on one logit scale.
# ---------------------------------------------------------------------------
#' Plot the person-item threshold distribution
#'
#' The targeting display: the person location distribution above the axis and
#' the calibration threshold distribution mirrored below it, on a shared
#' logit scale. Dashed lines mark the person and threshold means in their
#' distributions' colours. MFRM and EFRM thresholds belong to response cells.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param bins Number of histogram bins.
#' @param xlim Optional logit range for the shared scale; persons and
#'   thresholds outside it are omitted. By default the range is extended to
#'   labelled tick marks beyond the most extreme plotted estimate.
#' @param information Whether to overlay the test information function on a
#'   separate right-hand axis. Fits with more than one administrable design
#'   receive one curve per design.
#' @param group Optional person-group level: one level of a fitted person
#'   factor, restricting the person distribution to those persons. A level
#'   no fitted factor carries is an error.
#' @param items Optional item selection restricting the threshold
#'   distribution: item names, or one item-set name of an extended-frame
#'   fit, whose virtual item-by-group cells match through their underlying
#'   items. The selection is named in the legend, so a restricted map
#'   cannot be read as the whole instrument. The information curve follows
#'   the same item selection; for an extended-frame fit it also follows the
#'   response cells occupied by the selected person group. For EFRM and
#'   MFRM, only response patterns present in that person group are shown.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_pimap(rasch(X))
#' @export
plot_pimap <- function(fit, bins = 35, xlim = NULL, information = FALSE,
                       group = NULL, items = NULL) {
  .check_targeting_fit(fit)
  bins <- .check_whole(bins, "bins", 2)
  .check_flag(information, "information")
  .check_xlim(xlim)
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  # the two sides can be restricted independently: a person group against
  # the whole instrument, or the whole sample against one item set. The
  # subset is named on the plot, so a restricted map cannot be mistaken for
  # the full one.
  keep_p <- .pimap_persons(fit, group)
  keep_i <- .pimap_items(fit, items)
  info_idx <- NULL
  if (isTRUE(information)) {
    info_idx <- if (is.null(items)) NULL else
      unique(fit$thresholds$item[keep_i])
    if (!is.null(group) && inherits(fit, "rasch_efrm") &&
        !is.null(fit$virtual_map) &&
        fit$frame_group[1] %in% names(fit$factors)) {
      # Read the actual frames occupied by the selected persons.  A requested
      # item-by-frame cell from another group is contradictory, not a request
      # to substitute every response cell in the selected group.
      frames_of <- as.character(fit$factors[[fit$frame_group[1]]])
      glev <- unique(frames_of[keep_p])
      gcols <- which(as.character(fit$virtual_map$group) %in% glev)
      if (length(gcols)) {
        if (is.null(info_idx)) info_idx <- gcols
        else {
          info_idx <- intersect(info_idx, gcols)
          if (!length(info_idx))
            stop("`items` and `group` select no common EFRM response cells; ",
                 "select an underlying item or a response cell in the chosen ",
                 "person group", call. = FALSE)
        }
      }
    }
  }
  th <- fit$person$theta[keep_p]
  th <- th[!is.na(th)]
  tau <- fit$thresholds$tau[keep_i]
  scale <- .pimap_scale(c(th, tau), xlim)
  rng <- scale$range
  th <- th[th >= rng[1] & th <= rng[2]]
  tau <- tau[tau >= rng[1] & tau <= rng[2]]
  .check_map_range(length(tau), length(th))
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  hp <- hist(th, breaks = brk, plot = FALSE)
  hi <- hist(tau, breaks = brk, plot = FALSE)
  pp <- hp$counts / sum(hp$counts); pi <- hi$counts / sum(hi$counts)
  # both sides end on a labelled tick: the smallest pretty proportion above
  # the tallest bar on each side, mapped exactly so the axis cannot top out
  # between labels
  at <- pretty(c(0, max(c(pp, pi))))
  step <- diff(at)[1]
  top <- at[at > max(pp)][1]
  if (is.na(top)) top <- max(at) + step
  bot <- at[at > max(pi)][1]
  if (is.na(bot)) bot <- max(at) + step
  at <- unique(sort(c(at, top, bot)))
  ymax <- top; ymin <- -bot
  op2 <- par(yaxs = "i")
  on.exit(par(op2), add = TRUE)
  op <- .rr_canvas(rng, c(ymin, ymax), "Location (logits)", "Proportion",
                   "", grid_y = FALSE,
                   yaxis = FALSE, right = if (isTRUE(information)) 4.2 else 1.5,
                   xaxis = FALSE)
  on.exit(par(op), add = TRUE)
  axis(1, at = scale$ticks, col = .rr$grid, col.ticks = .rr$soft)
  axis(2, at = c(-rev(at[-1]), at), labels = c(rev(at[-1]), at),
       col = .rr$grid, col.ticks = .rr$soft, cex.axis = 0.85)
  # continuous spines meeting at the corner, instead of axis lines that
  # stop at their outermost tick marks
  usr <- par("usr")
  axis(1, at = usr[1:2], labels = FALSE, lwd.ticks = 0, col = .rr$grid)
  axis(2, at = usr[3:4], labels = FALSE, lwd.ticks = 0, col = .rr$grid)
  abline(h = 0, col = .rr$ink, lwd = 1)
  rect(brk[-length(brk)], 0, brk[-1], pp, col = .rr$blue, border = "white", lwd = 0.6)
  rect(brk[-length(brk)], -pi, brk[-1], 0, col = .rr$amber, border = "white", lwd = 0.6)
  segments(mean(th), 0, mean(th), ymax * 0.95, col = .rr$blue, lty = 2)
  segments(mean(tau), ymin * 0.95, mean(tau), 0, col = .rr$amber, lty = 2)
  map_labels <- c(
    if (is.null(group)) "Persons" else paste0("Persons: ", group),
    paste0(if (structural) "Calibration thresholds" else "Item thresholds",
           if (!is.null(items)) paste0(": ", .pimap_item_label(items)) else ""))
  .rr_legend("topleft", map_labels,
             fill = c(.rr$blue, .rr$amber), border = NA, cex = 0.76)
  if (isTRUE(information)) {
    grid <- seq(rng[1], rng[2], length.out = 241L)
    # the curve must describe the same design as the distributions below it:
    # a set-restricted map carrying the full instrument's information, or a
    # one-group EFRM map carrying every group's curves, would read as
    # precision the selection does not have
    information_fit <- fit
    if (structural && !is.null(group)) {
      # Retain the calibration, but identify administrations only among the
      # selected persons. Restricting frame columns alone cannot distinguish
      # subgroups who answered different items within the same frame.
      information_fit$X <- fit$X[keep_p %in% TRUE, , drop = FALSE]
    }
    ti <- test_information(information_fit, grid, items = info_idx)
    des <- if ("design" %in% names(ti)) unique(ti$design) else "Test information"
    # The palette limit belongs to the plot, not to the statistics: past it,
    # draw every design in one soft colour and count them in the legend
    # rather than cycling colours that imply distinctions no reader can
    # follow back to a curve.
    pal <- c(.rr$teal, .rr$purple, .rr$red, .rr$amber)
    named <- length(des) <= length(pal)
    cols <- if (named) rep_len(pal, length(des)) else
      rep(.rr$soft, length(des))
    imax <- max(ti$info, na.rm = TRUE)
    if (is.finite(imax) && imax > 0) {
      scl <- ymax * 0.92 / imax
      for (j in seq_along(des)) {
        z <- if ("design" %in% names(ti)) ti$design == des[j] else
          rep(TRUE, nrow(ti))
        lines(ti$theta[z], ti$info[z] * scl, lwd = 2.5, col = cols[j])
      }
      ticks <- pretty(c(0, imax))
      ticks <- ticks[ticks >= 0 & ticks * scl <= ymax]
      axis(4, at = ticks * scl, labels = ticks, col = .rr$grid,
           col.ticks = .rr$soft, col.axis = .rr$teal, cex.axis = 0.8)
      mtext("Test information", side = 4, line = 2.7,
            col = .rr$teal, cex = 0.85, las = 0)
      if (length(des) > 1L) {
        # Long crossed-frame labels can meet the distribution legend in a
        # narrow window. Use one title rather than repeating a prefix, scale
        # only as far as remains legible, and stack the legends when they do
        # not fit side by side.
        labs <- if (named) des else
          sprintf("%d administrations; more than the palette can name",
                  length(des))
        usr <- par("usr"); span <- diff(usr[1:2])
        info_cex <- 0.68
        info_w <- max(strwidth(c("Test information", labs), cex = info_cex,
                               units = "user"))
        if (info_w > 0.46 * span)
          info_cex <- max(0.50, info_cex * 0.46 * span / info_w)
        map_w <- max(strwidth(map_labels, cex = 0.76, units = "user"))
        info_w <- max(strwidth(c("Test information", labs), cex = info_cex,
                               units = "user"))
        down <- if (map_w + info_w > 0.82 * span) 0.17 else 0
        .rr_legend("topright", labs, title = "Test information",
                   lwd = 2.5, col = if (named) cols else .rr$soft,
                   cex = info_cex, inset = c(0, down))
      }
    }
  }
  invisible(NULL)
}

# Restricting the person side: `group` names one level of a fitted person
# factor. A level no factor carries is refused rather than silently drawing
# the whole sample.
.pimap_persons <- function(fit, group) {
  n <- nrow(fit$person)
  if (is.null(group)) return(rep(TRUE, n))
  if (!is.atomic(group) || !is.null(dim(group)) || length(group) != 1L ||
      is.na(group))
    stop("`group` must name exactly one person-group level", call. = FALSE)
  fac <- fit$factors
  if (is.null(fac) || !ncol(fac))
    stop("the fit carries no person factors, so it has no groups to select",
         call. = FALSE)
  group <- as.character(group)
  # the qualified form names the factor as well as the level, and is the
  # only unambiguous address when two factors share a level name. The split
  # is on a fitted factor's NAME followed by a colon -- not on any colon --
  # so a level (or factor) whose own name carries one still resolves; the
  # longest matching factor name wins when one name prefixes another.
  pref <- names(fac)[startsWith(group, paste0(names(fac), ":"))]
  if (length(pref)) {
    fname <- pref[which.max(nchar(pref))]
    level <- trimws(substring(group, nchar(fname) + 2L))
    v <- as.character(fac[[fname]])
    if (!level %in% unique(v))
      stop("'", level, "' is not a level of factor '", fname, "'",
           call. = FALSE)
    return(v == level)
  }
  if (grepl(":", group, fixed = TRUE) &&
      !group %in% unlist(lapply(fac, function(v) as.character(unique(v)))))
    stop("`group` '", group, "' does not match 'factor: level'; fitted ",
         "factors: ", paste(names(fac), collapse = ", "), call. = FALSE)
  hit <- vapply(fac, function(v) group %in% as.character(unique(v)), TRUE)
  if (!any(hit))
    stop("`group` value '", group, "' is not a level of any fitted person ",
         "factor; available: ",
         paste(unique(unlist(lapply(fac, function(v)
           as.character(unique(v))))), collapse = ", "), call. = FALSE)
  if (sum(hit) > 1L)
    stop("level '", group, "' belongs to several factors; name it as one of: ",
         paste(sprintf("'%s: %s'", names(fac)[hit], group), collapse = ", "),
         call. = FALSE)
  as.character(fac[[which(hit)]]) == group
}

# Restricting the item side: item names, or one item-set name of an EFRM fit.
# An EFRM calibrates item-by-group cells, so its fitted "items" are virtual
# keys; the set map and the underlying item names are recovered through the
# virtual map, and both name forms select.
.pimap_items <- function(fit, items) {
  fitted_nm <- fit$items$item[fit$thresholds$item]
  vm <- fit$virtual_map
  base_nm <- if (!is.null(vm))
    as.character(vm$item)[match(fitted_nm, as.character(vm$vkey))]
  else fitted_nm
  base_nm[is.na(base_nm)] <- fitted_nm[is.na(base_nm)]
  if (is.null(items)) return(rep(TRUE, length(fitted_nm)))
  if (!is.character(items) || !length(items) || anyNA(items) ||
      any(!nzchar(trimws(items))))
    stop("`items` must name item(s), or one item set", call. = FALSE)
  if (anyDuplicated(items))
    stop("`items` must not name the same item more than once", call. = FALSE)
  sets <- fit$set_of
  if (length(items) == 1L && !is.null(sets) &&
      items %in% as.character(sets))
    return(base_nm %in% names(sets)[as.character(sets) == items])
  if (length(items) == 1L && !is.null(vm) && items %in% as.character(vm$set))
    return(base_nm %in% as.character(vm$item)[as.character(vm$set) == items])
  known <- unique(c(fitted_nm, base_nm))
  miss <- setdiff(items, known)
  if (length(miss))
    stop("item(s) not in the fit: ", paste(miss, collapse = ", "),
         call. = FALSE)
  fitted_nm %in% items | base_nm %in% items
}

.pimap_item_label <- function(items)
  if (length(items) == 1L) items else paste0(length(items), " items")

.pimap_scale <- function(values, xlim = NULL) {
  if (is.null(xlim)) {
    padded <- range(values[is.finite(values)]) + c(-0.4, 0.4)
    w <- diff(padded)
    ticks <- if (w > 2.5 && w <= 14)
      seq(floor(padded[1] + 1e-9), ceiling(padded[2] - 1e-9))
    else pretty(padded, n = 6)
    return(list(range = range(ticks), ticks = ticks))
  }
  rng <- sort(xlim)
  ticks <- .rr_ticks(rng)
  ticks <- ticks[ticks >= rng[1] & ticks <= rng[2]]
  list(range = rng, ticks = ticks)
}

# ---------------------------------------------------------------------------
# Wright map: the conventional vertical person-item map (Wright & Stone
# 1979), person distribution left of a common logit axis, item threshold
# labels stacked to its right.
# ---------------------------------------------------------------------------
#' Plot a Wright map
#'
#' The conventional vertical person-item map (Wright and Stone 1979): the
#' person distribution to the left of a shared logit axis and the calibration
#' thresholds stacked to its right. Dashed lines mark the person and threshold
#' means in their distributions' colours. MFRM and EFRM labels identify
#' item-by-facet or item-by-frame response cells rather than additional items.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param bins Number of bins for the person distribution and the threshold
#'   label rows.
#' @param xlim Optional logit range for the shared scale; persons and
#'   thresholds outside it are omitted.
#' @param cex_labels Character expansion for the threshold labels.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @references Wright, B. D., & Stone, M. H. (1979). \emph{Best Test
#'   Design}. Chicago: MESA Press.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_wright(rasch(X))
#' @export
plot_wright <- function(fit, bins = 35, xlim = NULL, cex_labels = 0.8) {
  .check_targeting_fit(fit)
  bins <- .check_whole(bins, "bins", 2)
  .check_xlim(xlim)
  .check_cex(cex_labels)
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  th <- fit$person$theta[!is.na(fit$person$theta)]
  thr <- fit$thresholds
  rng <- if (is.null(xlim)) range(c(th, thr$tau)) + c(-0.4, 0.4) else sort(xlim)
  th <- th[th >= rng[1] & th <= rng[2]]
  keep <- thr$tau >= rng[1] & thr$tau <= rng[2]
  tv <- thr$tau[keep]
  lab <- if (all(fit$m == 1L)) fit$items$item[thr$item[keep]] else
    paste0(fit$items$item[thr$item[keep]], ".", thr$k[keep])
  .check_map_range(length(tv), length(th))
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  hp <- hist(th, breaks = brk, plot = FALSE)
  pp <- hp$counts / max(1L, max(hp$counts))
  bin <- pmin(findInterval(tv, brk, rightmost.closed = TRUE), bins)
  split <- 0.42
  op <- par(mar = c(1.2, 4.4, 1.6, 0.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = rng, xlab = "", ylab = "Location (logits)",
       axes = FALSE, main = "")
  abline(h = .rr_ticks(rng), col = .rr$grid, lwd = 0.8)
  .rr_axis(2)
  segments(split, rng[1], split, rng[2], col = .rr$ink, lwd = 1)
  segments(split, mean(tv), 0.99, mean(tv), col = .rr$amber, lty = 2)
  rect(split - pp * (split - 0.02), brk[-length(brk)], split, brk[-1],
       col = .rr$blue, border = "white", lwd = 0.6)
  segments(0.02, mean(th), split, mean(th), col = .rr$blue, lty = 2)
  cexl <- cex_labels
  colw <- max(strwidth(lab, cex = cexl)) * 1.2
  x0 <- split + 0.015
  kmax <- max(1L, floor((1 - x0) / colw))
  for (b in unique(bin)) {
    ls <- lab[bin == b][order(tv[bin == b])]
    if (length(ls) > kmax)
      ls <- c(ls[seq_len(kmax - 1L)], paste0("+", length(ls) - kmax + 1L))
    text(x0 + (seq_along(ls) - 1L) * colw, (brk[b] + brk[b + 1L]) / 2, ls,
         cex = cexl, adj = 0, col = .rr$ink)
  }
  text(0.02, rng[2], sprintf("persons: mean %.2f, SD %.2f",
                             mean(th), stats::sd(th)),
       cex = 0.78, adj = c(0, 1), col = .rr$blue, font = 2)
  text(0.99, rng[2], sprintf("%sthresholds: mean %.2f, SD %.2f",
                             if (structural) "calibration " else "",
                             mean(tv), stats::sd(tv)),
       cex = 0.78, adj = c(1, 1), col = .rr$amber, font = 2)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Threshold map: every item's thresholds on the logit scale, by location.
# ---------------------------------------------------------------------------
#' Plot the threshold map
#'
#' Each fitted column's threshold locations on a common logit scale, ordered
#' by location, with disordered thresholds highlighted. The columns are
#' response cells for MFRM and EFRM fits.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param order_by_location Order items by their location (the default)
#'   rather than their original sequence.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_threshold_map(rasch(X))
#' @export
plot_threshold_map <- function(fit, order_by_location = TRUE) {
  .check_response_display_fit(fit, "threshold plots")
  .check_flag(order_by_location, "order_by_location")
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  L <- length(fit$tau_list)
  ord <- if (order_by_location) order(fit$items$location) else seq_len(L)
  rng <- range(fit$thresholds$tau); rng <- rng + c(-0.4, 0.4)
  op <- par(mar = c(4.2, 7.5, 3.2, 1.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = rng, ylim = c(0.5, L + 0.5), xlab = "Location (logits)",
       ylab = "", axes = FALSE, main = "")
  abline(h = seq_len(L), col = .rr$grid, lwd = 0.8)
  abline(v = 0, lty = 2, col = .rr$soft)
  .rr_axis(1)
  axis(2, at = seq_len(L), labels = fit$items$item[ord], cex.axis = 0.75,
       col = .rr$grid, col.ticks = NA)
  for (row in seq_len(L)) {
    i <- ord[row]; tau_i <- fit$tau_list[[i]]
    disord <- !(all(diff(tau_i) > 0) || length(tau_i) == 1L)
    segments(min(tau_i), row, max(tau_i), row, col = .rr$soft, lwd = 1.4)
    points(tau_i, rep(row, length(tau_i)), pch = 21, cex = 1.25,
           bg = if (disord) .rr$red else .rr$amber, col = "white", lwd = 1)
    points(mean(tau_i), row, pch = 23, bg = .rr$blue, col = "white", cex = 1.3)
  }
  .rr_legend("bottomright", c(if (structural) "calibration threshold" else
    "threshold", "disordered", if (structural) "response-cell location" else
      "item location"),
             pch = c(21, 21, 23), pt.bg = c(.rr$amber, .rr$red, .rr$blue),
             col = "white", pt.cex = 1.2)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Test characteristic curve and test information.
# ---------------------------------------------------------------------------
#' Plot the test characteristic curve
#'
#' Expected total score against person location. Structural fits draw one
#' curve for each observed item pattern within the frame or facet design,
#' as defined by \code{\link{test_information}}.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param grid Logit grid.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_tcc(rasch(X))
#' @export
plot_tcc <- function(fit, grid = NULL) {
  grid <- .model_grid(fit, grid)
  # the same administrable design blocks as test_information(), so the two
  # displays can never disagree about which items form a curve
  blocks <- .design_blocks(fit)
  curves <- lapply(blocks, function(ii) vapply(grid, function(th)
    sum(vapply(ii, function(i)
      item_moments(th, fit$tau_list[[i]], disc = .disc_of(fit, i))$E, 0)), 0))
  Smax <- max(vapply(blocks, function(ii) sum(fit$m[ii]), 0))
  op <- .rr_canvas(range(grid), c(0, Smax), "Person location (logits)",
                   "Expected total score")
  on.exit(par(op))
  # past the palette limit, one soft colour and a count: see plot_pimap()
  named <- length(curves) <= length(.rr$pal)
  cols <- if (named) rep_len(.rr$pal, length(curves)) else
    rep(.rr$soft, length(curves))
  for (j in seq_along(curves)) lines(grid, curves[[j]], lwd = 3, col = cols[j])
  abline(h = c(0, Smax), lty = 3, col = .rr$soft)
  if (length(curves) > 1L)
    .rr_legend("topleft", if (named) names(curves) else
      sprintf("%d administrations; more than the palette can name",
              length(curves)),
      lwd = 3, col = if (named) cols else .rr$soft)
  invisible(NULL)
}

#' Plot the test information function
#'
#' Test information across the logit scale with the standard error of
#' measurement overlaid on a second axis.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param grid Logit grid.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_tif(rasch(X))
#' @export
plot_tif <- function(fit, grid = NULL) {
  grid <- .model_grid(fit, grid)
  ti <- test_information(fit, grid)
  if (!any(is.finite(ti$info) & ti$info > 0))
    stop("the requested grid contains no positive test information; use locations nearer the calibrated range")
  if ("design" %in% names(ti) && length(unique(ti$design)) > 1L) {
    des <- unique(ti$design)
    # past the palette limit, one soft colour and a count: see plot_pimap()
    named <- length(des) <= length(.rr$pal)
    cols <- if (named) rep_len(.rr$pal, length(des)) else
      rep(.rr$soft, length(des))
    ymax <- max(ti$info, na.rm = TRUE) * 1.1
    op <- .rr_canvas(range(grid), c(0, ymax), "Person location (logits)",
                     "Test information", right = 3.6)
    on.exit(par(op))
    sem <- ti$sem
    sem[!is.finite(sem)] <- NA_real_
    grid_mid <- mean(range(ti$theta))
    inr <- abs(ti$theta - grid_mid) < 4
    if (!any(inr & is.finite(sem))) inr <- is.finite(sem)
    sem_max <- max(sem[inr], na.rm = TRUE)
    scl <- if (is.finite(sem_max) && sem_max > 0)
      max(ti$info, na.rm = TRUE) * 1.05 / sem_max else NA_real_
    for (j in seq_along(des)) {
      z <- ti$design == des[j]
      lines(ti$theta[z], ti$info[z], lwd = 3, col = cols[j])
      if (is.finite(scl))
        lines(ti$theta[z], sem[z] * scl, lwd = 2.2, lty = 5,
              col = cols[j])
    }
    if (is.finite(scl)) {
      sem_ticks <- pretty(c(0, sem_max))
      sem_ticks <- sem_ticks[sem_ticks * scl <= ymax]
      axis(4, at = sem_ticks * scl, labels = sem_ticks,
           col = .rr$grid, col.ticks = .rr$soft, col.axis = .rr$red,
           cex.axis = 0.8)
      mtext("SEM", side = 4, line = 2.3, col = .rr$red, cex = 0.85)
    }
    .rr_legend("topleft", if (named) des else
      sprintf("%d administrations; more than the palette can name",
              length(des)),
      title = if (is.finite(scl)) "Solid: information; dashed: SEM" else NULL,
      lwd = 3, col = if (named) cols else .rr$soft)
    return(invisible(NULL))
  }
  imax <- max(ti$info[is.finite(ti$info) & ti$info > 0])
  op <- .rr_canvas(range(grid), c(0, imax * 1.1),
                   "Person location (logits)", "Test information",
                   right = 3.6)
  on.exit(par(op))
  polygon(c(ti$theta, rev(ti$theta)), c(ti$info, rep(0, nrow(ti))),
          col = paste0(.rr$blue, "22"), border = NA)
  lines(ti$theta, ti$info, lwd = 3, col = .rr$blue)
  sem <- ti$sem; sem[!is.finite(sem)] <- NA
  # the SEM axis is scaled over the central range; if the plotting grid lies
  # entirely outside it, fall back to the whole grid rather than max() over
  # an empty selection (which returns -Inf)
  grid_mid <- mean(range(ti$theta))
  inr <- abs(ti$theta - grid_mid) < 4
  if (!any(inr & is.finite(sem))) inr <- is.finite(sem)
  sem_max <- max(sem[inr], na.rm = TRUE)
  scl <- imax * 1.05 / sem_max
  lines(ti$theta, sem * scl, lwd = 2.2, col = .rr$red, lty = 5)
  sem_ticks <- pretty(c(0, sem_max))
  sem_ticks <- sem_ticks[sem_ticks * scl <= imax * 1.1]
  axis(4, at = sem_ticks * scl, labels = sem_ticks,
       col = .rr$grid, col.ticks = .rr$soft, col.axis = .rr$red, cex.axis = 0.8)
  mtext("SEM", side = 4, line = 2.3, col = .rr$red, cex = 0.85)
  .rr_legend("topleft", c("Information", "SEM"), lwd = c(3, 2.2),
             lty = c(1, 5), col = c(.rr$blue, .rr$red))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Item and person fit maps.
# ---------------------------------------------------------------------------
#' Plot the item map (location against fit residual)
#'
#' Fitted columns plotted by location and a fit statistic, with the
#' conventional acceptance band at +/- 2.5. The default statistic is the
#' log-of-mean-square fit residual; \code{"infit"} and \code{"outfit"} display
#' the Wilson--Hilferty standardised mean squares, to which the same band
#' convention applies. MFRM and EFRM points are response cells; ordinary
#' Rasch points are items.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param statistic \code{"residual"} (the default fit residual),
#'   \code{"infit"}, or \code{"outfit"}.
#' @param band Acceptance band for the standardised statistic.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_item_map(rasch(X))
#' @export
plot_item_map <- function(fit, statistic = c("residual", "infit", "outfit"),
                          band = 2.5) {
  .check_response_display_fit(fit, "item-fit plots")
  .check_band(band)
  statistic <- match.arg(statistic)
  d <- fit$items
  y <- switch(statistic, residual = d$fit_resid,
              infit = d$infit_z, outfit = d$outfit_z)
  ylab <- switch(statistic, residual = "Fit residual",
                 infit = "Infit (standardised)",
                 outfit = "Outfit (standardised)")
  ok <- is.finite(y) & is.finite(d$location)
  if (is.null(y) || !any(ok))
    .refuse("the fitted object does not carry a finite ", statistic,
            " statistic with a finite item location")
  structural <- inherits(fit, c("rasch_mfrm", "rasch_efrm"))
  ylim <- range(c(y[ok], -band, band)) * 1.2
  op <- .rr_canvas(range(d$location[ok]) + c(-0.5, 0.5), ylim,
                   if (structural) "Response-cell location (logits)" else
                     "Item location (logits)", ylab,
                   grid_x = TRUE)
  on.exit(par(op))
  rect(par("usr")[1], -band, par("usr")[2], band,
       col = paste0(.rr$teal, "11"), border = NA)
  abline(h = c(-band, band), lty = 2, col = .rr$soft)
  out <- ok & abs(y) > band
  points(d$location[ok], y[ok], pch = 21, cex = 1.7,
         bg = ifelse(out[ok], .rr$red, .rr$blue), col = "white", lwd = 1.2)
  text(d$location[ok], y[ok], d$item[ok], pos = 3, offset = 0.5, cex = 0.7,
       col = ifelse(out[ok], .rr$red, .rr$soft))
  mtext(sprintf("%d of %d %s beyond +/-%.1f (%.1f%%)", sum(out), sum(ok),
                if (structural) "response cells" else "items", band,
                100 * sum(out) / sum(ok)),
        side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$ink)
  invisible(NULL)
}

#' Plot person fit
#'
#' Person locations against a person fit statistic with the +/- 2.5 band;
#' persons beyond the band respond erratically (positive) or too
#' deterministically (negative). The default statistic is the
#' log-of-mean-square fit residual; \code{"infit"} and \code{"outfit"} display
#' the Wilson--Hilferty standardised mean squares, to which the same band
#' convention applies.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param statistic \code{"residual"} (the default fit residual),
#'   \code{"infit"}, or \code{"outfit"}.
#' @param band Acceptance band for the standardised statistic.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_person_fit(rasch(X))
#' @export
plot_person_fit <- function(fit, statistic = c("residual", "infit", "outfit"),
                            band = 2.5) {
  .check_response_display_fit(fit, "person-fit plots")
  .check_band(band)
  statistic <- match.arg(statistic)
  p <- fit$person
  y <- switch(statistic, residual = p$fit_resid,
              infit = p$infit_z, outfit = p$outfit_z)
  ylab <- switch(statistic, residual = "Fit residual",
                 infit = "Infit (standardised)",
                 outfit = "Outfit (standardised)")
  ok <- is.finite(p$theta) & is.finite(y)
  if (is.null(y) || !any(ok))
    .refuse("the fitted object does not carry a finite ", statistic,
            " statistic with a finite person location")
  ylim <- range(c(y[ok], -band - 0.5, band + 0.5))
  op <- .rr_canvas(range(p$theta[ok]) + c(-0.3, 0.3), ylim,
                   "Person location (logits)", ylab,
                   grid_x = TRUE)
  on.exit(par(op))
  rect(par("usr")[1], -band, par("usr")[2], band,
       col = paste0(.rr$teal, "11"), border = NA)
  abline(h = c(-band, band), lty = 2, col = .rr$soft)
  out <- abs(y[ok]) > band
  points(p$theta[ok], y[ok], pch = 21, cex = 0.9,
         bg = ifelse(out, .rr$red, paste0(.rr$blue, "99")), col = "white", lwd = 0.5)
  mtext(sprintf("%d of %d persons beyond +/-%.1f (%.1f%%)", sum(out), sum(ok),
                band, 100 * sum(out) / sum(ok)),
        side = 3, line = 0.2, adj = 0, cex = 0.8, col = .rr$ink)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Residual structure.
# ---------------------------------------------------------------------------
#' Plot the residual-correlation heatmap
#'
#' Only the lower triangle is drawn -- the matrix is symmetric, so each pair
#' is shown once. With \code{stat = "q3star"} (the default) cells are coloured
#' by Q3* -- each pair's residual correlation minus the average off-diagonal
#' correlation -- so white marks the value expected under local independence
#' and warm colour marks dependence; with \code{stat = "q3"} the raw residual
#' correlation is coloured, white at zero. The scale saturates at \code{cap}
#' rather than the +/-1 of an ordinary correlation: a residual correlation
#' seldom reaches even 0.5 under a fitting model, so the colour is spent where
#' the values actually discriminate. A Q3* value of 0.2 is sometimes used as
#' a heuristic screen, but it is not a universal critical value (Christensen,
#' Makransky and Horton 2017).
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param stat Which statistic to colour: \code{"q3star"} (adjusted Q3, the
#'   default) or \code{"q3"} (the raw residual correlation).
#' @param cap Value at which the colour saturates (default 0.5).
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_resid_cor(rasch(X))
#' @export
plot_resid_cor <- function(fit, stat = c("q3star", "q3"), cap = 0.5) {
  stat <- match.arg(stat)
  # the cap divides the colour-key coordinates, so it must be strictly
  # positive as well as finite
  if (length(cap) != 1L || !is.numeric(cap) || is.complex(cap) ||
      !is.null(dim(cap)) || !is.null(oldClass(cap)) ||
      !is.finite(cap) || cap <= 0)
    stop("`cap` must be one positive finite correlation")
  rc <- residual_correlations(fit)
  R <- rc$matrix; L <- ncol(R); avg <- rc$average
  # Q3* colours the excess over the average off-diagonal correlation, so white
  # marks the local-independence baseline and warm colour marks dependence;
  # Q3 colours the raw residual correlation, white at zero. Only the lower
  # triangle is drawn -- the matrix is symmetric, so each pair is read once.
  star <- stat == "q3star"
  S <- if (star) R - avg else R
  diag(S) <- NA                    # self-correlations carry no dependence
  # keep the visual lower-left triangle once image() flips the rows, so the
  # heatmap matches the lower-triangular table beside it
  S[lower.tri(S)] <- NA
  lab <- if (star) "Q3*" else "Q3"
  pal <- colorRampPalette(c("#1d4ed8", "#93c5fd", "#f8fafc",
                            "#f59e0b", "#dc2626"))(128)
  Sc <- pmax(pmin(S, cap), -cap)             # saturate at +/- cap
  op <- par(mar = c(5.5, 5.5, 1.6, 6.4), las = 1, col.axis = .rr$ink,
            col.main = .rr$ink, font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  image(1:L, 1:L, Sc[, L:1, drop = FALSE], col = pal, zlim = c(-cap, cap),
        axes = FALSE, xlab = "", ylab = "", main = "")
  cx <- if (L > 25) 0.5 else if (L > 15) 0.62 else 0.75
  axis(1, 1:L, colnames(R), las = 2, cex.axis = cx, col = NA, col.ticks = NA)
  axis(2, 1:L, rev(colnames(R)), cex.axis = cx, col = NA, col.ticks = NA)
  abline(h = 0.5 + 0:L, v = 0.5 + 0:L, col = "white", lwd = 0.6)
  # colour key: the bar is titled by the statistic and captioned so its
  # meaning stands on its own; ticks at -cap, 0, +cap
  ky <- seq(0.12, 0.80, length.out = 129) * L
  rect(L + 1.1, ky[-129], L + 1.5, ky[-1], col = pal, border = NA, xpd = TRUE)
  yof <- function(v) (v + cap) / (2 * cap) * (max(ky) - min(ky)) + min(ky)
  for (v in c(-cap, 0, cap))
    text(L + 1.6, yof(v), sprintf("%+.1f", v), cex = 0.62, xpd = TRUE,
         adj = 0, col = .rr$ink)
  text(L + 1.05, max(ky) + 0.55, lab, cex = 0.78, xpd = TRUE, adj = 0,
       col = .rr$ink, font = 2)
  text(L + 1.62, yof(cap * 0.72), "more\ndependence", cex = 0.55, xpd = TRUE,
       adj = 0, col = .rr$soft, font = 3)
  invisible(NULL)
}

#' Plot residual principal-component loadings
#'
#' Residual-component loadings against item location. Opposing clusters suggest
#' a further dimension. Any returned component may be plotted.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param component Which residual principal component to plot (default the
#'   first component).
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_pca(rasch(X))
#' @export
plot_pca <- function(fit, component = 1) {
  component <- .check_whole(component, "component", 1)
  k <- as.integer(component)
  pc <- residual_pca(fit, n_components = k)
  cn <- paste0("PC", k)
  if (!cn %in% names(pc$loadings_matrix))
    stop("component ", k, " is not available (", ncol(pc$loadings_matrix) - 1L,
         " returned)")
  ld <- pc$loadings_matrix[[cn]][match(fit$items$item, pc$loadings_matrix$item)]
  loc <- fit$items$location
  yr <- range(c(0, ld), finite = TRUE)
  if (!all(is.finite(yr)) || diff(yr) <= 0) yr <- c(-0.5, 0.5)
  else yr <- yr + c(-0.08, 0.08) * diff(yr)
  op <- .rr_canvas(range(loc) + c(-0.5, 0.5), yr,
                   "Item location (logits)", paste0("PC", k, " loading"),
                   sprintf("PC%d: eigenvalue %.3f  (%.1f%% of residual variance)",
                           k, pc$eigen_table$eigenvalue[k],
                           100 * pc$eigen_table$proportion[k]), grid_x = TRUE)
  on.exit(par(op))
  abline(h = 0, lty = 2, col = .rr$soft)
  points(loc, ld, pch = 21, cex = 1.6,
         bg = ifelse(ld > 0, .rr$blue, .rr$amber), col = "white", lwd = 1.2)
  text(loc, ld, fit$items$item, pos = 3, offset = 0.45, cex = 0.7, col = .rr$soft)
  invisible(NULL)
}

#' Biplot of the first two residual components
#'
#' Item loadings on the first two residual principal components -- the pair
#' that usually carries any interpretable second dimension -- plotted against
#' one another on equal (isometric) axes. Items far from the origin with
#' opposing signs on PC1 mark a possible contrast, and PC2 separates them
#' further. Point colour follows the sign of the PC1 loading, the split the
#' unidimensionality t-test uses by default.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 8)
#' X <- matrix(rbinom(500 * 8, 1, plogis(outer(rnorm(500), d, "-"))), 500, 8)
#' colnames(X) <- paste0("I", 1:8)
#' plot_pca_biplot(rasch(X))
#' @export
plot_pca_biplot <- function(fit) {
  pc <- residual_pca(fit); lm <- pc$loadings_matrix
  if (!all(c("PC1", "PC2") %in% names(lm))) {
    op <- par(mar = c(0, 0, 0, 0)); on.exit(par(op)); plot.new()
    text(0.5, 0.5, "At least two residual components are needed for a biplot.",
         col = .rr$soft, cex = 0.9)
    return(invisible(NULL))
  }
  x <- lm$PC1; y <- lm$PC2; items <- lm$item
  L <- max(abs(c(x, y)), na.rm = TRUE) * 1.2
  if (!is.finite(L) || L <= 0) L <- 1
  pct <- 100 * pc$eigen_table$proportion
  # equal (asp = 1) axes, so distances between items are read faithfully
  op <- par(mar = c(4.2, 4.4, 1.6, 1.6), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink)
  on.exit(par(op))
  plot(NA, xlim = c(-L, L), ylim = c(-L, L), asp = 1, axes = FALSE,
       xlab = sprintf("PC1 loading  (%.1f%% of residual variance)", pct[1]),
       ylab = sprintf("PC2 loading  (%.1f%%)", pct[2]))
  abline(h = pretty(c(-L, L)), v = pretty(c(-L, L)), col = .rr$grid, lwd = 0.8)
  abline(h = 0, v = 0, lty = 2, col = .rr$soft)
  .rr_axis(1)
  .rr_axis(2)
  points(x, y, pch = 21, cex = 1.6, bg = ifelse(x > 0, .rr$blue, .rr$amber),
         col = "white", lwd = 1.2)
  text(x, y, items, pos = 3, offset = 0.4, cex = 0.7, col = .rr$soft)
  invisible(NULL)
}

#' Plot category frequencies
#'
#' Observed response distribution over the categories of one item.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(1:4, function(i)
#'   sapply(th, function(t)
#'     sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
#' colnames(X) <- sprintf("P%02d", 1:4)
#' plot_catfreq(rasch(X), "P01")
#' @export
plot_catfreq <- function(fit, item) {
  .check_response_display_fit(fit, "category-frequency plots")
  if (length(item) != 1L) stop("`item` must name exactly one item")
  i <- .item_idx(fit, item)
  cnt <- fit$thresholds_diag[[i]]$category_counts
  cats <- seq_along(cnt) - 1L
  op <- .rr_canvas(c(-0.6, max(cats) + 0.6), c(0, max(cnt) * 1.12),
                   "Category", "Count",
                   fit$items$item[i],
                   grid_x = FALSE)
  on.exit(par(op))
  rect(cats - 0.38, 0, cats + 0.38, cnt, col = .rr$blue, border = "white")
  text(cats, cnt, cnt, pos = 3, cex = 0.8, col = .rr$ink)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Person characteristic curve: one person against the item difficulties.
# ---------------------------------------------------------------------------
#' Plot a person characteristic curve
#'
#' The person characteristic curve: the modelled expectation at the
#' person's estimated measure against item location, with the person's
#' observed responses overlaid, grouped into item-difficulty intervals
#' (proportion of maximum score per interval). A wholly dichotomous
#' unit-discrimination fit draws the exact logistic curve. Polytomous,
#' rating scale, many-facet, and frame fits have no single curve in the
#' item location alone, so the model is displayed as its expected
#' proportion of maximum for the actual items within each interval, under
#' the fitted thresholds, response cells, and frame units. Erratic
#' responding (for example lucky guessing on hard items by a
#' low-proficiency person) shows as observed points far from the model,
#' complementing the person fit residual.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param person Row number of the person, or an ID matching
#'   \code{fit$person$id}.
#' @param n_groups Number of item-difficulty intervals for the observed
#'   points (capped by the number of observed items).
#' @param grid Item-location grid over which the dichotomous curve is
#'   drawn; interval displays ignore it.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 12)
#' X <- matrix(rbinom(300 * 12, 1, plogis(outer(rnorm(300), d, "-"))), 300, 12)
#' colnames(X) <- paste0("I", 1:12)
#' plot_pcc(rasch(X), person = 1)
#' @export
plot_pcc <- function(fit, person, n_groups = 5, grid = NULL) {
  n_groups <- .check_whole(n_groups, "n_groups", 2)
  grid <- .model_grid(fit, grid, half_width = 5)
  n <- .person_row(fit, person)
  th <- fit$person$theta[n]
  if (!is.finite(th)) stop("no finite estimate for this person")
  x <- fit$X[n, ]
  loc <- fit$items$location; mm <- fit$m
  ok <- !is.na(x) & is.finite(loc) & is.finite(mm) & mm > 0
  if (sum(ok) < 3) stop("fewer than 3 observed responses for this person")
  # the fitted per-cell expectations carry the thresholds, response cells,
  # and frame units of every model family; plogis(theta - location) is the
  # model only for a wholly dichotomous unit-discrimination fit
  E_n <- fit$moments$E[n, ]
  dich <- all(mm[ok] == 1L) &&
    (is.null(fit$disc) || all(fit$disc[ok] == 1))
  xr <- if (dich) range(grid) else range(loc[ok]) + c(-0.5, 0.5)
  op <- .rr_canvas(xr, c(0, 1),
                   if (inherits(fit, c("rasch_mfrm", "rasch_efrm")))
                     "Response-cell location (logits)" else
                     "Item location (logits)",
                   if (dich) "Probability of success" else
                     "Proportion of maximum score",
                   sprintf("%s  (location %.3f, fit residual %s)",
                           fit$person$id[n], th,
                           ifelse(is.na(fit$person$fit_resid[n]), "NA",
                                  sprintf("%.2f", fit$person$fit_resid[n]))))
  on.exit(par(op))
  abline(v = th, lty = 3, col = .rr$soft)
  g <- .class_intervals(loc[ok], rep(FALSE, sum(ok)),
                        min(n_groups, max(2, floor(sum(ok) / 2))))
  obsL <- tapply(loc[ok], g, mean)
  obsP <- tapply((x[ok] / mm[ok]), g, mean)
  if (dich) {
    lines(grid, plogis(th - grid), lwd = 3, col = .rr$ink)
  } else {
    modP <- tapply((E_n[ok] / mm[ok]), g, mean)
    oo <- order(obsL)
    lines(obsL[oo], modP[oo], lwd = 3, col = .rr$ink)
    points(obsL, modP, pch = 22, bg = .rr$ink, col = "white", cex = 1.3,
           lwd = 1)
  }
  points(obsL, obsP, pch = 21, bg = .rr$blue, col = "white", cex = 1.6, lwd = 1.2)
  .rr_legend("topright", c("Model", "Observed"),
             lwd = c(3, NA), pch = c(if (dich) NA else 22, 21),
             pt.bg = c(if (dich) NA else .rr$ink, .rr$blue),
             col = c(.rr$ink, "white"), pt.cex = 1.4)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Kidmap: the person diagnostic map (Wright, Mead & Ludlow 1980).
# ---------------------------------------------------------------------------
#' Plot a kidmap
#'
#' The person diagnostic map (Wright, Mead and Ludlow 1980): item thresholds
#' the person achieved to the right of a vertical logit axis and thresholds
#' not achieved to the left, with the person's location drawn as a dashed
#' line inside its confidence band. Achieved thresholds above the band
#' (unexpected successes) and unachieved thresholds below it (unexpected
#' failures) are highlighted; a clean response pattern shows achieved
#' thresholds below the band and unachieved ones above it.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param person Row number of the person, or an ID matching
#'   \code{fit$person$id}.
#' @param level Confidence level of the band around the person location used
#'   to mark unexpected responses. The selected person must have a positive,
#'   finite standard error; without it the band and the unexpected-response
#'   classification are unavailable.
#' @param bins Number of vertical bins used to stack the threshold labels.
#' @param xlim Optional logit range; thresholds outside it are omitted.
#' @param cex_labels Character expansion for the threshold labels.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @references Wright, B. D., Mead, R. J., & Ludlow, L. H. (1980).
#'   \emph{KIDMAP: person-by-item interaction mapping} (Research Memorandum
#'   No. 29). Chicago: University of Chicago, MESA Psychometric Laboratory.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 12)
#' X <- matrix(rbinom(300 * 12, 1, plogis(outer(rnorm(300), d, "-"))), 300, 12)
#' colnames(X) <- paste0("I", 1:12)
#' plot_kidmap(rasch(X), person = 1)
#' @export
plot_kidmap <- function(fit, person, level = 0.95, bins = 35, xlim = NULL,
                        cex_labels = 0.8) {
  .check_response_display_fit(fit, "person diagnostic maps")
  .check_prob(level, "level")
  bins <- .check_whole(bins, "bins", 2)
  .check_xlim(xlim)
  .check_cex(cex_labels)
  n <- .person_row(fit, person)
  th <- fit$person$theta[n]
  if (!is.finite(th)) stop("no finite estimate for this person")
  se <- fit$person$se[n]
  if (!is.finite(se) || se <= 0)
    stop("no positive finite standard error for this person; the kidmap confidence band is unavailable",
         call. = FALSE)
  z <- stats::qnorm(1 - (1 - level) / 2)
  band <- z * se
  x <- fit$X[n, ]
  thr <- fit$thresholds
  obs <- !is.na(x[thr$item])
  ii <- thr$item[obs]; kk <- thr$k[obs]; tv <- thr$tau[obs]
  att <- x[ii] >= kk
  lab <- if (all(fit$m == 1L)) fit$items$item[ii] else
    paste0(fit$items$item[ii], ".", kk)
  rng <- if (is.null(xlim))
    range(c(th - band, th + band, tv)) + c(-0.5, 0.5) else sort(xlim)
  keep <- tv >= rng[1] & tv <= rng[2]
  tv <- tv[keep]; att <- att[keep]; lab <- lab[keep]
  .check_map_range(length(tv))
  unexp <- (att & tv > th + band) | (!att & tv < th - band)
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  bin <- pmin(findInterval(tv, brk, rightmost.closed = TRUE), bins)
  split <- 0.5
  op <- par(mar = c(2.6, 4.4, 1.6, 0.5), mgp = c(2.5, 0.7, 0), tcl = -0.25,
            las = 1, col.axis = .rr$ink, col.lab = .rr$ink, col.main = .rr$ink,
            font.main = 2, cex.main = 1.15)
  on.exit(par(op))
  plot(NA, xlim = c(0, 1), ylim = rng, xlab = "", ylab = "Location (logits)",
       axes = FALSE, main = "")
  abline(h = .rr_ticks(rng), col = .rr$grid, lwd = 0.8)
  .rr_axis(2)
  if (band > 0)
    rect(0, th - band, 1, th + band, col = adjustcolor(.rr$blue, 0.10),
         border = NA)
  segments(split, rng[1], split, rng[2], col = .rr$ink, lwd = 1)
  segments(0, th, 1, th, col = .rr$blue, lty = 2, lwd = 1.4)
  cexl <- cex_labels
  colw <- max(strwidth(lab, cex = cexl)) * 1.2
  for (side in c(TRUE, FALSE)) {           # TRUE = achieved (right)
    sel <- att == side
    for (b in unique(bin[sel])) {
      pick <- sel & bin == b
      ls <- lab[pick][order(tv[pick])]
      cols <- ifelse(unexp[pick][order(tv[pick])], .rr$red, .rr$ink)
      kmax <- max(1L, floor((split - 0.015) / colw))
      if (length(ls) > kmax) {
        cols <- c(cols[seq_len(kmax - 1L)], .rr$soft)
        ls <- c(ls[seq_len(kmax - 1L)], paste0("+", length(ls) - kmax + 1L))
      }
      xs <- if (side) split + 0.015 + (seq_along(ls) - 1L) * colw
            else split - 0.015 - (seq_along(ls) - 1L) * colw
      text(xs, (brk[b] + brk[b + 1L]) / 2, ls, cex = cexl,
           adj = if (side) 0 else 1, col = cols, font = 2)
    }
  }
  fr <- fit$person$fit_resid[n]
  text(0.01, rng[2], sprintf("%s: location %.2f (SE %s), fit residual %s",
                             fit$person$id[n], th,
                             ifelse(is.na(se), "NA", sprintf("%.2f", se)),
                             ifelse(is.na(fr), "NA", sprintf("%.2f", fr))),
       cex = 0.78, adj = c(0, 1), col = .rr$ink, font = 2)
  mtext("not achieved", side = 1, line = 0.4, at = split / 2, cex = 0.8,
        col = .rr$soft)
  mtext("achieved", side = 1, line = 0.4, at = (1 + split) / 2, cex = 0.8,
        col = .rr$soft)
  if (any(unexp))
    legend("bottomleft", "unexpected response", bty = "n",
           text.col = .rr$red, cex = 0.85)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Residual statistics distributions of the fit residuals.
# ---------------------------------------------------------------------------
#' Plot the fit residual distribution
#'
#' A histogram of the
#' item or person fit residuals -- the log-transformed statistic or its
#' untransformed natural form -- against the standard normal density they
#' should approximate under fit (Andrich and Marais 2019, ch. 15). The
#' natural residual is visibly skewed
#' (that is why the log transform is reported); both are available.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param what \code{"items"} or \code{"persons"}.
#' @param statistic \code{"fit_resid"} (log-transformed, default) or
#'   \code{"natural"}.
#' @param bins Number of histogram bins.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 10)
#' X <- matrix(rbinom(400 * 10, 1, plogis(outer(rnorm(400), d, "-"))), 400, 10)
#' colnames(X) <- paste0("I", 1:10)
#' plot_resid_dist(rasch(X), what = "persons")
#' @export
plot_resid_dist <- function(fit, what = c("items", "persons"),
                            statistic = c("fit_resid", "natural"), bins = 25) {
  .check_response_display_fit(fit, "fit-residual plots")
  what <- match.arg(what); statistic <- match.arg(statistic)
  bins <- .check_whole(bins, "bins", 2)
  v <- if (what == "items") {
    if (statistic == "fit_resid") fit$items$fit_resid else fit$items$natural_resid
  } else {
    if (statistic == "fit_resid") fit$person$fit_resid else fit$person$natural_resid
  }
  v <- v[is.finite(v)]
  if (length(v) < 3) stop("fewer than 3 residuals to display")
  lab <- if (statistic == "fit_resid") "fit residual (log-transformed)"
         else "natural fit residual"
  rng <- range(c(v, -3, 3)); rng <- rng + c(-0.5, 0.5) * diff(rng) * 0.05
  brk <- seq(rng[1], rng[2], length.out = bins + 1)
  h <- hist(v, breaks = brk, plot = FALSE)
  ymax <- max(h$density, dnorm(0)) * 1.15
  op <- .rr_canvas(rng, c(0, ymax), lab, "Density",
                   sprintf("n = %d, mean %.2f, SD %.2f", length(v), mean(v),
                           sd(v)), grid_x = FALSE)
  on.exit(par(op))
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density,
       col = adjustcolor(.rr$blue, 0.45), border = "white")
  xs <- seq(rng[1], rng[2], length.out = 200)
  lines(xs, dnorm(xs), lwd = 2.6, col = .rr$red)
  abline(v = c(-2.5, 2.5), lty = 3, col = .rr$soft)
  .rr_legend("topright", c("Observed", "Standard normal"),
             lwd = c(NA, 2.6), pch = c(22, NA), pt.bg = c(adjustcolor(.rr$blue, 0.45), NA),
             col = c("white", .rr$red), pt.cex = 1.6)
  invisible(NULL)
}
