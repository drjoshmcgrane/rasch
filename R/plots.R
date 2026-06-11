# RaschR plots :: Rasch diagnostic plots in base R graphics (no dependencies).
# House palette: navy / ochre / terracotta.
.NAVY <- "#1f3b5c"; .OCHRE <- "#d8a13a"; .TERRA <- "#c1622d"; .GREY <- "#8a8f98"
.PAL  <- c("#1f3b5c", "#c1622d", "#d8a13a", "#4f7d68", "#7a5c8e", "#9c4a4a")

.item_idx <- function(fit, item) if (is.character(item)) match(item, fit$items$item) else item

# Item characteristic curve: model expected score with observed class-interval means.
#' Plot an item characteristic curve
#'
#' Draws the model expected-score curve with the observed class-interval means
#' overlaid.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @param n_groups Number of class intervals for the observed means.
#' @param grid Logit grid over which to draw the model curve.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- sprintf("I%02d", 1:6)
#' plot_icc(rasch(X), "I03")
#' @export
plot_icc <- function(fit, item, n_groups = fit$n_groups, grid = seq(-5, 5, 0.05)) {
  i <- .item_idx(fit, item); tau_i <- fit$tau_list[[i]]; mmax <- length(tau_i)
  Ecurve <- vapply(grid, function(th) item_moments(th, tau_i)$E, 0)
  th <- fit$theta_person; x <- fit$X[, i]; ok <- !is.na(th) & !is.na(x)
  ci <- cut(rank(th[ok], ties.method = "first"), n_groups, labels = FALSE)
  obsTh <- tapply(th[ok], ci, mean); obsX <- tapply(x[ok], ci, mean)
  plot(grid, Ecurve, type = "l", lwd = 2.5, col = .NAVY, ylim = c(0, mmax),
       xlab = "Person location (logits)", ylab = "Expected score",
       main = paste0("ICC: ", fit$items$item[i]), bty = "n")
  abline(h = seq(0, mmax, 0.5), col = "#eeeeee")
  lines(grid, Ecurve, lwd = 2.5, col = .NAVY)
  points(obsTh, obsX, pch = 19, col = .TERRA, cex = 1.2)
  legend("topleft", c("Model", "Observed (class intervals)"), lwd = c(2.5, NA),
         pch = c(NA, 19), col = c(.NAVY, .TERRA), bty = "n")
}

# Category probability curves with threshold markers.
#' Plot category probability curves
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @param item Item name or column index.
#' @param grid Logit grid over which to draw the curves.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' simP <- function(th, t) { x <- 0:length(t); p <- exp(x * th - c(0, cumsum(t))); p / sum(p) }
#' th <- rnorm(400)
#' X <- sapply(1:4, function(i) sapply(th, function(t) sample(0:3, 1, prob = simP(t, c(-1, 0, 1)))))
#' colnames(X) <- sprintf("P%02d", 1:4)
#' plot_ccc(rasch(X), "P01")
#' @export
plot_ccc <- function(fit, item, grid = seq(-6, 6, 0.05)) {
  i <- .item_idx(fit, item); tau_i <- fit$tau_list[[i]]; mmax <- length(tau_i)
  P <- vapply(grid, function(th) item_moments(th, tau_i)$P, numeric(mmax + 1))
  cols <- rep(.PAL, length.out = mmax + 1)
  plot(NA, xlim = range(grid), ylim = c(0, 1), xlab = "Person location (logits)",
       ylab = "Category probability", main = paste0("Category curves: ", fit$items$item[i]), bty = "n")
  for (cat in 0:mmax) lines(grid, P[cat + 1, ], lwd = 2.2, col = cols[cat + 1])
  abline(v = tau_i, lty = 3, col = .GREY)
  ord <- if (all(diff(tau_i) > 0)) "thresholds ordered" else "THRESHOLDS DISORDERED"
  mtext(ord, side = 3, line = 0.2, cex = 0.8,
        col = if (all(diff(tau_i) > 0)) .NAVY else .TERRA)
  legend("topright", legend = paste0("cat ", 0:mmax), lwd = 2.2, col = cols, bty = "n", cex = 0.8)
}

# Person-item threshold map (Wright map): person distribution over item thresholds.
#' Plot the person-item threshold map
#'
#' Shows the person location distribution above the item threshold locations on
#' a shared logit scale, with disordered items highlighted.
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_pimap(rasch(X))
#' @export
plot_pimap <- function(fit) {
  th <- fit$theta_person[!is.na(fit$theta_person)]
  thr <- fit$thresholds; rng <- range(c(th, thr$tau), na.rm = TRUE)
  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(1:2, 2, 1), heights = c(1, 1.5)); par(mar = c(0.5, 4, 2, 2))
  h <- hist(th, breaks = 30, plot = FALSE)
  barplot(h$counts, space = 0, col = .NAVY, border = NA, axes = FALSE,
          main = "Person-item threshold map", xlim = c(0, length(h$counts)))
  mtext("Persons", side = 2, line = 2, cex = 0.9)
  par(mar = c(4, 4, 0.5, 2))
  plot(NA, xlim = rng, ylim = c(0.5, length(fit$tau_list) + 0.5),
       xlab = "Location (logits)", ylab = "", yaxt = "n", bty = "n")
  axis(2, at = seq_along(fit$tau_list), labels = fit$items$item, las = 1, cex.axis = 0.7)
  for (i in seq_along(fit$tau_list)) {
    tau_i <- fit$tau_list[[i]]
    disord <- !all(diff(tau_i) > 0)
    points(tau_i, rep(i, length(tau_i)), pch = 19, cex = 1.1,
           col = if (disord) .TERRA else .OCHRE)
    points(mean(tau_i), i, pch = 4, col = .NAVY, lwd = 2)   # item location
  }
  legend("topright", c("threshold", "disordered", "item location"),
         pch = c(19, 19, 4), col = c(.OCHRE, .TERRA, .NAVY), bty = "n", cex = 0.75)
}

# Residual-correlation heatmap (blues to reds, per house preference).
#' Plot the residual-correlation heatmap
#'
#' @param fit A fitted object from \code{\link{rasch}}.
#' @return Called for its plotting side effect; invisibly \code{NULL}.
#' @examples
#' set.seed(1)
#' d <- seq(-2, 2, length.out = 6)
#' X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
#' colnames(X) <- paste0("I", 1:6)
#' plot_resid_cor(rasch(X))
#' @export
plot_resid_cor <- function(fit) {
  R <- residual_correlations(fit)$matrix; L <- ncol(R)
  pal <- colorRampPalette(c("#2c5f8a", "#dfe6ec", "#b2342c"))(64)
  op <- par(no.readonly = TRUE); on.exit(par(op)); par(mar = c(5, 5, 3, 2))
  image(1:L, 1:L, R[, L:1], col = pal, zlim = c(-1, 1), axes = FALSE,
        xlab = "", ylab = "", main = "Residual correlations")
  axis(1, 1:L, colnames(R), las = 2, cex.axis = 0.6)
  axis(2, 1:L, rev(colnames(R)), las = 1, cex.axis = 0.6)
}
