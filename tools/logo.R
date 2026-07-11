#!/usr/bin/env Rscript
# tools/logo.R
# ---------------------------------------------------------------------------
# Reproducible hex-sticker logo for the 'rasch' package.
#
# Draws, with base graphics only (no external design tools, no font files):
#   * a pointy-top regular hexagon (community-standard 2:sqrt(3) aspect),
#   * one clean logistic item characteristic curve (ICC / ogive) in blue,
#   * a dashed red reference marking the item location where p = 0.5,
#     dropped to the theta axis,
#   * a bold lowercase "rasch" wordmark.
#
# Writes:  man/figures/logo.svg   (vector, via svg() / cairo)
#          man/figures/logo.png   (>= 1024 px tall, transparent background)
#
# Run from the package root:  Rscript tools/logo.R
# ---------------------------------------------------------------------------

## ---- palette (from the pkgdown theme) ------------------------------------
navy  <- "#1d3557"  # primary  -- hexagon fill
ink   <- "#1e293b"  # ink
blue  <- "#2563eb"  # blue accent
red   <- "#dc2626"  # red accent -- reference line
light <- "#f8fafc"  # very light slate
white <- "#ffffff"

## Theme: navy-filled hexagon with light artwork.
hex_fill   <- navy
hex_border <- "#4d7cc7"   # a brighter steel-blue rim for crisp definition
ogive_col  <- "#7fb2ff"   # luminous blue so the ogive reads on navy
axis_col   <- "#8aa2c4"   # muted slate-blue axes
ref_col    <- "#f26d6d"   # brightened red so the reference reads on navy
word_col   <- white       # wordmark

## ---- geometry ------------------------------------------------------------
R     <- 1                         # circumradius (centre -> vertex)
halfW <- sqrt(3) / 2 * R           # 0.8660254 : flat-to-flat half width
pad   <- 1.045                     # margin so the border stroke isn't clipped
xlim  <- c(-halfW, halfW) * pad
ylim  <- c(-R, R) * pad

# Pointy-top hexagon vertices (a point at top and bottom).
hex_xy <- function(r = R) {
  ang <- c(90, 150, 210, 270, 330, 30) * pi / 180
  list(x = r * cos(ang), y = r * sin(ang))
}

## ---- item characteristic curve -------------------------------------------
# Logistic ogive p = plogis(theta), mapped into a chart box inside the hex.
theta  <- seq(-4.4, 4.4, length.out = 512)
p      <- plogis(theta)

cx0    <- -0.50;  cx1 <- 0.50       # curve horizontal extent (centred on 0)
axisY  <- -0.05                     # baseline, i.e. p = 0
curveH <- 0.64                      # vertical span for p in [0, 1]

px <- cx0 + (theta - min(theta)) / (max(theta) - min(theta)) * (cx1 - cx0)
py <- axisY + p * curveH

# Item location b: theta where p = 0.5 -> mid of the mapped x-range.
bx <- (cx0 + cx1) / 2               # = 0
by <- axisY + 0.5 * curveH          # curve height at p = 0.5

## ---- one drawing routine, used for every device --------------------------
draw_sticker <- function() {
  op <- par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i",
            lend = "round", ljoin = "round")
  on.exit(par(op))

  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1)

  ## hexagon: fill first, then a crisp border on top
  h <- hex_xy(R)
  polygon(h$x, h$y, col = hex_fill, border = NA)
  polygon(h$x, h$y, col = NA, border = hex_border, lwd = 6.5)

  ## axes (subtle L-shape so the curve reads as a plot)
  segments(cx0 - 0.02, axisY, cx1 + 0.04, axisY, col = axis_col, lwd = 2)
  segments(cx0, axisY, cx0, axisY + curveH + 0.03, col = axis_col, lwd = 2)

  ## p = 0.5 reference: horizontal to the curve, then dropped to the axis
  segments(cx0, by, bx, by, col = ref_col, lwd = 3.2, lty = "22")
  segments(bx, by, bx, axisY, col = ref_col, lwd = 3.2, lty = "22")
  points(bx, by, pch = 19, cex = 1.6, col = ref_col)

  ## the ogive
  lines(px, py, col = ogive_col, lwd = 9)

  ## wordmark
  text(0, -0.50, "rasch", col = word_col, font = 2,
       family = "sans", cex = 6.2)

  invisible(NULL)
}

## ---- write outputs -------------------------------------------------------
out_dir <- "man/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Device aspect matches the hexagon's bounding box exactly (sqrt(3) : 2),
# so a regular hexagon is drawn without distortion.
asp_wh <- (halfW * pad) / (R * pad)   # width / height = sqrt(3)/2

## SVG (vector)
svg(file.path(out_dir, "logo.svg"),
    width = 6.06, height = 6.06 / asp_wh, bg = "transparent")
draw_sticker()
dev.off()

## PNG (>= 1024 px tall, transparent background)
png(file.path(out_dir, "logo.png"),
    width = round(1400 * asp_wh), height = 1400, units = "px",
    res = 200, bg = "transparent", type = "cairo")
draw_sticker()
dev.off()

message("Wrote ", file.path(out_dir, "logo.svg"), " and ",
        file.path(out_dir, "logo.png"),
        "  (", round(1400 * asp_wh), " x ", 1400, " px)")
