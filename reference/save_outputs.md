# Save every output of a Rasch analysis to a folder

Writes all tables (item statistics, thresholds with standard errors,
person estimates including ID and factors, the score-to-measure table,
residual correlations, principal-component loadings, category
frequencies, and DIF results for every nominated factor) as CSV; every
plot, including the per-item characteristic, category, threshold, and
frequency plots, as PNG and optionally PDF; and a plain-text analysis
summary.

## Usage

``` r
save_outputs(
  fit,
  dir,
  formats = c("png", "pdf"),
  width = 9,
  height = 6,
  dpi = 300,
  item_plots = TRUE
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- dir:

  Output directory; created if absent.

- formats:

  Plot formats, any of `"png"` and `"pdf"`.

- width, height:

  Plot size in inches.

- dpi:

  PNG resolution; the default 300 is publication quality.

- item_plots:

  Also write the per-item plot set (one ICC, category curve, threshold
  curve, and frequency chart per item).

## Value

Invisibly, the vector of files written.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
out <- file.path(tempdir(), "rasch-out")
save_outputs(rasch(X), out, formats = "png", item_plots = FALSE)
```
