# Save a plot for every item

Writes one plot per item – the item characteristic curve, category
probability curves, threshold probability curves, or category
frequencies – to a single multi-page PDF or a ZIP archive of PNGs,
chosen by the extension of `file`.

## Usage

``` r
save_item_plots(
  fit,
  what = c("icc", "ccc", "tpc", "cfreq"),
  file,
  items = NULL,
  n_groups = fit$n_groups,
  grid = seq(-5, 5, 0.05),
  observed = TRUE,
  width = 8,
  height = 5.5,
  dpi = 300
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- what:

  Which plot: `"icc"`, `"ccc"`, `"tpc"`, or `"cfreq"`.

- file:

  Output path ending in `.pdf` (one page per item) or `.zip` (one PNG
  per item).

- items:

  Item names or indices; all items by default.

- n_groups:

  Class intervals for observed overlays.

- grid:

  Logit grid for the curves.

- observed:

  Overlay observed proportions on the category and threshold probability
  curves.

- width, height, dpi:

  Device size in inches and PNG resolution.

## Value

Invisibly, the output path.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
f <- rasch(X)
save_item_plots(f, "icc", file.path(tempdir(), "icc_all.pdf"))
```
