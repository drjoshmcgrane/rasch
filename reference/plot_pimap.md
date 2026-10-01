# Plot the person-item threshold distribution

The targeting display: the person location distribution above the axis
and the calibration threshold distribution mirrored below it, on a
shared logit scale. Dashed lines mark the person and threshold means in
their distributions' colours. MFRM and EFRM thresholds belong to
response cells.

## Usage

``` r
plot_pimap(
  fit,
  bins = 35,
  xlim = NULL,
  information = FALSE,
  group = NULL,
  items = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- bins:

  Number of histogram bins.

- xlim:

  Optional logit range for the shared scale; persons and thresholds
  outside it are omitted. By default the range is extended to labelled
  tick marks beyond the most extreme plotted estimate.

- information:

  Whether to overlay the test information function on a separate
  right-hand axis. Fits with more than one administrable design receive
  one curve per design, a fit with split items one per test form.

- group:

  Optional person-group level: one level of a fitted person factor,
  restricting the person distribution to those persons. A level no
  fitted factor carries is an error.

- items:

  Optional item selection restricting the threshold distribution: item
  names, or one item-set name of an extended-frame fit, whose virtual
  item-by-group cells match through their underlying items. The
  selection is named in the legend, so a restricted map cannot be read
  as the whole instrument. The information curve follows the same item
  selection; for an extended-frame fit it also follows the response
  cells occupied by the selected person group. For EFRM and MFRM, only
  response patterns present in that person group are shown, and for a
  fit with split items only the forms its persons answered.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_pimap(rasch(X))
```
