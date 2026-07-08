# Plot a test-equating comparison

Scatter of the two calibrations' common-item locations with the shifted
identity line and per-item 95 per cent bands; drifting items
(BH-adjusted) are highlighted and labelled.

## Usage

``` r
plot_equate(fit, reference, shift = c("mean", "none"))
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- reference:

  A second
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  fit, or a data frame with columns `item`, `location`, and optionally
  `se`.

- shift:

  Passed to
  [`equate_tests`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md).

## Value

Called for its plotting side effect; invisibly the
[`equate_tests`](https://drjoshmcgrane.github.io/rasch/reference/equate_tests.md)
result.

## Examples

``` r
set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
mk <- function() {
  X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
  colnames(X) <- paste0("I", 1:8); rasch(X)
}
plot_equate(mk(), mk())
```
