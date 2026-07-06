# Plot the Guttman scalogram

Displays the Guttman-ordered response matrix as a heatmap (dark for high
categories), with persons sorted by location down the rows and items by
location across the columns. The coefficient of reproducibility is shown
in the subtitle.

## Usage

``` r
plot_guttman(fit, max_persons = 80)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- max_persons:

  Persons are thinned to at most this many evenly spaced rows for
  legibility on large samples.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), d, "-"))), 200, 6)
colnames(X) <- paste0("I", 1:6)
plot_guttman(rasch(X))
```
