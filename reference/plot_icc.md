# Plot an item characteristic curve

Draws the model expected-score curve with observed class-interval means
overlaid. With `group` supplied, observed means are drawn separately per
group, the conventional graphical DIF display.

## Usage

``` r
plot_icc(fit, item, group = NULL, n_groups = NULL, grid = seq(-5, 5, 0.05))
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- item:

  Item name or column index.

- group:

  Optional person grouping vector, or one or more names of factors
  nominated in the fit, for a DIF overlay; several names give the
  factor-combination cells (the factorial display).

- n_groups:

  Number of class intervals for the observed means; by default the fit's
  own count, or – with a `group` overlay – a count adapted to keep the
  smallest group's interval cells adequately filled.

- grid:

  Logit grid over which to draw the model curve.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- sprintf("I%02d", 1:6)
plot_icc(rasch(X), "I03")
```
