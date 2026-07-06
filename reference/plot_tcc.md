# Plot the test characteristic curve

Expected total score against person location for the whole instrument.

## Usage

``` r
plot_tcc(fit, grid = seq(-6, 6, 0.05))
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- grid:

  Logit grid.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_tcc(rasch(X))
```
