# Plot residual principal-component loadings

Residual-component loadings against item location; opposing clusters at
top and bottom suggest a further dimension. Any leading component may be
shown, not only the first.

## Usage

``` r
plot_pca(fit, component = 1)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- component:

  Which residual principal component to plot (default the first
  component).

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_pca(rasch(X))
```
