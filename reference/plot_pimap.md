# Plot the person-item threshold distribution

The targeting display: the person location distribution above the axis
and the item threshold distribution mirrored below it, on a shared logit
scale.

## Usage

``` r
plot_pimap(fit, bins = 35, xlim = NULL)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- bins:

  Number of histogram bins.

- xlim:

  Optional logit range for the shared scale; persons and thresholds
  outside it are omitted.

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
