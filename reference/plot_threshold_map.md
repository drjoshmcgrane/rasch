# Plot the threshold map

Each item's threshold locations on a common logit scale, ordered by item
location, with disordered thresholds highlighted.

## Usage

``` r
plot_threshold_map(fit, order_by_location = TRUE)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- order_by_location:

  Order items by their location (the default) rather than their original
  sequence.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_threshold_map(rasch(X))
```
