# Plot the fit residual distribution

A histogram of the item or person fit residuals – the log-transformed
statistic or its untransformed natural form – against the standard
normal density they should approximate under fit (Andrich and Marais
2019, ch. 15). The natural residual is visibly skewed (that is why the
log transform is reported); both are available.

## Usage

``` r
plot_resid_dist(
  fit,
  what = c("items", "persons"),
  statistic = c("fit_resid", "natural"),
  bins = 25
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- what:

  `"items"` or `"persons"`.

- statistic:

  `"fit_resid"` (log-transformed, default) or `"natural"`.

- bins:

  Number of histogram bins.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 10)
X <- matrix(rbinom(400 * 10, 1, plogis(outer(rnorm(400), d, "-"))), 400, 10)
colnames(X) <- paste0("I", 1:10)
plot_resid_dist(rasch(X), what = "persons")
```
