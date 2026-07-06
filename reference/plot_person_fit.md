# Plot person fit

Person locations against person fit residuals with the +/- 2.5 band;
persons beyond the band respond erratically (positive) or too
deterministically (negative).

## Usage

``` r
plot_person_fit(fit, band = 2.5)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- band:

  Fit residual acceptance band.

## Value

Called for its plotting side effect; invisibly `NULL`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
plot_person_fit(rasch(X))
```
