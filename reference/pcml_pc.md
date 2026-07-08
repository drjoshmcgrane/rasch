# Estimate Rasch thresholds via the Andrich principal-components reparameterisation

An optional alternative to
[`pcml`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md)'s
free-threshold estimation, useful when some response categories are
sparsely populated. Each item's thresholds are re-expressed as up to
four orthogonal polynomial components in the category score: location,
spread, skewness, and kurtosis (Andrich 1978, 1985; Pedler 1987).
Location is always estimated; spread, skewness, and kurtosis are added
in turn as an item's number of thresholds and `n_components` allow.
Estimation uses the same pairwise conditional likelihood as `pcml`
(Andrich and Luo 2003), so it inherits the same missing-data handling
and sandwich standard errors. The component family stops at the quartic
(kurtosis) term, so the reparameterisation is exact, matching `pcml`'s
free partial credit thresholds and log-likelihood, only while every item
has at most 3 thresholds (4 categories); from 4 thresholds on `pcml_pc`
is necessarily a reduced-rank smoothing of the thresholds to a
polynomial trend across categories, however large `n_components` is set,
trading flexibility for the stability that comes from pooling
information across all of an item's categories – useful when a category
has low or zero frequency.

## Usage

``` r
pcml_pc(X, n_components = 4, maxit = 60, tol = 1e-08)
```

## Arguments

- X:

  Persons-by-items integer score matrix (categories from 0). Missing
  values are handled by pairwise deletion.

- n_components:

  Maximum number of components per item: 1 (location only) up to 4
  (location, spread, skewness, kurtosis; the highest derived by Pedler
  1987). Capped per item at its own number of thresholds, and further
  wherever a component would be collinear with lower-order ones for that
  item's threshold count (kurtosis is unidentified, and dropped, at
  exactly 4 thresholds).

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

A list with the threshold table `thr` (columns `id`, `item`, `k`, `tau`,
`se`), the component table `components` (one row per item, with
`location`, `spread`, `skewness`, `kurtosis` and their standard errors,
`NA` where an item's rank does not support that component), the
threshold covariance matrix `cov_tau`, the pairwise conditional
log-likelihood, the iteration count, a convergence flag, and the
max-score vector `m`.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
pcml_pc(X)$components
#>   item   location location_se spread spread_se skewness skewness_se kurtosis
#> 1   I1 -1.5174835   0.1297013     NA        NA       NA          NA       NA
#> 2   I2 -0.7730904   0.1121427     NA        NA       NA          NA       NA
#> 3   I3 -0.2093836   0.1023242     NA        NA       NA          NA       NA
#> 4   I4  0.3667512   0.1047318     NA        NA       NA          NA       NA
#> 5   I5  0.7889434   0.1105680     NA        NA       NA          NA       NA
#> 6   I6  1.3442630   0.1221786     NA        NA       NA          NA       NA
#>   kurtosis_se
#> 1          NA
#> 2          NA
#> 3          NA
#> 4          NA
#> 5          NA
#> 6          NA
```
