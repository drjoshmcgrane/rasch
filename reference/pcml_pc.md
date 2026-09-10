# Estimate Rasch thresholds using a principal-component parameterisation

Re-expresses each item's thresholds as orthogonal polynomial components:
location, spread, skewness, and kurtosis (Andrich and Luo 2003).
Estimation uses the same pairwise conditional likelihood as
[`pcml`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md). With
at most four thresholds per item the full parameterisation is exact.
Items with five or more thresholds are fitted by a reduced-rank
polynomial trend, which can stabilise sparse categories at the cost of
restricting the threshold pattern.

## Usage

``` r
pcml_pc(X, n_components = 4, maxit = 60, tol = 1e-08)
```

## Arguments

- X:

  Persons-by-items integer score matrix. Each item must have at least
  two observed categories, numbered consecutively from 0. Missing values
  are handled by pairwise deletion.

- n_components:

  Maximum number of components per item: 1 (location only) up to 4
  (location, spread, skewness, kurtosis). Capped per item at its own
  number of thresholds. Kurtosis requires at least four thresholds (five
  response categories).

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

A list with the threshold table `thr` (columns `id`, `item`, `k`, `tau`,
`se`), the component table `components` (one row per item, with
`location`, `spread`, `skewness`, `kurtosis` and their standard errors,
`NA` where an item's rank does not support that component), the
threshold covariance matrix `cov_tau`, the pairwise conditional
log-likelihood, the iteration count, a convergence flag, and the
max-score vector `m`. If estimation does not converge, the function
warns and all standard errors and covariance entries are `NA`. The
independent-person and effective-support conditions described for
[`pcml`](https://drjoshmcgrane.github.io/rasch/reference/pcml.md) also
apply.

## Details

For scores \\x=0,\ldots,m\\, the cumulative threshold function is
\$\$C(x)=x\omega_1-x(m-x)\omega_2-x(m-x)(2x-m)\omega_3
-x(m-x)(5x^2-5mx+m^2+1)\omega_4.\$\$ Threshold \\k\\ is \\C(k)-C(k-1)\\.
The four coefficients are location, spread, skewness and kurtosis,
respectively.

## References

Andrich, D. and Luo, G. (2003). Conditional pairwise estimation in the
Rasch model for ordered response categories using principal components.
Journal of Applied Measurement, 4(3), 205–221.

Zwinderman, A. H. (1995). Pairwise parameter estimation in Rasch models.
Applied Psychological Measurement, 19(4), 369–375.

Andrich, D. (1978). A rating formulation for ordered response
categories. Psychometrika, 43(4), 561–573.

Andrich, D. (1985). An elaboration of Guttman scaling with Rasch models
for measurement. In N. B. Tuma (Ed.), Sociological Methodology 1985 (pp.
33–80). Jossey-Bass.

Pedler, P. J. (1987). Accounting for psychometric dependence with a
class of latent trait models. PhD thesis, University of Western
Australia.

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
