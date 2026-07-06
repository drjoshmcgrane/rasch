# Estimate Rasch thresholds by pairwise conditional maximum likelihood

Maximises the pairwise conditional likelihood, in which the person
parameter cancels within every item pair, by Newton-Raphson (Andrich and
Luo 2003; Zwinderman 1995). The partial credit model estimates every
threshold freely; the rating scale model constrains
`tau_ik = delta_i + kappa_k` through the design matrix.

## Usage

``` r
pcml(X, model = c("PCM", "RSM"), anchors = NULL, maxit = 60, tol = 1e-08)
```

## Arguments

- X:

  Persons-by-items integer score matrix (categories from 0). Missing
  values are handled by pairwise deletion.

- model:

  `"PCM"` or `"RSM"`.

- anchors:

  Optional anchor table for equating: a data frame with columns `item`
  (name or column index), `k`, and `tau` (the fixed value). A numeric
  `k` fixes that single threshold (individual anchoring); `k = NA` fixes
  the item's mean location at `tau` while its thresholds remain free
  (average anchoring). The remaining parameters are estimated on the
  anchored scale and no recentring is applied. PCM only.

- maxit, tol:

  Newton-Raphson iteration cap and convergence tolerance.

## Value

A list with the threshold table `thr` (columns `id`, `item`, `k`, `tau`,
`se`, `anchored`), the threshold covariance matrix `cov_tau`, the
pairwise conditional log-likelihood, the iteration count, a convergence
flag, and the max-score vector `m`.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(400 * 6, 1, plogis(outer(rnorm(400), d, "-"))), 400, 6)
colnames(X) <- paste0("I", 1:6)
pcml(X)$thr
#>   id item k        tau        se anchored
#> 1  1    1 1 -1.5174835 0.1297013    FALSE
#> 2  2    2 1 -0.7730904 0.1121427    FALSE
#> 3  3    3 1 -0.2093836 0.1023242    FALSE
#> 4  4    4 1  0.3667512 0.1047318    FALSE
#> 5  5    5 1  0.7889434 0.1105680    FALSE
#> 6  6    6 1  1.3442630 0.1221786    FALSE
# anchor two items at fixed values (equating)
pcml(X, anchors = data.frame(item = c("I1", "I6"), k = 1, tau = c(-1.5, 1.5)))$thr
#>   id item k        tau        se anchored
#> 1  1    1 1 -1.5000000 0.0000000     TRUE
#> 2  2    2 1 -0.6880959 0.1613982    FALSE
#> 3  3    3 1 -0.1211423 0.1463167    FALSE
#> 4  4    4 1  0.4586830 0.1532659    FALSE
#> 5  5    5 1  0.8835376 0.1568160    FALSE
#> 6  6    6 1  1.5000000 0.0000000     TRUE
```
