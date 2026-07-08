# Spread-parameter test for dependence within subtests

Andrich's (1985) least-upper-bound screen: the spread component
\\\lambda\\ of a polytomous item (half the distance between successive
thresholds in the principal-components parameterisation, estimated here
by
[`pcml_pc`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md))
cannot fall below the value implied by the binomial distribution when
the item is a subtest of equally difficult, independent dichotomous
items; different difficulties only raise it. A spread estimate below the
bound therefore indicates response dependence among the members (Andrich
and Marais 2019, Table 24.1). Typically applied after
[`combine_items`](https://drjoshmcgrane.github.io/rasch/reference/combine_items.md),
whose super-items are exactly such subtests.

## Usage

``` r
spread_test(fit, maxit = 60, tol = 1e-08)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- maxit, tol:

  Passed to the
  [`pcml_pc`](https://drjoshmcgrane.github.io/rasch/reference/pcml_pc.md)
  refit.

## Value

A data frame with one row per polytomous item: `item`, `m`, the `spread`
estimate and its `se`, the bound `lub` (available for maximum scores 2
to 8), `z` = (spread - lub)/se, and `dependent` = spread below the
bound. Dichotomous items carry no spread and are omitted.

## Examples

``` r
set.seed(1); N <- 600
d0 <- seq(-1.5, 1.5, length.out = 8)
X <- matrix(rbinom(N * 8, 1, plogis(outer(rnorm(N), d0, "-"))), N, 8)
X[, 5] <- X[, 4]; X[, 6] <- X[, 4]                 # a dependent triple
colnames(X) <- paste0("I", 1:8)
fit2 <- combine_items(rasch(X), list(c("I4", "I5", "I6"), c("I1", "I2", "I3")))
spread_test(fit2)
#> Spread-parameter screen (Andrich 1985): spread below the binomial bound indicates dependence
#>      item m spread    se bound     z dependent
#>  I1+I2+I3 3  0.653 0.114 0.550 0.899          
```
