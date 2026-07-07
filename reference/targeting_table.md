# Targeting and reliability summary as a table

The person and item location moments (with and without extreme persons),
the threshold range and coverage, and the reliability indices – PSI, PSI
without extremes, item separation, and coefficient alpha – as a
two-column table suitable for saving and reporting.

## Usage

``` r
targeting_table(fit)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

## Value

A data frame with columns `statistic` and `value`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
targeting_table(rasch(X))
#>                             statistic  value
#> 1                Person location mean  0.015
#> 2                  Person location SD  1.376
#> 3            Person location skewness  -0.12
#> 4            Person location kurtosis   0.23
#> 5  Person location mean (no extremes)  0.052
#> 6    Person location SD (no extremes)  1.113
#> 7    Item location mean (constrained) -0.000
#> 8                    Item location SD  1.476
#> 9              Item location skewness   0.04
#> 10             Item location kurtosis  -1.40
#> 11                  Threshold minimum -1.877
#> 12                  Threshold maximum  1.920
#> 13  Persons below threshold range (%)   13.3
#> 14  Persons above threshold range (%)    2.7
#> 15                                PSI  0.371
#> 16                         Separation   0.77
#> 17                      Person strata    1.4
#> 18               PSI without extremes  0.140
#> 19                 n without extremes    281
#> 20        Item separation reliability  0.991
#> 21                  Coefficient alpha  0.470
#> 22                 n complete (alpha)    300
```
