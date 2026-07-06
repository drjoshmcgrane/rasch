# Warm's weighted likelihood estimates by raw score

Computes the weighted likelihood estimate (WLE) of person location for
every possible raw score on a set of items, with standard errors. WLE
estimates are finite at the extreme (zero and maximum) scores, unlike
the maximum likelihood estimate.

## Usage

``` r
person_wle(tau_list, disc = 1)
```

## Arguments

- tau_list:

  List of per-item threshold vectors.

- disc:

  Common discrimination (frame unit) of the items; with a constant
  discrimination the raw score remains sufficient.

## Value

A list with `theta` and `se`, each named by raw score.

## Examples

``` r
person_wle(list(c(-1, 0), c(-0.5, 0.5), c(0, 1)))
#> $theta
#>             0             1             2             3             4 
#> -2.444406e+00 -1.237127e+00 -5.611339e-01  3.552714e-15  5.611339e-01 
#>             5             6 
#>  1.237127e+00  2.444406e+00 
#> 
#> $se
#>         0         1         2         3         4         5         6 
#> 1.5455221 0.9731988 0.8347598 0.8003924 0.8347598 0.9731988 1.5455221 
#> 
```
