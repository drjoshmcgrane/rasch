# Equate two test calibrations through their common items

Compares the item locations of a fit with those of a second fit (or a
reference table such as an item bank), matched by item name. A scale
shift between the two origins is estimated by the precision-weighted
mean difference, and each common item is then tested against the shifted
identity line; flagged items show drift and weaken the equating link.

## Usage

``` r
equate_tests(fit, reference, shift = c("mean", "none"))
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- reference:

  A second
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  fit, or a data frame with columns `item`, `location`, and optionally
  `se`.

- shift:

  `"mean"` (default) allows a scale shift between the two analyses;
  `"none"` compares raw locations, appropriate when both analyses are
  already on a shared (anchored) scale.

## Value

A list with the comparison `table` (locations, standard errors,
difference, t, raw and BH-adjusted p, drift flag), the estimated
`shift`, the location `correlation`, the root mean square difference
after shifting (`rmsd`), and the number of common items `n`.

## Examples

``` r
set.seed(1); d <- seq(-1.5, 1.5, length.out = 8)
mk <- function() {
  X <- matrix(rbinom(400 * 8, 1, plogis(outer(rnorm(400), d, "-"))), 400, 8)
  colnames(X) <- paste0("I", 1:8); rasch(X)
}
eq <- equate_tests(mk(), mk())
eq$table
#>   item location_1      se_1 location_2      se_2  difference adj_difference
#> 1   I1 -1.4780195 0.1246379 -1.5760281 0.1292289  0.09800863     0.10727215
#> 2   I2 -1.1006867 0.1168127 -1.1398494 0.1155055  0.03916271     0.04842624
#> 3   I3 -0.8227470 0.1123970 -0.5314218 0.1083620 -0.29132523    -0.28206171
#> 4   I4 -0.2085726 0.1093047 -0.2287656 0.1063021  0.02019291     0.02945643
#> 5   I5  0.1917743 0.1077741  0.3229450 0.1053953 -0.13117072    -0.12190720
#> 6   I6  0.7088231 0.1124051  0.5366178 0.1097245  0.17220529     0.18146881
#> 7   I7  1.1279739 0.1189861  1.1874751 0.1180524 -0.05950128    -0.05023776
#> 8   I8  1.5814546 0.1284329  1.4290269 0.1263670  0.15242768     0.16169120
#>            t          p     p_adj drift
#> 1  0.5974823 0.55018539 0.8468081 FALSE
#> 2  0.2947852 0.76815795 0.8468081 FALSE
#> 3 -1.8066264 0.07082055 0.5665644 FALSE
#> 4  0.1931926 0.84680810 0.8468081 FALSE
#> 5 -0.8087095 0.41868227 0.8373645 FALSE
#> 6  1.1552580 0.24798478 0.8373645 FALSE
#> 7 -0.2997250 0.76438690 0.8468081 FALSE
#> 8  0.8974040 0.36950338 0.8373645 FALSE
```
