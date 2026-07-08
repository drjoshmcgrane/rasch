# Guttman-ordered response matrix and reproducibility

Orders persons by location (descending) and items by location
(ascending), the Guttman scalogram arrangement, and computes the
coefficient of reproducibility against the deterministic Guttman pattern
implied by each person's total score.

## Usage

``` r
guttman_table(fit)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

## Value

A list with the ordered score matrix `matrix` (persons by items, row and
column names carrying the ID and item labels), the person and item
orderings, and the coefficient of reproducibility `CR`.

## Examples

``` r
set.seed(1)
d <- seq(-2, 2, length.out = 6)
X <- matrix(rbinom(200 * 6, 1, plogis(outer(rnorm(200), d, "-"))), 200, 6)
colnames(X) <- paste0("I", 1:6)
guttman_table(rasch(X))$CR
#> [1] 0.8116667
```
