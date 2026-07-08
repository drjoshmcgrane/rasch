# Test information function

Fisher information of the whole test over a grid of person locations,
with the corresponding standard error of measurement.

## Usage

``` r
test_information(fit, grid = seq(-6, 6, by = 0.1))
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- grid:

  Logit grid over which to evaluate the information.

## Value

A data frame with `theta`, `info`, and `sem`.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
head(test_information(rasch(X)))
#>   theta       info      sem
#> 1  -6.0 0.02364631 6.503068
#> 2  -5.9 0.02609577 6.190346
#> 3  -5.8 0.02879468 5.893101
#> 4  -5.7 0.03176749 5.610590
#> 5  -5.6 0.03504089 5.342105
#> 6  -5.5 0.03864388 5.086976
```
