# Enumerate item-category thresholds

Builds the index mapping each item-category threshold to a global id,
given the maximum score of each item.

## Usage

``` r
threshold_index(m)
```

## Arguments

- m:

  Integer vector of maximum scores per item (1 for dichotomous items).

## Value

A data frame with columns `id`, `item`, and `k` (the within-item
threshold number).

## Examples

``` r
threshold_index(c(1, 3, 2))
#>   id item k
#> 1  1    1 1
#> 2  2    2 1
#> 3  3    2 2
#> 4  4    2 3
#> 5  5    3 1
#> 6  6    3 2
```
