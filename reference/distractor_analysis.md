# Distractor analysis for multiple-choice items

For every keyed item and response option: the count and proportion
choosing it, the mean location of those persons, and the point-biserial
correlation between choosing the option and the person measure.
Locations and correlations use the rest measure (the person estimate
from the other items), so the analysed item cannot credit its own
takers. The keyed option should attract the ablest persons and carry the
only positive point-biserial; a distractor whose takers are abler than
the keyed option's (with at least `min_n` takers) is flagged as a
possible miskey.

## Usage

``` r
distractor_analysis(fit, items = NULL, min_n = 10)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  run with a `key`.

- items:

  Optional subset of item names; defaults to every keyed item.

- min_n:

  Minimum takers for an option to be eligible for the miskey flag.

## Value

A data frame with one row per item-option: `item`, `option`, its
assigned `score`, `keyed` (full credit), `n`, `prop`, `mean_location`,
`point_biserial`, and `flag`.

## Examples

``` r
set.seed(1); Np <- 400
th <- rnorm(Np)
raw <- sapply(seq(-1, 1, length.out = 6), function(d) {
  ok <- rbinom(Np, 1, plogis(th - d))
  ifelse(ok == 1, "A", sample(c("B", "C", "D"), Np, replace = TRUE))
})
colnames(raw) <- paste0("M", 1:6)
fit <- rasch(raw, key = setNames(rep("A", 6), colnames(raw)))
head(distractor_analysis(fit))
#>   item option score keyed   n   prop mean_location point_biserial  flag
#> 1   M1      A     1  TRUE 277 0.6925     0.2457337     0.23025114 FALSE
#> 2   M1      B     0 FALSE  40 0.1000    -0.3836142    -0.12294026 FALSE
#> 3   M1      C     0 FALSE  38 0.0950    -0.2619386    -0.08678212 FALSE
#> 4   M1      D     0 FALSE  45 0.1125    -0.4096550    -0.13900664 FALSE
#> 5   M2      A     1  TRUE 245 0.6125     0.3067212     0.26204987 FALSE
#> 6   M2      B     0 FALSE  60 0.1500    -0.3834449    -0.15309278 FALSE
```
