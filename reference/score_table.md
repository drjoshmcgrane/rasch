# Raw score to measure conversion table

The score-to-logit conversion for complete responders: every possible
raw score with its location, standard error, and the frequency and
cumulative percentage of complete responders at that score (the
complete-data estimates table of Andrich and Marais 2019, ch. 10).

## Usage

``` r
score_table(
  fit,
  method = c("wle", "mle"),
  extremes = c("model", "extrapolated")
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- method:

  `"wle"` (Warm, default) or `"mle"`.

- extremes:

  `"model"` keeps the estimator's own extreme-score values (`NA` for
  MLE); `"extrapolated"` applies the geometric extrapolation.

## Value

A data frame with `score`, `theta`, `se`, `freq`, `cum_pct` (omitted
when no complete responders exist), and `extrapolated`; `NULL` for fits
without a common raw-score metric (EFRM).

## Details

Two estimators are available. `"wle"` (the default) is Warm's weighted
likelihood estimate, finite at the extreme scores. `"mle"` is the plain
maximum likelihood estimate, infinite at the extremes.
`extremes = "extrapolated"` replaces the extreme-score entries by the
geometric extrapolation described in Andrich and Marais (2019, ch. 10):
successive score-to-score differences grow towards the extremes, so the
last difference is continued geometrically – the extrapolated top
difference \\d\\ solves \\b = \sqrt{a d}\\ where \\a, b\\ are the two
preceding differences (equivalently \\d = b^2/a\\), and symmetrically at
zero. The standard error at an extrapolated location is
\\1/\sqrt{I(\theta)}\\ evaluated there. With `method = "wle"` the
extrapolation replaces the finite Warm estimates at the extremes, giving
the extrapolated form of the conversion table from a WLE analysis.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
score_table(rasch(X), method = "mle", extremes = "extrapolated")
#>   score       theta        se extrapolated freq    cum_pct
#> 1     0 -3.27788845 1.8246799         TRUE   13   4.333333
#> 2     1 -1.95393671 1.1597821        FALSE   35  16.000000
#> 3     2 -0.88068333 0.9582186        FALSE   63  37.000000
#> 4     3 -0.01065703 0.9221359        FALSE   73  61.333333
#> 5     4  0.86803622 0.9667528        FALSE   66  83.333333
#> 6     5  1.96255223 1.1702696        FALSE   41  97.000000
#> 7     6  3.32590095 1.8484082         TRUE    9 100.000000
```
