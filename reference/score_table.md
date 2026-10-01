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

  Treatment of the extreme scores. `"model"` keeps the estimator's own
  values; these are `NA` for MLE. `"extrapolated"` applies the geometric
  extrapolation.

## Value

A data frame with `score`, `theta`, `se`, `freq`, `cum_pct` (omitted
when no complete responders exist), and `extrapolated`, preceded by
`form` when the fit has split items; `NULL` when the fitted items do not
share one discrimination or an item is represented by several MFRM or
EFRM response cells.

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

A fit with split items
([`split_items`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md),
[`resolve_dif`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md))
has no common raw score: each person answers one copy of a split item,
so a score summed over every calibrated column belongs to nobody. The
table then holds one conversion per test form, the copies that a set of
persons answered together with every unsplit item, under a leading
`form` column that names those copies. A person who answered no copy of
a split item, because the copy was dropped, the person's level was not
split or the response is missing, takes the form without that item,
named so. Each form's frequencies count its own complete responders.

## References

Andrich, D. and Marais, I. (2019). A Course in Rasch Measurement Theory:
Measuring in the Educational, Social and Health Sciences. Springer.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. Psychometrika, 54(3), 427–450.

## Examples

``` r
set.seed(1)
d <- seq(-1.5, 1.5, length.out = 6)
X <- matrix(rbinom(300 * 6, 1, plogis(outer(rnorm(300), d, "-"))), 300, 6)
colnames(X) <- paste0("I", 1:6)
score_table(rasch(X), method = "mle", extremes = "extrapolated")
#>  score  theta    se extrapolated freq cum_pct
#>      0 -3.278 1.825            *   13   4.333
#>      1 -1.954 1.160                35  16.000
#>      2 -0.881 0.958                63  37.000
#>      3 -0.011 0.922                73  61.333
#>      4  0.868 0.967                66  83.333
#>      5  1.963 1.170                41  97.000
#>      6  3.326 1.848            *    9 100.000
```
