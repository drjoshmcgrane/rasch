# Propose polytomous option scores from the distractor evidence

Multiple-choice items can be rescored polytomously so that a distractor
carrying information about the trait receives partial credit (Andrich
and Styles 2011). This function proposes such a scoring from the
rest-measure distractor analysis: within each keyed item, a distractor
qualifies for credit when it attracts at least `min_n` takers, its
takers' mean rest location exceeds that of the uncredited distractors by
more than `z` standard errors of the difference, and it remains below
the keyed option. Qualifying distractors are ranked by mean location and
scored 1, 2, ... below the keyed option's top score. The result is a
proposal for substantive review, not an automatic decision: inspect
[`plot_distractors`](https://drjoshmcgrane.github.io/rmt/reference/plot_distractors.md)
and the item content, edit as needed, then refit with
`rasch(raw_data, key = proposal$option_scores)`.

## Usage

``` r
distractor_rescore(fit, items = NULL, min_n = 20, z = 1.96)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md) run
  with a `key`.

- items:

  Optional subset of keyed item names.

- min_n:

  Minimum takers for a distractor to be considered.

- z:

  Required separation, in standard errors, between a credited distractor
  and the uncredited ones.

## Value

A list of class `"rmt_rescore"`: `option_scores`, a data frame (`item`,
`option`, `score`) ready for `rasch(key = )` and covering every observed
option of the examined items, and `evidence`, the distractor analysis
with the proposed scores and the separation z per option.

## References

Andrich, D. and Styles, I. (2011). Distractors with information in
multiple choice items: A rationale based on the Rasch model. Journal of
Applied Measurement, 12, 67-95.

## Examples

``` r
set.seed(1); Np <- 600
th <- rnorm(Np)
raw <- sapply(seq(-0.5, 0.5, length.out = 4), function(d) {
  x <- vapply(th, function(b) sample(0:2, 1,
    prob = item_moments(b, c(d - 0.7, d + 0.7))$P), 0L)
  c("D", "B", "A")[x + 1]   # B is an informative distractor
})
colnames(raw) <- paste0("M", 1:4)
fit <- rasch(raw, key = setNames(rep("A", 4), colnames(raw)))
pr <- distractor_rescore(fit)
pr$option_scores
#>    item option score
#> 1    M1      A     2
#> 2    M1      B     1
#> 3    M1      D     0
#> 4    M2      A     2
#> 5    M2      B     1
#> 6    M2      D     0
#> 7    M3      A     2
#> 8    M3      B     1
#> 9    M3      D     0
#> 10   M4      A     2
#> 11   M4      B     1
#> 12   M4      D     0
```
