# Planned DIF contrasts derived from the factor structure

The confirmatory alternative to exhaustive post-hoc comparison: instead
of every pair of design cells, a small family of one-degree-of-freedom
questions is tested, so familywise control costs little power (Maxwell
and Delaney 2004, ch. 5). By default the family is derived from the
structure of the factors themselves – a two-level factor contributes its
difference; an ordered factor (declared ordered, or with numeric levels
such as ages or waves) contributes its linear and quadratic trends; a
nominal factor contributes all pairs when it has up to four levels and
each-level-against-the-rest otherwise; and every pair of factors with a
leading contrast (a difference or a linear trend) contributes the
product interaction. Print the returned object to see the family in
words before reading the results; a family endorsed in advance of the
results is what makes the contrasts planned.

## Usage

``` r
dif_contrasts(
  fit,
  factors = NULL,
  items = NULL,
  within = NULL,
  id = NULL,
  contrasts = "auto",
  p_adjust = "holm",
  alpha = 0.05,
  flag_logits = 0.5,
  min_n = 20
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rmt/reference/rasch.md).

- factors:

  A data frame of person factors, a character vector naming factors
  nominated in the fit, or a single grouping vector. Defaults to every
  factor stored in the fit.

- items:

  Item names or indices to test; all items by default.

- within:

  Names of factors that vary within person (for example time). Detected
  automatically when `id` is supplied and a factor varies within an id.

- id:

  Person identifier with one entry per row, or the name of a nominated
  factor holding it; required for stacked designs where the same person
  occupies several rows.

- contrasts:

  `"auto"` (derive the family from the factor structure) or a named list
  of numeric cell-weight vectors, each named by the design-cell labels
  (factor levels joined by `":"`). Weights are rescaled so the positive
  and negative parts each sum to one.

- p_adjust:

  Familywise adjustment over the whole family (items by contrasts);
  default `"holm"`.

- alpha:

  Significance level for the adjusted probabilities.

- flag_logits:

  Absolute estimate flagged as practically significant.

- min_n:

  Cells with fewer responders to an item are dropped from that item's
  resolution, with a note.

## Value

A list of class `"rmt_dif_contrasts"`: `table` (one row per item and
contrast: estimate in logits, SE, statistic, df where a t test was used,
raw and adjusted p, 95 per cent interval, `significant`, `practical`,
`within`), `family` (the derived questions with their cell weights), the
settings, and any `notes`.

## Details

Each contrast is estimated in logits from resolved item locations (the
item split into one copy per design cell and the model refitted, as in
[`dif_size`](https://drjoshmcgrane.github.io/rmt/reference/dif_size.md)),
with cell weights scaled so every estimate is a difference between two
weighted averages – directly comparable to the practical-significance
criterion. Because resolution is used, magnitudes are read from a
calibration in which compensating artificial DIF has been removed
(Andrich and Hagquist 2015).

When `id` shows that persons repeat across rows (a stacked
repeated-measures design), between-row independence fails and the usual
tests would be invalid. Significance is then computed from person-level
scores of the standardised residuals: a within-subject contrast (for
example a trend over time) becomes one contrast score per person, tested
against zero; a between-subjects contrast is tested on person-mean
residuals; and a between-by-within interaction tests the person contrast
scores across the between groups. Logit estimates are still reported
from the resolved locations; their standard errors treat rows as
independent and are conservative for within-subject differences.

## References

Maxwell, S. E., & Delaney, H. D. (2004). *Designing Experiments and
Analyzing Data* (2nd ed.). Mahwah, NJ: Erlbaum.

Andrich, D., & Hagquist, C. (2015). Real and artificial differential
item functioning in polytomous items. *Educational and Psychological
Measurement*, 75(2), 185-207.

Hagquist, C., & Andrich, D. (2017). Recent advances in analysis of
differential item functioning in health research using the Rasch model.
*Health and Quality of Life Outcomes*, 15, 181.

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
dif_contrasts(fit, items = c("I3", "I5"))
#> Planned DIF contrasts (1 questions x 2 items; holm over the family)
#>   grp: b - a
#> 
#>  item   contrast estimate    se statistic   p_adj significant practical
#>    I3 grp: b - a    0.907 0.204     4.441 < 0.001           *         *
#>    I5 grp: b - a   -0.570 0.200    -2.846   0.004           *         *
```
