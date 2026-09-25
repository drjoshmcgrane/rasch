# Conditional Wald test of DIF on the resolved calibration

Tests each item for uniform DIF by splitting it by a person factor,
recalibrating with the unsplit items as the anchor, and comparing the
locations of the split copies with a Wald test. The comparison
conditions on item-pair totals through the pairwise conditional
likelihood used by
[`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
with its Godambe sandwich covariance. No class intervals or person
estimates are needed. It addresses item-location invariance, as
Andersen's (1973) test does at the test level, but is not a full-score
conditional likelihood-ratio test. For polytomous items it compares mean
threshold locations, not every threshold separately.

## Usage

``` r
dif_wald(
  fit,
  factors = NULL,
  items = NULL,
  p_adjust = "holm",
  alpha = 0.05,
  min_n = 20L
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md),
  including one that already carries splits.

- factors:

  Person factors to test, as in
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md);
  defaults to every nominated factor. Each factor is tested on its own.

- items:

  Items to test; defaults to every item in the fit.

- p_adjust:

  Multiplicity adjustment over every item-by-factor test, as in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- alpha:

  Significance level for the adjusted probabilities.

- min_n:

  Levels with fewer distinct responders to an item are dropped from that
  item's comparison.

## Value

A list of class `"rasch_dif_wald"`. `summary` has one row per item and
factor: `n_levels` compared, `shift` (for two levels the location of the
second level minus the first, positive when the item is harder for the
second level; for more levels the range of the locations), its `se` (two
levels only), the `wald` statistic, its hypothesis `df` (levels minus
one), the reference `ref_df` (infinite for independent response rows,
the cluster count minus one for a supported repeated-person
calibration), `p`, `p_adj` and `significant`. `levels` gives each
compared level's resolved `location`, `se` and `n`. `notes` records what
was left out or withheld and why.

## Details

A split copy of an item that persons in one level of the factor only
answered has no group contrast and is left out with a note. Levels with
fewer than `min_n` distinct responders to an item are dropped from that
item's comparison. Items whose resolution is unavailable (an externally
anchored item, a group that does not observe every score category, a
refit that fails, a location that rests on a near-empty category) keep
their row with `NA` statistics and a note.

Every other item is treated as invariant while one item is tested. Under
pervasive DIF that assumption fails and the anchor carries artificial
DIF (Andrich and Hagquist 2012), so the test is best used iteratively,
as
[`resolve_dif`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
does with `criterion = "wald"`, or with a set of anchors chosen on other
grounds.

## References

Andersen, E. B. (1973). A goodness of fit test for the Rasch model.
Psychometrika, 38(1), 123–140.

Glas, C. A. W. and Verhelst, N. D. (1995). Testing the Rasch model. In
G. H. Fischer and I. W. Molenaar (eds), Rasch Models: Foundations,
Recent Developments, and Applications (pp. 69–95). Springer.

Kopf, J., Zeileis, A. and Strobl, C. (2015). Anchor selection strategies
for DIF analysis: review, assessment, and new approaches. Educational
and Psychological Measurement, 75(1), 22–56.

Andrich, D. and Hagquist, C. (2012). Real and artificial differential
item functioning. Journal of Educational and Behavioral Statistics,
37(3), 387–416.

## See also

[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
for the residual analysis by class interval,
[`dif_size`](https://drjoshmcgrane.github.io/rasch/reference/dif_size.md)
for the pairwise magnitudes of one item and
[`resolve_dif`](https://drjoshmcgrane.github.io/rasch/reference/resolve_dif.md)
for iterative resolution.

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 0.8
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
dif_wald(fit)
#> Conditional Wald test of DIF by grp (resolved locations, logits)
#>  item factor n_levels  shift    se   wald df ref_df       p   p_adj significant
#>    I1    grp        2 -0.087 0.265  0.108  1    Inf   0.742   1.000            
#>    I2    grp        2 -0.198 0.237  0.697  1    Inf   0.404   1.000            
#>    I3    grp        2  0.907 0.204 19.722  1    Inf < 0.001 < 0.001           *
#>    I4    grp        2 -0.133 0.200  0.442  1    Inf   0.506   1.000            
#>    I5    grp        2 -0.570 0.200  8.102  1    Inf   0.004   0.031           *
#>    I6    grp        2 -0.261 0.208  1.581  1    Inf   0.209   1.000            
#>    I7    grp        2  0.361 0.230  2.459  1    Inf   0.117   0.701            
#>    I8    grp        2  0.065 0.270  0.058  1    Inf   0.809   1.000            
#> p adjusted by holm over 8 item-by-factor test(s); chi-square reference
#> shift: second level minus first for two levels, range of the locations otherwise; positive means harder for the second level
```
