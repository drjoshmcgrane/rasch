# Resolve differential item functioning by iterative item splitting

Splits items with uniform DIF one at a time, beginning with the largest
estimated effect, and refits after each split. This order addresses the
artificial DIF that a large departure can induce in otherwise invariant
items (Andrich and Hagquist 2012, 2015). The DIF in each round is judged
by the residual analysis of
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
or, with `criterion = "wald"`, by the conditional Wald test of
[`dif_wald`](https://drjoshmcgrane.github.io/rasch/reference/dif_wald.md),
which compares each item's split locations on the calibration anchored
by the unsplit items and so needs no class intervals. Each split gives
the item a separate location in every factor cell. A PCM also estimates
the split copies' thresholds separately; an RSM retains its common
rating-scale threshold structure. A location split does not model a
group-specific discrimination, so items with non-uniform DIF are left
for review rather than being made untestable by a split. The procedure
stops when no resolvable uniform DIF remains or the remaining unsplit
reference set reaches `min_anchors`. Items fixed by external anchors are
not split.

## Usage

``` r
resolve_dif(
  fit,
  factors = NULL,
  alpha = 0.05,
  p_adjust = "holm",
  min_n = 20L,
  min_anchors = NULL,
  max_splits = NULL,
  effects = c("main", "factorial"),
  criterion = c("anova", "wald")
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md).

- factors:

  Person factors to test, as in
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md);
  defaults to every nominated factor.

- alpha:

  Significance level for the adjusted probabilities.

- p_adjust:

  Multiplicity adjustment for the DIF tests in each round.

- min_n:

  Minimum distinct responders required in every item-by-factor cell
  before an automatic split is allowed. Repeated response rows from one
  person count once within a cell. The omnibus DIF test determines
  whether a split is needed; pairwise follow-ups describe where the
  difference lies but are not a second significance gate.

- min_anchors:

  Minimum number of original items to leave unsplit as the internal
  reference set. The procedure stops before this set becomes smaller;
  pervasive DIF is not artificial DIF. Default `max(3, items / 4)`.

- max_splits:

  Hard cap on the number of splits. Default: the number of items.

- effects:

  `"main"` fits the factors additively; `"factorial"` also tests their
  interactions. The same model is used at every round and in the final
  DIF assessment.

- criterion:

  `"anova"` flags and ranks items by the residual analysis of
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
  (uniform DIF only is split; a significant non-uniform term leaves the
  item for review). `"wald"` flags items by the conditional Wald test of
  [`dif_wald`](https://drjoshmcgrane.github.io/rasch/reference/dif_wald.md)
  on each factor's main effect and ranks them by the resolved location
  shift; it requires `effects = "main"` and does not test non-uniform
  DIF.

## Value

A list of class `"rasch_resolve_dif"`: the final resolved `fit`, the
`splits` performed (`order`, `item`, `factor`, `base_item`, `eta2`,
`magnitude` in logits), the `stopped` reason, the residual `dif` table,
and the number of distinct source items that still show DIF in the final
fit. `n_untested` counts the uniform and non-uniform hypotheses the
final assessment could not estimate although the design could answer
them; those terms are reported as neither DIF nor no DIF, so the
remaining-DIF count is a lower bound whenever `n_untested` is positive.
A split copy answered in one level of its splitting factor only is not
counted: its term is structurally absent, not lost. A split's
`magnitude` is `NA` when the complete-design post-hoc comparison that
supplies it failed or returned no trustworthy finite estimate; the split
still stands, and `notes` names the item, the factor and the reason. An
`NA` magnitude is an unmeasured difference, not a zero one.
`n_remaining_dif` is `NA` when no hypothesis was estimable.
`n_nonuniform` counts significant non-uniform item-factor findings and
is `NA` if any answerable non-uniform hypothesis is unavailable, or no
hypothesis was estimable. `n_untested` is always a count. `effects`
records the factor model used and `criterion` the test. With the Wald
criterion `eta2` is `NA`, `magnitude` is the absolute resolved shift the
test compared (the range of the locations for more than two levels), the
`dif` table carries the `shift` and `p_adj` of each item still flagged,
and `n_nonuniform` is `NA` because non-uniform DIF is not tested.

## References

Andrich, D., & Hagquist, C. (2012). Real and artificial differential
item functioning. *Journal of Educational and Behavioral Statistics*,
37(3), 387-416.

## See also

[`split_items`](https://drjoshmcgrane.github.io/rasch/reference/split_items.md)
for a single split,
[`drop_items`](https://drjoshmcgrane.github.io/rasch/reference/drop_items.md)
to remove an item instead,
[`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md)
and
[`dif_wald`](https://drjoshmcgrane.github.io/rasch/reference/dif_wald.md)
for the tests it resolves.

## Examples

``` r
set.seed(1); n <- 300
d <- seq(-2, 2, length.out = 6); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 6); sh[g == "b", 3] <- 1.5      # one strong DIF item
X <- matrix(rbinom(n * 6, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 6)
colnames(X) <- paste0("I", 1:6)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
resolve_dif(fit)$splits
#>  order item factor base_item  eta2 magnitude
#>      1   I3    grp        I3 0.108     1.726
resolve_dif(fit, criterion = "wald")$splits
#>  order item factor base_item eta2 magnitude
#>      1   I3    grp        I3          1.726
```
