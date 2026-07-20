# Resolve differential item functioning by iterative item splitting

Splits DIF items one at a time, largest effect first, refitting after
each split, until no item shows significant DIF (or the anchor set would
fall too low). Splitting the item with the largest real DIF first
removes the artificial DIF it induces on other items (Andrich & Hagquist
2012, 2015), so the procedure resolves genuine DIF without chasing the
artificial DIF that a single simultaneous pass would flag. Each split
resolves the item into one copy per group (or per factor-combination
cell for an interaction), with independent locations and thresholds, so
both uniform and non-uniform DIF are resolved together.

## Usage

``` r
resolve_dif(
  fit,
  factors = NULL,
  alpha = 0.05,
  p_adjust = "BH",
  min_anchors = NULL,
  max_splits = NULL
)
```

## Arguments

- fit:

  A fitted object from
  [`rasch`](https://drjoshmcgrane.github.io/rasch/reference/rasch.md)
  carrying person factors.

- factors:

  Person factors to test, as in
  [`dif_anova`](https://drjoshmcgrane.github.io/rasch/reference/dif_anova.md);
  defaults to every nominated factor.

- alpha:

  Significance level for the adjusted probabilities.

- p_adjust:

  Multiplicity adjustment across items each round.

- min_anchors:

  Minimum number of original items to leave unsplit; the procedure stops
  before the anchor set falls below this (pervasive DIF is not
  artificial DIF). Default `max(3, items / 4)`.

- max_splits:

  Hard cap on the number of splits. Default: the number of items.

## Value

A list of class `"rasch_resolve_dif"`: the final resolved `fit`, the
`splits` performed (order, item, factor, partial eta-squared, DIF
magnitude in logits), the `stopped` reason, and the residual `dif` table
for the final fit.

## References

Andrich, D., & Hagquist, C. (2012). Real and artificial differential
item functioning. *Journal of Educational and Behavioral Statistics*,
37(3), 387-416.

## Examples

``` r
set.seed(1); n <- 600
d <- seq(-2, 2, length.out = 8); g <- rep(c("a", "b"), each = n / 2)
sh <- matrix(0, n, 8); sh[g == "b", 3] <- 1.2      # one strong DIF item
X <- matrix(rbinom(n * 8, 1, plogis(outer(rnorm(n), d, "-") - sh)), n, 8)
colnames(X) <- paste0("I", 1:8)
fit <- rasch(data.frame(X, grp = g), factors = "grp")
resolve_dif(fit)$splits
#>   order item factor       eta2 magnitude
#> 1     1   I3    grp 0.05389885  1.110898
```
