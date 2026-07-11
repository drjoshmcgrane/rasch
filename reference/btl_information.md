# Information and targeting of a paired-comparison design

The paired-comparison analogue of the test-information function. The
Fisher information a single comparison carries about the location
difference `d = beta_a - beta_b` is, in this exponential family, the
variance of its score – `P(1 - P)` for the dichotomous choice and the
graded response variance `V` for the ordinal extension (the score is the
sufficient statistic for `d`, so its variance is the information).
Weighted by each comparison's replication count and summed over the
comparisons the design actually contains, this gives a *design
information* for every object: how much the observed comparisons pin its
location down, the counterpart of an item's contribution to test
information. Because the information peaks at gap zero and falls away
with the location gap, near-neighbour contests are the informative ones.

## Usage

``` r
btl_information(fit)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

## Value

A list of class `"rasch_btl_info"`: `objects` (per object: `location`,
the fit's `se`, `n_comparisons`, the design `information`, and
`se_naive`); `pairs` (per observed pair: `n`, the mean location `gap`,
and the pair's `information`); `comparisons` (per comparison: the signed
`gap`, `weight`, and the single-comparison `information`); the scalar
`total` information; `m`; the `clustered` flag; and `notes`.

## Details

The design information inverts to `se_naive = 1 / sqrt(information)`:
the error the object's comparisons would give if its location were the
ONLY free parameter – a single-parameter lower bound, useful for reading
which objects the design serves well. It is not the model's standard
error even with independent comparisons (every location is estimated
jointly with the others, and the fit's own `se` is additionally the
judge-clustered Godambe sandwich), so `se` sits above `se_naive` as a
rule; treat their ratio as descriptive, not as a clustering test.

## References

Pollitt, A. (2012). The method of adaptive comparative judgement.
*Assessment in Education*, 19(3), 281-300.

## See also

[`plot_btl_targeting`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md),
[`btl_next_pairs`](https://drjoshmcgrane.github.io/rasch/reference/btl_next_pairs.md)

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_information(btl(d, "a", "b", "win"))
#> Paired-comparison design information: 4 objects, total 30.04
#> One-comparison Fisher information about the location gap (dichotomous: P(1 - P))
#>  object location    se n_comparisons information se_naive
#>       A   -1.238 0.214            90      12.487    0.283
#>       B   -0.354 0.186            90      17.100    0.242
#>       C    0.448 0.180            90      17.030    0.242
#>       D    1.144 0.209            90      13.463    0.273
#> Note: se is the Godambe sandwich standard error; se_naive = 1/sqrt(information) is a single-parameter lower bound (as if the object's location were the only free parameter), so se sits above it as a rule
```
