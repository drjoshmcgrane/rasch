# Recommend the next informative comparisons (adaptive step)

The adaptive comparative judgement step of Pollitt (2012): rank
candidate object pairs by the information one additional comparison
would carry at the current estimates. That information peaks when the
two objects are close in location, so at equal measurement the
recommender favours near-neighbour contests. With `weight_se = TRUE`
(the default) each pair's priority is the one-step reduction in TOTAL
location variance that one added comparison of the pair would deliver,
from a rank-one (Sherman-Morrison) update of the fit's stored covariance
with the comparison's information on the contrast, so pairs of poorly
measured (and correlated) objects are promoted. The update formula is
exact for a model-based information matrix; applied to the sandwich
covariance it is a scoring device, consistent with the ranking-heuristic
status described below.

## Usage

``` r
btl_next_pairs(fit, n = 10, weight_se = TRUE)
```

## Arguments

- fit:

  A paired-comparison fit from
  [`btl`](https://drjoshmcgrane.github.io/rasch/reference/btl.md).

- n:

  Number of pairs to return. The priority is a greedy one-step RANKING
  heuristic: it plugs the judge-clustered sandwich covariance into an
  information-update formula that is exact only for a model-based
  information matrix, so the ranking orders candidate pairs sensibly but
  the implied variance reductions are not exact sandwich updates. Treat
  the ordering, not the magnitudes, as the output.

- weight_se:

  Rank by the one-step total-variance reduction (default `TRUE`; falls
  back to `expected_information * (se_a^2 + se_b^2)` if the fit carries
  no covariance). When `FALSE`, pairs are ranked by expected information
  alone (pure closeness).

## Value

A data frame of the top `n` candidate pairs, each oriented to its
stronger object: `object_a`, `object_b`, the location `gap`,
`n_existing` (replications already observed for the pair),
`expected_information` (of one new comparison), and `priority`. Sorted
by `priority` (or by `expected_information` when `weight_se = FALSE`).

## Details

Two honest cautions. This is a *greedy* rule that scores each pair on
its own immediate one-step gain (an A-optimality step at the current
estimates, taking the clustered covariance as the state); it is not a
full optimal design and can be beaten by one that plans several
comparisons jointly. And adaptive selection is known to inflate a
separation (scale) reliability computed naively afterwards, because the
design concentrates comparisons where they shrink the errors most:
report reliability from an independent or non-adaptive subset, or treat
an adaptive reliability as an upper bound (Bramley 2015).

## References

Pollitt, A. (2012). The method of adaptive comparative judgement.
*Assessment in Education*, 19(3), 281-300. Bramley, T. (2015).
Investigating the reliability of Adaptive Comparative Judgment.
*Cambridge Assessment Research Report*.

## See also

[`btl_information`](https://drjoshmcgrane.github.io/rasch/reference/btl_information.md),
[`plot_btl_targeting`](https://drjoshmcgrane.github.io/rasch/reference/plot_btl_targeting.md)

## Examples

``` r
set.seed(1)
beta <- c(A = -1, B = -0.3, C = 0.4, D = 0.9)
pr <- t(combn(names(beta), 2))
d <- data.frame(a = rep(pr[, 1], each = 30), b = rep(pr[, 2], each = 30))
d$win <- ifelse(runif(nrow(d)) < plogis(beta[d$a] - beta[d$b]), d$a, d$b)
btl_next_pairs(btl(d, "a", "b", "win"), n = 5)
#>   object_a object_b       gap n_existing expected_information     priority
#> 1        B        A 0.8837971         30            0.2068987 0.0010120610
#> 2        D        C 0.6961084         30            0.2220026 0.0009330875
#> 3        D        B 1.4978719         30            0.1493481 0.0009279925
#> 4        C        A 1.6855607         30            0.1319119 0.0008269402
#> 5        C        B 0.8017636         30            0.2137663 0.0007881426
```
